# ----------------------------------------------------------------------------
#

# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import hashlib
import subprocess

import biom
import skbio
import qiime2.util
import pandas as pd

from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import (
    FastqGzFormat, SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)


def run_commands(cmds, verbose=True):
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


def _check_featureless_table(fp):
    with open(fp) as fh:
        # There is a header before the feature data
        for line_count, _ in zip(range(1, 3), fh):
            pass
    if line_count < 2:
        raise ValueError("No features remain after denoising. Try adjusting "
                         "your truncation and trim parameter settings.")


_WHOLE_NUM = (lambda x: x >= 0, 'non-negative')
_NAT_NUM = (lambda x: x > 0, 'greater than zero')
_POOL_STR = (lambda x: x in {'pseudo', 'independent'},
             'pseudo or independent')
_CHIM_STR = (lambda x: x in {'pooled', 'consensus', 'none'},
             'pooled, consensus or none')
_BOOLEAN = (lambda x: type(x) is bool, 'True or False')
# Better to choose to skip, than to implicitly ignore things that KeyError
_SKIP = (lambda x: True, '')
_valid_inputs = {
    'trunc_len': _WHOLE_NUM,
    'trim_left': _WHOLE_NUM,
    'max_ee': _NAT_NUM,
    'trunc_q': _WHOLE_NUM,
    'min_overlap': _WHOLE_NUM,
    'pooling_method': _POOL_STR,
    'chimera_method': _CHIM_STR,
    'min_fold_parent_over_abundance': _NAT_NUM,
    'allow_one_off': _BOOLEAN,
    'n_threads': _WHOLE_NUM,
    # 0 is technically allowed, but we don't want to support it because it only
    # takes all reads from the first sample (alphabetically by sample id)
    'n_reads_learn': _NAT_NUM,
    # Skipped because they are valid for whole domain of type
    'hashed_feature_ids': _SKIP,
    'demultiplexed_seqs': _SKIP,
    'homopolymer_gap_penalty': _SKIP,
    'band_size': _SKIP,
    'front': _SKIP,
    'adapter': _SKIP,
    'indels': _SKIP,
}


# TODO: Replace this with Range predicates when interfaces support them better
def _check_inputs(**kwargs):
    for param, arg in kwargs.items():
        check_is_valid, explanation = _valid_inputs[param]
        if not check_is_valid(arg):
            raise ValueError('Argument to %r was %r, should be %s.'
                             % (param, arg, explanation))


def _filepath_to_sample(fp):
    return fp.rsplit('_', 4)[0]


# Since `denoise-single` and `denoise-pyro` are almost identical, break out
# the bulk of the functionality to this helper util. Typechecking is assumed
# to have occurred in the calling functions, this is primarily for making
# sure that DADA2 is able to do what it needs to do.
def _denoise_helper(biom_fp, track_fp, hashed_feature_ids):
    _check_featureless_table(biom_fp)
    with open(biom_fp) as fh:
        table = biom.Table.from_tsv(fh, None, None, None)

    df = pd.read_csv(track_fp, sep='\t', index_col=0)
    df.index.name = 'sample-id'
    df = df.rename(index=_filepath_to_sample)

    PASSED_FILTER = 'percentage of input passed filter'
    NON_CHIMERIC = 'percentage of input non-chimeric'

    round_cols = {PASSED_FILTER: 2, NON_CHIMERIC: 2}

    df[PASSED_FILTER] = df['filtered'] / df['input'] * 100
    df[NON_CHIMERIC] = df['non-chimeric'] / df['input'] * 100

    col_order = ['input', 'filtered', PASSED_FILTER, 'denoised',
                 'non-chimeric', NON_CHIMERIC]

    # only calculate percentage of input merged if paired end
    if 'merged' in df:
        MERGED = 'percentage of input merged'
        round_cols[MERGED] = 2
        df[MERGED] = df['merged'] / df['input'] * 100
        col_order.insert(4, 'merged')
        col_order.insert(5, MERGED)

    # only calculate percentage of input primer-removed if ccs
    if 'primer-removed' in df:
        PASSED_PRIMERREMOVE = 'percentage of input primer-removed'
        round_cols[PASSED_PRIMERREMOVE] = 2
        df[PASSED_PRIMERREMOVE] = df['primer-removed'] / df['input'] * 100
        col_order.insert(1, 'primer-removed')
        col_order.insert(2, PASSED_PRIMERREMOVE)

    df = df[col_order]
    df.fillna(0, inplace=True)
    df = df.round(round_cols)
    metadata = qiime2.Metadata(df)

    # Currently the sample IDs in DADA2 are the file names. We make
    # them the sample id part of the filename here.
    sid_map = {id_: _filepath_to_sample(id_)
               for id_ in table.ids(axis='sample')}
    table.update_ids(sid_map, axis='sample', inplace=True)
    # The feature IDs in DADA2 are the sequences themselves.
    if hashed_feature_ids:
        # Make feature IDs the md5 sums of the sequences.
        fid_map = {id_: hashlib.md5(id_.encode('utf-8')).hexdigest()
                   for id_ in table.ids(axis='observation')}
        table.update_ids(fid_map, axis='observation', inplace=True)

        rep_sequences = DNAIterator((skbio.DNA(k, metadata={'id': v})
                                     for k, v in fid_map.items()))
    else:
        rep_sequences = DNAIterator(
            (skbio.DNA(id_, metadata={'id': id_})
             for id_ in table.ids(axis='observation')))
    return table, rep_sequences, metadata


def _denoise_single(demultiplexed_seqs, trunc_len, trim_left, max_ee, trunc_q,
                     pooling_method, chimera_method,
                    min_fold_parent_over_abundance, allow_one_off,
                    n_threads, n_reads_learn, hashed_feature_ids,
                    homopolymer_gap_penalty, band_size):
    _check_inputs(**locals())
    if trunc_len != 0 and trim_left >= trunc_len:
        raise ValueError("trim_left (%r) must be smaller than trunc_len (%r)"
                         % (trim_left, trunc_len))

    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        track_fp = os.path.join(temp_dir_name, 'track.tsv')

        cmd = ['run_dada.R',
               '--input_directory', str(demultiplexed_seqs),
               '--output_path', biom_fp,
               '--output_track', track_fp,
               '--filtered_directory', temp_dir_name,
               '--truncation_length', str(trunc_len),
               '--trim_left', str(trim_left),
               '--max_expected_errors', str(max_ee),
               '--truncation_quality_score', str(trunc_q),
               '--pooling_method', str(pooling_method),
               '--chimera_method', str(chimera_method),
               '--min_parental_fold', str(min_fold_parent_over_abundance),
               '--allow_one_off', str(allow_one_off),
               '--num_threads', str(n_threads),
               '--learn_min_reads', str(n_reads_learn),
               '--homopolymer_gap_penalty', str(homopolymer_gap_penalty),
               '--band_size', str(band_size)]
        try:
            run_commands([cmd])
        except subprocess.CalledProcessError as e:
            if e.returncode == 2:
                raise ValueError(
                    "No reads passed the filter. trunc_len (%r) may be longer"
                    " than read lengths, or other arguments (such as max_ee"
                    " or trunc_q) may be preventing reads from passing the"
                    " filter." % trunc_len)
            else:
                raise Exception("An error was encountered while running DADA2"
                                " in R (return code %d), please inspect stdout"
                                " and stderr to learn more." % e.returncode)
        return _denoise_helper(biom_fp, track_fp, hashed_feature_ids)


def denoise_single(demultiplexed_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
                   trunc_len: int, trim_left: int = 0, max_ee: float = 2.0,
                   trunc_q: int = 2, pooling_method: str = 'independent',
                   chimera_method: str = 'consensus',
                   min_fold_parent_over_abundance: float = 1.0,
                   allow_one_off: bool = False,
                   n_threads: int = 1, n_reads_learn: int = 1000000,
                   hashed_feature_ids: bool = True
                   ) -> (biom.Table, DNAIterator, qiime2.Metadata):
    return _denoise_single(
        demultiplexed_seqs=demultiplexed_seqs,
        trunc_len=trunc_len,
        trim_left=trim_left,
        max_ee=max_ee,
        trunc_q=trunc_q,
        pooling_method=pooling_method,
        chimera_method=chimera_method,
        min_fold_parent_over_abundance=min_fold_parent_over_abundance,
        allow_one_off=allow_one_off,
        n_threads=n_threads,
        n_reads_learn=n_reads_learn,
        hashed_feature_ids=hashed_feature_ids,
        homopolymer_gap_penalty='NULL',
        band_size='16')