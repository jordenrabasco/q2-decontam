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
from qiime2.plugin.util import transform
from ._stats import ScoreTable, ScoreTableDirFmt, ScoreTableFormat

import biom
import skbio
import qiime2.util
import pandas as pd
import numpy as np
from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
                           Choices, MetadataColumn, Categorical, List,
                           Citations, TypeMatch)

from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import (
    FastqGzFormat, SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
from q2_types.feature_table import FeatureTable, Frequency


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
_COL_STR = (lambda x: x in { 'column_name', 'column_number'},
             'sample_name or column_name or column_number')
_DECON_METHOD_STR = (lambda x: x in {'frequency', 'prevalence', 'combined'},
             'freqeuncy, prevalence, combined')
_BOOLEAN = (lambda x: type(x) is bool, 'True or False')
# Better to choose to skip, than to implicitly ignore things that KeyError
_SKIP = (lambda x: True, '')
_valid_inputs = {
    'm_control_file': _SKIP,
    'threshold': _NAT_NUM,
    'control_sample_id_method': _COL_STR,
    'decon_method': _DECON_METHOD_STR,
    'control_column_id': _SKIP,
    'control_sample_indicator': _SKIP,
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
def _identify_helper(track_fp):

    df = pd.read_csv(track_fp, sep='\t', index_col=0)
    df.index.name = '#OTU ID'
    df=df.drop(df.columns[[(len(df.columns)-1)]], axis=1)
    temp_transposed_table = df.transpose()
    temp_transposed_table = temp_transposed_table.dropna()
    df = temp_transposed_table.transpose()

    metadata = transform(df, from_type=pd.DataFrame, to_type=ScoreTableFormat)

    return metadata


def prevalence_identify(asv_or_otu_table: pd.DataFrame, meta_data: qiime2.Metadata, threshold: float=0.1,
                   control_sample_id_method: str='column_name', control_column_id: str = 'NULL',control_sample_indicator: str='NULL'
                   ) -> (ScoreTableFormat):
    #_check_inputs(**locals())
    with tempfile.TemporaryDirectory() as temp_dir_name:
        biom_fp = os.path.join(temp_dir_name, 'output.tsv.biom')
        track_fp = os.path.join(temp_dir_name,'track.tsv')
        ASV_dest = os.path.join(temp_dir_name,'temp_ASV_table.csv')
        transposed_table =  asv_or_otu_table.transpose()
        transposed_table.to_csv(os.path.join(ASV_dest))

        metadata = meta_data.to_dataframe()
        meta_dest = os.path.join(temp_dir_name,'temp_metadata.csv')
        metadata.to_csv(os.path.join(meta_dest))

        cmd = ['run_decontam.R',
                   '--asv_table_path', str(ASV_dest),
                   '--threshold', str(threshold),
                   '--decon_method', str('prevalence'),
                   '--output_path', biom_fp,
                   '--output_track', track_fp,
                   '--temp_dir_name', temp_dir_name,
                   '--meta_table_path', str(meta_dest),
                   '--control_sample_id_method', str(control_sample_id_method),
                   '--control_column_id', str(control_column_id),
                   '--control_sample_indicator', str(control_sample_indicator)]
        try:
            run_commands([cmd])
        except subprocess.CalledProcessError as e:
            if e.returncode == 2:
                raise ValueError(
                        "There was an issue running run_decontam.R please check your inputs")
            else:
                raise Exception("An error was encountered while running Decontam"
                                    " in R (return code %d), please inspect stdout"
                                    " and stderr to learn more." % e.returncode)
        return _identify_helper(track_fp)

