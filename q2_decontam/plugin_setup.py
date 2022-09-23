# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

import qiime2.plugin
from q2_types.per_sample_sequences import (
    SequencesWithQuality, PairedEndSequencesWithQuality)
from q2_types.sample_data import SampleData
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency

import q2_decontam
from q2_decontam import DecontamStats, DecontamStatsFormat, DecontamStatsDirFmt

_POOL_OPT = {'pseudo', 'independent'}
_CHIM_OPT = {'pooled', 'consensus', 'none'}

plugin = qiime2.plugin.Plugin(
    name='decontam',
    version=q2_decontam.__version__,
    website='http://benjjneb.github.io/decontam/',
    package='q2_decontam',
    description=('Place holder'),
    short_description='Place holder',
    citations=qiime2.plugin.Citations.load('citations.bib', package='q2_decontam')
)


plugin.methods.register_function(
    function=q2_decontam.denoise_single,
    inputs={'demultiplexed_seqs': SampleData[SequencesWithQuality |
                                             PairedEndSequencesWithQuality]},
    parameters={'trunc_len': qiime2.plugin.Int,
                'trim_left': qiime2.plugin.Int,
                'max_ee': qiime2.plugin.Float,
                'trunc_q': qiime2.plugin.Int,
                'pooling_method': qiime2.plugin.Str %
                qiime2.plugin.Choices(_POOL_OPT),
                'chimera_method': qiime2.plugin.Str %
                qiime2.plugin.Choices(_CHIM_OPT),
                'min_fold_parent_over_abundance': qiime2.plugin.Float,
                'allow_one_off': qiime2.plugin.Bool,
                'n_threads': qiime2.plugin.Int,
                'n_reads_learn': qiime2.plugin.Int,
                'hashed_feature_ids': qiime2.plugin.Bool},
    outputs=[('table', FeatureTable[Frequency]),
             ('representative_sequences', FeatureData[Sequence]),
             ('denoising_stats', SampleData[DecontamStats])],
    input_descriptions={
        'demultiplexed_seqs': ('The single-end demultiplexed sequences to be '
                               'denoised.')
    },
    parameter_descriptions={
        'trunc_len': ('Position at which sequences should be truncated due to '
                      'decrease in quality. This truncates the 3\' end of the '
                      'of the input sequences, which will be the bases that '
                      'were sequenced in the last cycles. Reads that are '
                      'shorter than this value will be discarded. If 0 is '
                      'provided, no truncation or length filtering will be '
                      'performed'),
        'trim_left': ('Position at which sequences should be trimmed due to '
                      'low quality. This trims the 5\' end of the '
                      'of the input sequences, which will be the bases that '
                      'were sequenced in the first cycles.'),
        'max_ee': ('Reads with number of expected errors higher than this '
                   'value will be discarded.'),
        'trunc_q': ('Reads are truncated at the first instance of a quality '
                    'score less than or equal to this value. If the resulting '
                    'read is then shorter than `trunc_len`, it is discarded.'),
        'pooling_method': (
            'The method used to pool samples for denoising. '
            '"independent": Samples are denoised independently. '
            '"pseudo": The pseudo-pooling method is used to '
            'approximate pooling of samples. In short, samples '
            'are denoised independently once, ASVs detected '
            'in at least 2 samples are recorded, and samples '
            'are denoised independently a second time, but '
            'this time with prior knowledge of the recorded '
            'ASVs and thus higher sensitivity to those ASVs.'
        ),
        'chimera_method': ('The method used to remove chimeras. '
                           '"none": No chimera removal is performed. '
                           '"pooled": All reads are pooled prior to chimera '
                           'detection. "consensus": Chimeras are detected in '
                           'samples individually, and sequences found '
                           'chimeric in a sufficient fraction of samples are '
                           'removed.'),
        'min_fold_parent_over_abundance': (
            'The minimum abundance of potential parents of a sequence being '
            'tested as chimeric, expressed as a fold-change versus the '
            'abundance of the sequence being tested. Values should be greater '
            'than or equal to 1 (i.e. parents should be more abundant than '
            'the sequence being tested). This parameter has no effect if '
            'chimera_method is "none".'),
        'allow_one_off': (
            'Bimeras that are one-off from exact are also '
            'identified if the `allow_one_off` argument is True.'
            'If True, a sequence will be identified as bimera if it is one '
            'mismatch or indel away from an exact bimera.'),
        'n_threads': ('The number of threads to use for multithreaded '
                      'processing. If 0 is provided, all available cores will '
                      'be used.'),
        'n_reads_learn': ('The number of reads to use when training the '
                          'error model. Smaller numbers will result in a '
                          'shorter run time but a less reliable error '
                          'model.'),
        'hashed_feature_ids': ('If true, the feature ids in the resulting '
                               'table will be presented as hashes of the '
                               'sequences defining each feature. The hash '
                               'will always be the same for the same sequence '
                               'so this allows feature tables to be merged '
                               'across runs of this method. You should only '
                               'merge tables if the exact same parameters are '
                               'used for each run.')
    },
    output_descriptions={
        'table': 'The resulting feature table.',
        'representative_sequences': ('The resulting feature sequences. Each '
                                     'feature in the feature table will be '
                                     'represented by exactly one sequence.')
    },
    name='Denoise and dereplicate single-end sequences',
    description=('This method denoises single-end sequences, dereplicates '
                 'them, and filters chimeras.')
)


plugin.register_formats(DecontamStatsFormat, DecontamStatsDirFmt)
plugin.register_semantic_types(DecontamStats)
plugin.register_semantic_type_to_format(
    SampleData[DecontamStats], DecontamStatsDirFmt)
importlib.import_module('q2_decontam._transformer')
