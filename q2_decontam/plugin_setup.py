# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

import qiime2.plugin
from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
                           Choices, MetadataColumn, Categorical, List,
                           Citations, TypeMatch)
from q2_types.sample_data import SampleData
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency

import q2_decontam
from q2_decontam import DecontamScore, DecontamScoreFormat, DecontamScoreDirFmt

_DECON_METHOD_OPT = {'frequency', 'prevalence', 'combined'}

plugin = qiime2.plugin.Plugin(
    name='decontam',
    version=q2_decontam.__version__,
    website='https://github.com/jordenrabasco/q2-decontam',
    package='q2_decontam',
    description=('Identify and/or removes contamination sequences from a seq by sample table'),
    short_description='Plugin for removal of contamination sequences',
    citations=qiime2.plugin.Citations.load('citations.bib', package='q2_decontam')
)


plugin.methods.register_function(
    function=q2_decontam.decontam_identify,
    inputs={'asv_or_otu_table': FeatureTable[Frequency]},
    parameters={ 'meta_data': Metadata,
                'decon_method': qiime2.plugin.Str %
                qiime2.plugin.Choices(_DECON_METHOD_OPT),
                'freq_concentration_column': qiime2.plugin.Str,
                'prev_control_or_exp_sample_column': qiime2.plugin.Str,
                'prev_control_sample_indicator': qiime2.plugin.Str,},
    outputs=[('score_table', FeatureData[DecontamScore])],
    input_descriptions={
        'asv_or_otu_table': ('Table with presence counts in the matrix '
                             'rownames are sample id and column names are'
                             'seqeunce id')
    },
    parameter_descriptions={
        'meta_data': ('metadata file indicating which samples in the '
                           'experiment are control samples, '
                           'assumes sample names in file correspond '
                           'to ASV_or_OTU_table'),
        'decon_method': ('Select how to which method to id contaminants with'),
        'freq_concentration_column': ('Input column name that has concentration information for the samples'),
        'prev_control_or_exp_sample_column': ('Input column name containing experimental or control sample metadata'),
        'prev_control_sample_indicator': ('indicate the control sample identifier')
    },
    output_descriptions={
        'score_table': ('The resulting table of scores from the input ASV table')

    },
    name='Identify contaminants',
    description=('This method identifies contaminant sequences from an '
                 'OTU or ASV table and reports them to the user')
)

plugin.methods.register_function(
    function=q2_decontam.decontam_remove,
    inputs={'decon_identify_table': FeatureData[DecontamScore],
            'asv_or_otu_table': FeatureTable[Frequency]},
    parameters={'threshold': qiime2.plugin.Float},
    outputs=[('no_contaminant_asv_table', FeatureTable[Frequency])],
    input_descriptions={
        'decon_identify_table': ('Output table from decontam identify'),
        'asv_or_otu_table': ('Table with presence counts in the matrix '
                             'rownames are sample id and column names are'
                             'seqeunce id')
    },
    parameter_descriptions={
        'threshold': ('Select threshold cutoff for decontam algorithm scores')
    },
    output_descriptions={
        'no_contaminant_asv_table': ('The resulting table of scores once contaminants are removed')

    },
    name='Removes contaminant',
    description=('This method identifies contaminant sequences from an '
                 'OTU or ASV table and reports them to the user')
)


plugin.visualizers.register_function(
    function=q2_decontam.decontam_score_viz,
    inputs={
        'decon_identify_table': FeatureData[DecontamScore],
            'asv_or_otu_table': FeatureTable[Frequency]
    },
    parameters={
        'threshold':  qiime2.plugin.Float,
        'weighted': qiime2.plugin.Bool,
        'bin_size': qiime2.plugin.Float
    },
    name='Generate a histogram representation of the scores',
    description='Creates histogram based on the output of decontam identify',
    input_descriptions={
        'decon_identify_table': 'Output from decontam identify to be vizualized',
        'asv_or_otu_table': 'Raw OTU/ASV table that was used as input to identify'
    },
    parameter_descriptions={
        'threshold': ('Select threshold cutoff for decontam algorithm scores'),
        'weighted': ('weight the decontam scores by their assoicated read number'),
        'bin_size': ('Select bin size for the histogram')
    }
)



plugin.register_formats(DecontamScoreFormat, DecontamScoreDirFmt)
plugin.register_semantic_types(DecontamScore)
plugin.register_semantic_type_to_format(
    FeatureData[DecontamScore], DecontamScoreDirFmt)
importlib.import_module('q2_decontam._transformer')
