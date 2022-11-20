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
from q2_decontam import ScoreTable, ScoreTableFormat, ScoreTableDirFmt

_COL_OPT = {'column_name', 'column_number'}
_DECON_METHOD_OPT = {'frequency', 'prevalence', 'combined'}

plugin = qiime2.plugin.Plugin(
    name='decontam',
    version=q2_decontam.__version__,
    website='http://benjjneb.github.io/decontam/',
    package='q2_decontam',
    description=('Identify and/or removes contamination sequences from a seq by sample table'),
    short_description='Plugin for removal of contamination sequences',
    citations=qiime2.plugin.Citations.load('citations.bib', package='q2_decontam')
)


plugin.methods.register_function(
    function=q2_decontam.identify,
    inputs={'asv_or_otu_table': FeatureTable[Frequency]},
    parameters={ 'meta_data': Metadata,
                 'threshold': qiime2.plugin.Float,
                'decon_method': qiime2.plugin.Str %
                qiime2.plugin.Choices(_DECON_METHOD_OPT),
                'freq_concentration_column': qiime2.plugin.Str,
                'prev_control_or_exp_sample_column': qiime2.plugin.Str,
                'prev_control_sample_indicator': qiime2.plugin.Str,},
    outputs=[('score_table', FeatureData[ScoreTable])],
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
        'threshold': ('Select threshold cutoff for decontam algorithm scores'),
        'freq_concentration_column': ('Input column name that has concentration information for the samples'),
        'prev_control_or_exp_sample_column': ('Input column name containing experimental or control sample metadata'),
        'prev_control_sample_indicator': ('indicate the control sample identifier')
    },
    output_descriptions={
        'score_table': ('The resulting table of scores from the input ASV table')
    },
    name='Identify contaminants via the prevelance method',
    description=('This method identifies contaminant sequences from an '
                 'OTU or ASV table and reports them to the user')
)











plugin.visualizers.register_function(
    function=q2_decontam.graph,
    inputs={
        'table': FeatureTable[Frequency]
    },
    parameters={
        'metadata': Metadata,
    },
    name='Generate a heatmap representation of a feature table',
    description='Generate a heatmap representation of a feature table with '
                'optional clustering on both the sample and feature axes.\n\n'
                'Tip: To generate a heatmap containing taxonomic annotations, '
                'use `qiime taxa collapse` to collapse the feature table at '
                'the desired taxonomic level.',
    input_descriptions={
        'table': 'The feature table to visualize.'
    },
    parameter_descriptions={
        'metadata': 'Annotate the sample IDs with these metadata values. '
                    'When metadata is present and `cluster`=\'feature\', '
                    'samples will be sorted by the metadata values.',
    }
)



plugin.register_formats(ScoreTableFormat, ScoreTableDirFmt)
plugin.register_semantic_types(ScoreTable)
plugin.register_semantic_type_to_format(
    FeatureData[ScoreTable], ScoreTableDirFmt)
importlib.import_module('q2_decontam._transformer')
