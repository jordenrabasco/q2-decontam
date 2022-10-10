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
from q2_decontam import DecontamStats, DecontamStatsFormat, DecontamStatsDirFmt

_CONTROL_OPT = {'column_name', 'column_number'}

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
    function=q2_decontam.contaminant_prevelance,
    inputs={'asv_or_otu_table': FeatureTable[Frequency]},
    parameters={ 'meta_data': Metadata,
                 'control_sample_id_method': qiime2.plugin.Str %
                qiime2.plugin.Choices(_CONTROL_OPT),
                'control_column_id': qiime2.plugin.Str,
                'control_sample_indicator': qiime2.plugin.Str,},
    outputs=[('table', Metadata)],
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
        'control_sample_id_method': ('Select how to id experiment control'
                                     'sequences'),
        'control_column_id': ('Input control column identification'),
        'control_sample_indicator': ('indicate the control sample identifier')
    },
    output_descriptions={
        'table': ('The resulting metadata file indicating contaminant seqs')
    },
    name='Identify contaminants via the prevelance method',
    description=('This method identifies contaminant sequences from an '
                 'OTU or ASV table and reports them to the user')
)

plugin.register_formats(DecontamStatsFormat, DecontamStatsDirFmt)
plugin.register_semantic_types(DecontamStats)
plugin.register_semantic_type_to_format(
    SampleData[DecontamStats], DecontamStatsDirFmt)

#plugin.register_formats(DecontamStatsFormat, DecontamStatsDirFmt)
#plugin.register_semantic_types(DecontamStats)
#plugin.register_semantic_type_to_format(
#    SampleData[DecontamStats], DecontamStatsDirFmt)
importlib.import_module('q2_decontam._transformer')
