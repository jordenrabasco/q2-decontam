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
from q2_decontam import DADA2Stats, DADA2StatsFormat, DADA2StatsDirFmt

_POOL_OPT = {'pseudo', 'independent'}
_CHIM_OPT = {'pooled', 'consensus', 'none'}

plugin = qiime2.plugin.Plugin(
    name='decontam',
    version=q2_decontam.__version__,
    website='http://benjjneb.github.io/decontam/',
    package='q2_decontam',
    description=('This QIIME 2 plugin wraps DADA2 and supports '
                 'sequence quality control for single-end and paired-end '
                 'reads using the DADA2 R library.'),
    short_description='Plugin for sequence quality control with DADA2.',
    citations=qiime2.plugin.Citations.load('citations.bib', package='q2_decontam')
)


plugin.methods.register_function(
    function=q2_decontam.denoise_single,
    inputs={},
    parameters={},
    outputs=[('table', FeatureTable[Frequency])],
    input_descriptions={
    },
    parameter_descriptions={

    },
    output_descriptions={'table': 'The resulting feature table.'

    },
    name='Denoise and dereplicate single-end sequences',
    description=('This method denoises single-end sequences, dereplicates '
                 'them, and filters chimeras.')
)


plugin.register_formats(DADA2StatsFormat, DADA2StatsDirFmt)
plugin.register_semantic_types(DADA2Stats)
plugin.register_semantic_type_to_format(
    SampleData[DADA2Stats], DADA2StatsDirFmt)
importlib.import_module('q2_decontam._transformer')
