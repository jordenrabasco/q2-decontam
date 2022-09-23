# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import SemanticType, model
from q2_types.sample_data import SampleData


DecontamStats = SemanticType('DecontamStats', variant_of=SampleData.field['type'])


class DecontamStatsFormat(model.TextFileFormat):
    def validate(*args):
        pass


DecontamStatsDirFmt = model.SingleFileDirectoryFormat(
    'DecontamStatsDirFmt', 'stats.tsv', DecontamStatsFormat)