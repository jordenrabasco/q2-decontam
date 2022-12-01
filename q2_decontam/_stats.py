# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import SemanticType, model
from q2_types.feature_data import FeatureData

#defines types for scoretable
ScoreTable = SemanticType('ScoreTable', variant_of=FeatureData.field['type'])
class ScoreTableFormat(model.TextFileFormat):
    def validate(*args):
        pass

ScoreTableDirFmt = model.SingleFileDirectoryFormat(
    'ScoreTableDirFmt', 'stats.tsv', ScoreTableFormat)