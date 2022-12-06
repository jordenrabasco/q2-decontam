# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._decontamination import identify, remove
from ._version import get_versions
from ._stats import ScoreTable, ScoreTableDirFmt, ScoreTableFormat
from ._threshold_graph import (score_viz)


__version__ = get_versions()['version']
del get_versions

__all__ = ['identify','remove',
           'ScoreTable', 'ScoreTableFormat', 'ScoreTableDirFmt',
           'score_viz']
