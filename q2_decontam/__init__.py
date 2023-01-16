# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._decontamination import decontam_identify, decontam_remove
from ._version import get_versions
from ._stats import DecontamScore, DecontamScoreDirFmt, DecontamScoreFormat
from ._threshold_graph import (decontam_score_viz)


__version__ = get_versions()['version']
del get_versions

__all__ = ['decontam_identify','decontam_remove',
           'DecontamScore', 'DecontamScoreFormat', 'DecontamScoreDirFmt',
           'decontam_score_viz']
