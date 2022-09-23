# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._decontamination import denoise_single
from ._version import get_versions
from ._stats import DecontamStats, DecontamStatsDirFmt, DecontamStatsFormat


__version__ = get_versions()['version']
del get_versions

__all__ = ['denoise_single',
           'DecontamStats', 'DecontamStatsFormat', 'DecontamStatsDirFmt']
