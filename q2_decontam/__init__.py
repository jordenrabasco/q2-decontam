# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._decontamination import denoise_single
from ._version import get_versions
from ._stats import DADA2Stats, DADA2StatsDirFmt, DADA2StatsFormat


__version__ = get_versions()['version']
del get_versions

__all__ = ['denoise_single',
           'DADA2Stats', 'DADA2StatsFormat', 'DADA2StatsDirFmt']
