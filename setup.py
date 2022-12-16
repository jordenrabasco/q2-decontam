# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer


setup(
    name="q2-decontam",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url="https://qiime2.org",
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Jorden Rabasco and Benjamin Callahan",
    author_email="jrabasc@ncsu.edu",
    description="Apply decontam to present or remove potential contmaination ASVs. ",
    scripts=['q2_decontam/assets/run_decontam.R'],
    package_data={
        'q2_decontam': ['citations.bib'],
        'q2_decontam._threshold_graph': ['assets/index.html'],
        'q2_decontam.tests': ['data/*',
                           'data/expected/*',
                              'data/tutorial_data/*']
    },
    entry_points={
        "qiime2.plugins":
        ["q2-decontam=q2_decontam.plugin_setup:plugin"]
    },
    zip_safe=False,
)
