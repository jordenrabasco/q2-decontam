# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import hashlib
import subprocess

import biom
import skbio
import qiime2.util
import pandas as pd

from q2_types.feature_data import DNAIterator
from q2_types.per_sample_sequences import (
    FastqGzFormat, SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)


def run_commands(cmds, verbose=True):
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


def _check_featureless_table(fp):
    with open(fp) as fh:
        # There is a header before the feature data
        for line_count, _ in zip(range(1, 3), fh):
            pass
    if line_count < 2:
        raise ValueError("No features remain after denoising. Try adjusting "
                         "your truncation and trim parameter settings.")


_WHOLE_NUM = (lambda x: x >= 0, 'non-negative')
_NAT_NUM = (lambda x: x > 0, 'greater than zero')
_POOL_STR = (lambda x: x in {'pseudo', 'independent'},
             'pseudo or independent')
_CHIM_STR = (lambda x: x in {'pooled', 'consensus', 'none'},
             'pooled, consensus or none')
_BOOLEAN = (lambda x: type(x) is bool, 'True or False')
# Better to choose to skip, than to implicitly ignore things that KeyError
_SKIP = (lambda x: True, '')
_valid_inputs = {

}


# TODO: Replace this with Range predicates when interfaces support them better
def _check_inputs(**kwargs):
    for param, arg in kwargs.items():
        check_is_valid, explanation = _valid_inputs[param]
        if not check_is_valid(arg):
            raise ValueError('Argument to %r was %r, should be %s.'
                             % (param, arg, explanation))


def _filepath_to_sample(fp):
    return fp.rsplit('_', 4)[0]



def denoise_single():
    # create an Empty DataFrame object
    df = pd.DataFrame()
    # append columns to an empty DataFrame
    df['Name'] = ['Anna', 'Pete', 'Tommy']
    df['Scores'] = [97, 600, 200]
    df['Questions'] = [2200, 75, 100]
    metadata = qiime2.Metadata(df)
    cmd = ['run_decontam.R']
    run_commands([cmd])
    return df
