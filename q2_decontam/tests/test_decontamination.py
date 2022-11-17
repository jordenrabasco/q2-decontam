# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
#from pandas.testing import assert_frame_equal
import skbio
import os
import tempfile
import qiime2
import biom
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
from q2_decontam._stats import ScoreTable, ScoreTableDirFmt, ScoreTableFormat
from qiime2.plugin.util import transform

from q2_decontam import prevalence_identify
from q2_decontam._decontamination import _check_featureless_table


def _sort_feature_index(df):
    temp_sorted_df = df.sort_values(by=['#OTU ID'])
    return temp_sorted_df


class TestPrevalenceIdentify(TestPluginBase):
    package = 'q2_decontam.tests'

    def setUp(self):
        super().setUp()
        table = qiime2.Artifact.load(self.get_data_path('expected/decon_default_ASV_table.qza'))
        self.asv_table = table.view(qiime2.Metadata).to_dataframe()
        self.metadata_input = qiime2.Metadata.load(self.get_data_path('expected/test_metadata.tsv'))

    def test_defaults(self):
        exp_table = pd.read_csv(self.get_data_path('expected/prevalence-score-table.tsv'), sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table=temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()

        output_feature_table = prevalence_identify(asv_or_otu_table=self.asv_table, meta_data=self.metadata_input,
                                        threshold=0.1, control_sample_id_method='column_name',
                                        control_column_id='Sample_or_ConTrol', control_sample_indicator='Control')
        df_output_feature_table = transform(output_feature_table, from_type=ScoreTableFormat, to_type=pd.DataFrame)
        df_output_feature_table=df_output_feature_table.round(decimals=6)
        exp_table=exp_table.round(decimals=6)

        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name, 'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)

            self.assertEqual(test_table,expecter_table)
        #self.assertEqual(feature_stats, exp_md)


if __name__ == '__main__':
    unittest.main()
