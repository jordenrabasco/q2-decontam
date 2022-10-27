import qiime2
import pandas as pd
#from q2_types.feature_data import (TaxonomyFormat, HeaderlessTSVTaxonomyFormat, TSVTaxonomyFormat)
from q2_decontam import DecontamStatsFormat
from q2_decontam.plugin_setup import plugin


def _dataframe_to_tsv_taxonomy_format(df):
    ff = DecontamStatsFormat()
    df.index.name = '#OTU ID'
    df['contaminant'] = df['contaminant'].astype('bool')
    df.to_csv(str(ff), sep='\t', header=True, index=True)
    return ff


@plugin.register_transformer
def _1(ff: DecontamStatsFormat) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))


@plugin.register_transformer
def _2(obj: qiime2.Metadata) -> DecontamStatsFormat:
    ff = DecontamStatsFormat()
    obj.save(str(ff))
    return ff

@plugin.register_transformer
def _3(df: pd.DataFrame) -> DecontamStatsFormat:

    return _dataframe_to_tsv_taxonomy_format(df)
