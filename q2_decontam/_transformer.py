import qiime2
import pandas as pd
from q2_decontam import DecontamScoreFormat
from q2_decontam.plugin_setup import plugin
import collections


_score_table_header = collections.OrderedDict([
     ('#OTU ID', str),
     ('freq', float),
     ('prev', float),
     ('p.freq', float),
     ('p.prev', float),
     ('p', float)])


def _dataframe_to_tsv_DecontamScore_format(df):
    ff = DecontamScoreFormat()
    df.to_csv(str(ff), sep='\t', header=True, index=True)
    return ff

def _DecontamScore_to_df(ff):
    temp_meta = qiime2.Metadata.load(str(ff))
    df = temp_meta.to_dataframe()
    return df


@plugin.register_transformer
def _1(ff: DecontamScoreFormat) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))


@plugin.register_transformer
def _2(obj: qiime2.Metadata) -> DecontamScoreFormat:
    ff = DecontamScoreFormat()
    obj.save(str(ff))
    return ff

@plugin.register_transformer
def _3(df: pd.DataFrame) -> DecontamScoreFormat:

    return _dataframe_to_tsv_DecontamScore_format(df)

@plugin.register_transformer
def _4(ff: DecontamScoreFormat) -> pd.DataFrame:
    return _DecontamScore_to_df(ff)

