import qiime2
import pandas as pd
from q2_decontam import ScoreTableFormat, SummaryTableFormat
from q2_decontam.plugin_setup import plugin
import collections


_score_table_header = collections.OrderedDict([
     ('#OTU ID', str),
     ('freq', float),
     ('prev', float),
     ('p.freq', float),
     ('p.prev', float),
     ('p', float)])


def _dataframe_to_tsv_ScoreTable_format(df):
    ff = ScoreTableFormat()
    df.to_csv(str(ff), sep='\t', header=True, index=True)
    return ff

def _ScoreTable_to_df(ff):
    temp_meta = qiime2.Metadata.load(str(ff))
    df = temp_meta.to_dataframe()
    return df


@plugin.register_transformer
def _1(ff: ScoreTableFormat) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))


@plugin.register_transformer
def _2(obj: qiime2.Metadata) -> ScoreTableFormat:
    ff = ScoreTableFormat()
    obj.save(str(ff))
    return ff

@plugin.register_transformer
def _3(df: pd.DataFrame) -> ScoreTableFormat:

    return _dataframe_to_tsv_ScoreTable_format(df)

@plugin.register_transformer
def _4(ff: ScoreTableFormat) -> pd.DataFrame:
    return _ScoreTable_to_df(ff)

