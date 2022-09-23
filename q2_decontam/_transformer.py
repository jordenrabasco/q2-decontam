import qiime2

from q2_decontam import DecontamStatsFormat
from q2_decontam.plugin_setup import plugin


@plugin.register_transformer
def _1(ff: DecontamStatsFormat) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))


@plugin.register_transformer
def _2(obj: qiime2.Metadata) -> DecontamStatsFormat:
    ff = DecontamStatsFormat()
    obj.save(str(ff))
    return ff
