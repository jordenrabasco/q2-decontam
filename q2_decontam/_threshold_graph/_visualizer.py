# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import pkg_resources

import q2templates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import qiime2
import q2_decontam


TEMPLATES = pkg_resources.resource_filename('q2_decontam._threshold_graph',
                                            'assets')
def score_viz(output_dir, decon_identify_table: qiime2.Metadata, asv_or_otu_table: pd.DataFrame, threshold: float=0.1):

    df = decon_identify_table.to_dataframe()
    values = df['p'].tolist()
    plt.xlim(0.0, 1.0)
    plt.xlabel('score value')
    plt.ylabel('number of ASVs')

    values = np.array(values)
    discreetlevel=10

    # Manually create `discreetlevel` bins anchored to  `threshold`
    binsAbove = round(discreetlevel * np.count_nonzero(values  > threshold) / values.size)
    binsBelow = discreetlevel - binsAbove
    binwidth = max((values.max() - threshold) / binsAbove,
                   (threshold - values.min()) / binsBelow)
    bins = np.concatenate([
        np.arange(threshold - binsBelow * binwidth, threshold, binwidth),
        np.arange(threshold, threshold + (binsAbove + 0.5) * binwidth, binwidth)
    ])

    h, bins, patches = plt.hist(values, bins)
    plt.setp([p for p, b in zip(patches, bins) if b < threshold], color='r', edgecolor="white", label="Contaminant ASVs")
    plt.setp([p for p, b in zip(patches, bins) if b >= threshold], color='b', edgecolor="white", label="True ASVs")
    plt.axvline(threshold, ymin=-.1,ymax=1.1 ,color='k', linestyle='dashed', linewidth=1, label="Threshold")

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc="upper left", framealpha=1)

    # Add a table at the bottom of the axes
    #cell_text = [[1,2,3]]
    #rows = ["test"]
    #columns = ('Contaminant ASVs', 'True ASVs', '% Contaminant')
    #plt.table(cellText=cell_text,
    #                colLabels=columns,
    #                cellLoc = 'center', rowLoc = 'center',
    #                loc='bottom', bbox=[0.25, -0.3, 0.3, 0.3])

    # Adjust layout to make room for the table:
    #plt.subplots_adjust(left=0.4, bottom=0.4)

    contam_asvs = 0
    true_asvs = 0
    for val in values:
        if val < threshold:
            contam_asvs = contam_asvs + 1
        else:
            true_asvs = true_asvs + 1

    percent_asvs = float(contam_asvs)/float((contam_asvs+true_asvs))



    for ext in ['png', 'svg']:
        img_fp = os.path.join(output_dir, 'identify-table-histogram.%s' % ext)
        plt.savefig(img_fp)
    index_fp = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index_fp, output_dir, context={'contamer': str(contam_asvs), 'truer': str(true_asvs), 'percenter': str(percent_asvs)})

