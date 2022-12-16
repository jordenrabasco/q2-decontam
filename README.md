# q2-decontam


This is a QIIME 2 plugin. For details on QIIME 2, see https://qiime2.org.

Decontam Tutorial: 

This short install and orientation tutorial is intended for those that already have qiime installed on their computer.

Installation:
1) Clone this github repo (https://github.com/jordenrabasco/q2-decontam.git)
2) Per the Qiime development guidance provided here (https://dev.qiime2.org/latest/tutorials/first-plugin-tutorial/) add the downloaded repo to your python path to make it discoverable 
   1) Navigate to your .bashrc file (.zshrc for mac users) which should be located in your primary user account directory (for me it's at /Users/jrabasc) and add this line to the end of the file: export PYTHONPATH=~PATHTODOWNLOADEDFOLDER  
My PATHTODOWNLOADEDFOLDER is /Users/jrabasc/Desktop/github/q2-decontam
3) Testart your terminal and navigate to the setup.py file in the downloaded repository folder
4) run command: "python setup.py install" followed by command "qiime dev refresh-cache"
   1) decontam should not be integrated into your local qiime deployment
   2) run command: "qiime decontam" to confirm successful integration

Running q2-decontam:

This plugin has 3 actions identify, score-viz, and remove, which are designed to be used in the given order
1) identify
   1) Inputs are an ASV/OTU table, and an associated metadate file with at least a column listing which samples are control/experimental and a column for concentration reads of all the samples.
   2) output is the standard output from the decontam R package minus the column indicating which samples were designated as contaminants
2) score-viz
   1) inputs are the decontam identify output and the OTU/ASV input table also used as input in decontam identify
   2) output is a histogram showing the distribution of ASVs decontam scores
   3) if the weighted option is used it will show the number of reads for the ASVs at each decontam score level
   4) this action is used to ascertain the most appropriate threshold for your data
3) remove
   1) inputs are the decontam identify output and the OTU/ASV input table (same as the score-viz action)
   2) output is an OTU/ASV table with those samples identified as contaminants by your set threshold removed

Examples of Commands:

Example data used in the commands below, sourced from the decontam oral contamination vignette, can be found in the downloaded repo in the data folder at ~/q2-decontam/q2_decontam/tests/data/
1)  qiime decontam identify --i-asv-or-otu-table feature-table-1.qza --m-meta-data-file test_metadata.tsv --o-score-table score_table.qza --p-freq-concentration-column quant_reading --p-prev-control-or-exp-sample-column Sample_or_ConTrol --p-prev-control-sample-indicator Control  --p-decon-method combined
2) qiime decontam score-viz --i-decon-identify-table score_table.qza --i-asv-or-otu-table feature-table-1.qza --p-threshold 0.01 --o-visualization vizualize_test.qzv --p-weighted
3) qiime decontam remove --i-decon-identify-table score_table.qza --i-asv-or-otu-table feature-table-1.qza --p-threshold 0.1 --o-no-contaminant-asv-table no_contam.qza