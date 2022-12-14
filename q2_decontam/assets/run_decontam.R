#!/usr/bin/env Rscript

#install.packages("optparse", repos='http://cran.us.r-project.org')

library("decontam")
library("optparse")

cat(R.version$version.string, "\n")
errQuit <- function(mesg, status=1) { message("Error: ", mesg); q(status=status) }

option_list = list(
  make_option(c("--asv_table_path"), action="store", default='NULL', type='character',
              help="File path to table in .csv format "),
  make_option(c("--threshold"), action="store", default='NULL', type='character',
              help="threshold value for decontam algorithm"),
  make_option(c("--decon_method"), action="store", default='NULL', type='character',
              help="algoithm mode"),
  make_option(c("--output_track"), action="store", default='NULL', type='character',
              help="File path to tracking tsv file. If already exists, will be overwritten"),
  make_option(c("--meta_table_path"), action="store", default='NULL', type='character',
              help="File path to metadata in .tsv format"),
  make_option(c("--freq_con_column"), action="store", default='NULL', type='character',
              help="Name of column for frequency method"),
  make_option(c("--prev_control_or_exp_sample_column"), action="store", default='NULL', type='character',
              help="Name of column for prevelance method"),
  make_option(c("--prev_control_sample_indicator"), action="store", default='NULL', type='character',
              help="Indicator to identify control samples")
)
opt = parse_args(OptionParser(option_list=option_list))


#--asv_table_path /Users/jrabasc/Desktop/temp_ASV_table.csv --meta_table_path /Users/jrabasc/Desktop/test_metadata.tsv --control_sample_indicator Control  --control_sample_id_method column_name --control_column_id Sample_or_ConTrol

# Assign each of the arguments, in positional order, to an appropriately named R variable

inp.loc <- opt$asv_table_path
threshold <- if(opt$threshold=='NULL') NULL else as.numeric(opt$threshold)
out.track <- opt$output_track
metadata.loc<-opt$meta_table_path
decon.mode<-opt$decon_method
prev.control.col <- opt$prev_control_or_exp_sample_column
prev.id.controls<-opt$prev_control_sample_indicator
freq.control.col<-opt$freq_con_column

#testing variables

#inp.loc <- '/Users/jrabasc/Desktop/temp_ASV_table.csv' 
#out.path <- '/Users/jrabasc/Desktop/temp_summary_decontam.tsv'
#threshold<-0.1
#out.track <- '/Users/jrabasc/Desktop/temp_decontam.tsv'
#metadata.loc<-'/Users/jrabasc/Desktop/temp_metadata.csv'
#decon.mode<-'combined'
#prev.control.col <- 'Sample_or_ConTrol'
#prev.id.controls<-'Control'
#freq.control.col<-'quant_reading'

if(!file.exists(inp.loc)) {
  errQuit("Input ASV table does not exist.")
}else if(!file.exists(metadata.loc)) {
  errQuit("Input metadata file does not exist.")
}else{
  print("Congrats your files exist")
}

meta_data_cols <-function(metadata_df, control.col){
  control_vec<-c()
  index<-0
  for (id in colnames(metadata_df)) {
    index=index+1
    if(tolower(id) == tolower(control.col)){
      control_vec<-metadata_df[,c(index)]
    }
  }
  return(control_vec)
}

outputer<-function(decon_output, out.track,asv_df, out.path){
  ### WRITE OUTPUT AND QUIT ###
  cat("7) Write output\n")
  write.table(decon_output, out.track, sep="\t",
              row.names=TRUE, col.names=NA, quote=FALSE)
  
  q(status=0)
}

asv_df <- read.csv(file = inp.loc)
rownames(asv_df) <- asv_df[, 1]  ## set rownames
asv_df <- asv_df[, -1]
numero_df <- as.matrix(sapply(asv_df, as.numeric)) 
metadata_df<-read.csv(file = metadata.loc)

if(decon.mode == 'prevalence'){
  control_vec <- meta_data_cols(metadata_df, prev.control.col)
  #genretates true/false vec for is contamination
  true_false_control_vec<-grepl(prev.id.controls,control_vec)
  # Prevalence-based contaminant classification
  prev_contam <- isContaminant(numero_df, neg=true_false_control_vec, threshold=threshold, detailed=TRUE, normalize=TRUE, method='prevalence')
  outputer(prev_contam, out.track,asv_df)
}else if(decon.mode == 'frequency'){
  control_vec <- meta_data_cols(metadata_df, freq.control.col)
  #genretates numeric vector for contamination analysis
  quant_vec<-as.numeric(control_vec)
  # Prevalence-based contaminant classification
  freq_contam <- isContaminant(numero_df, conc=quant_vec, threshold=threshold, detailed=TRUE, normalize=TRUE, method='frequency')
  outputer(freq_contam, out.track,asv_df)
}else{
  prev_control_vec <- meta_data_cols(metadata_df, prev.control.col)
  quant_control_vec <- meta_data_cols(metadata_df, freq.control.col)

  quant_vec<-as.numeric(quant_control_vec)
  true_false_control_vec<-grepl(prev.id.controls, prev_control_vec)
  
  comb_contam <- isContaminant(numero_df, neg=true_false_control_vec, conc=quant_vec, threshold=threshold, detailed=TRUE, normalize=TRUE, method='combined')
  outputer(comb_contam, out.track,asv_df)
}







  
