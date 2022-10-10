#!/usr/bin/env Rscript

#install.packages("optparse", repos='http://cran.us.r-project.org')

library("decontam")
library("optparse")

cat(R.version$version.string, "\n")
errQuit <- function(mesg, status=1) { message("Error: ", mesg); q(status=status) }

option_list = list(
  make_option(c("--asv_table_path"), action="store", default='NULL', type='character',
              help="File path to table in .csv format "),
  make_option(c("--output_path"), action="store", default='NULL', type='character',
              help="File path to output tsv file. If already exists, will be overwritten"),
  make_option(c("--output_track"), action="store", default='NULL', type='character',
              help="File path to tracking tsv file. If already exists, will be overwritten"),
  make_option(c("--temp_dir_name"), action="store", default='NULL', type='character',
              help="temporary directory location address"),
  make_option(c("--meta_table_path"), action="store", default='NULL', type='character',
              help="File path to metadata in .tsv format"),
  make_option(c("--control_sample_id_method"), action="store", default='NULL', type='character',
              help="Method to identify control samples"),
  make_option(c("--control_column_id"), action="store", default='NULL', type='character',
              help="Method to identify control samples"),
  make_option(c("--control_sample_indicator"), action="store", default='NULL', type='character',
              help="Indicator to identify control samples")
)
opt = parse_args(OptionParser(option_list=option_list))


#--asv_table_path /Users/jrabasc/Desktop/temp_ASV_table.csv --meta_table_path /Users/jrabasc/Desktop/test_metadata.tsv --control_sample_indicator Control  --control_sample_id_method column_name --control_column_id Sample_or_ConTrol

# Assign each of the arguments, in positional order, to an appropriately named R variable

inp.loc <- opt$asv_table_path
#inp.loc <- '/Users/jrabasc/Desktop/temp_ASV_table.csv' 
out.path <- opt$output_path
out.track <- opt$output_track
temp.dir.name<-opt$temp_dir_name
metadata.loc<-opt$meta_table_path
#metadata.loc<-'/Users/jrabasc/Desktop/temp_metadata.csv'

how.id.controls<-opt$control_sample_id_method
#how.id.controls<-'column_name'

control.col <- opt$control_column_id
#control.col <- 'Sample_or_ConTrol'

id.controls<-opt$control_sample_indicator
#id.controls<-'Control'

if(!file.exists(inp.loc)) {
  errQuit("Input ASV table does not exist.")
}else if(!file.exists(metadata.loc)) {
  errQuit("Input metadata file does not exist.")
}else{
  print("Congrats your files exist")
}

prevelance <- function(asv_df, control_vec, id.controls) {
  #genretates true/false vec for is contamination
  true_false_control_vec<-grepl(id.controls,control_vec)
  # Prevalence-based contaminant classification
  numero_df <- as.matrix(sapply(asv_df, as.numeric)) 
  prev_contam <- isContaminant(numero_df, neg=true_false_control_vec, threshold=0.1, detailed=TRUE, normalize=TRUE, method='prevalence')
  return(prev_contam)
}

asv_df <- read.csv(file = inp.loc)
rownames(asv_df) <- asv_df[, 1]  ## set rownames
asv_df <- asv_df[, -1]  
metadata_df<-read.csv(file = metadata.loc)

if(how.id.controls == 'column_name'){
  index<-0
  for (id in colnames(metadata_df)) {
    index=index+1
    if(tolower(id) == tolower(control.col)){
      control_vec<-metadata_df[,c(index)]
      prev_contam <<-prevelance(asv_df,control_vec,id.controls)
    }
  }
}else if(how.id.controls == 'column_number'){
  control.col <- if(opt$control_column_id=='NULL') NULL else as.numeric(opt$control_column_id)
  if(ncol(metadata_df) >= int(control.col)){
    control_vec<-metadata_df[,c(int(control.col))]
    prev_contam <<-prevelance(asv_df,control_vec,id.controls)
  }else{
    print("Not a valid column number")
  }
}else{
  print("Something has gone horribly wrong")
}

### WRITE OUTPUT AND QUIT ###
cat("7) Write output\n")
write.table(seqtab.nochim, out.path, sep="\t",
            row.names=TRUE, col.names=col.names, quote=FALSE)
write.table(df,file='/Users/admin/new_file.csv',col.names=FALSE)

q(status=0)




  
