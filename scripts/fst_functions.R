#!/usr/bin/env Rscript
# Andrew Borgman
# 2/15/2013
# Helper functions for Fst and d(i) analysis

# Takes a directory of Fst or d(i) files, returns
# a single data frame
df_from_dir <- function(directory=""){

  # Grab all file names from directory
  files = sort(paste(directory, dir(directory), sep=""))

  # Start off the DF with the first set of vals
  big = read.table(files[1], stringsAsFactors=F, sep="\t")
  f = strsplit(files[1], split="/")[[1]][length(strsplit(files[1], split="/")[[1]])]
  names(big) = c("SNP", "CHR",  "BP", f)
  
  # Loop through rest of files, append values and name column
  for (i in 2:length(files)){
    f = strsplit(files[i], split="/")[[1]][length(strsplit(files[i], split="/")[[1]])]
    temp = read.table(files[i], stringsAsFactors=F, sep="\t")[,4]
    big$temp = temp
    names(big)[length(names(big))] = f
  }  
  
  return(big)
}