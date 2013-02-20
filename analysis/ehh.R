#!/usr/bin/env Rscript
# Andrew Borgman
# 2/19/2013
# EHH analysis from fastPHASE output
library(rehh)
setwd("~/projects/whats_the_point/data/")

#dat = data2haplohh("pointer_haps.inp", )

dat = haplohh_cgu_bta12

#computing EHH statisitics for the focal SNP at position 456
#which displays a strong signal of selection
ehh <- calc_ehh(haplohh_cgu_bta12, mrk=456)

# Will this save to github
