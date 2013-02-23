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


#plotting bifurcation diagram for both ancestral and derived allele
#from the focal SNP at position 456
#which display a strong signal of selection
layout(matrix(1:2,2,1))
#ancestral allele
bifurcation.diagram(haplohh_cgu_bta12,mrk_foc=456,all_foc=1,nmrk_l=20,nmrk_r=20)
#derived allele
bifurcation.diagram(haplohh_cgu_bta12,mrk_foc=456,all_foc=2,nmrk_l=20,nmrk_r=20)