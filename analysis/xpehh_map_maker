#!/usr/bin/env python
# Andrew Borgman
# XPEHH attempt
import os,sys
from point_functions import *
from operator import itemgetter
from itertools import chain

# Storage
MAP_DIR = "/home/andrew/point/hap_analysis/map_files/"
POINT_DIR = "/home/andrew/point/hap_analysis/point_haps/"
NON_DIR = "/home/andrew/point/hap_analysis/non_haps/"

# Grabbing SNP order
map_data = []
with open("/share/canFam3.map") as map_file:
	for line in map_file:
		spl = line.replace("\n", "").split()
		if spl[0] == 'X' or spl[0] == 'Y':
			spl[0] = 39
		tup = [spl[1], int(spl[0]), int(spl[3])]
		map_data.append(tup)

map_data = sorted(map_data, key=itemgetter(1,2))

# Loop though all chromosomes and create XPEHH input files
chroms = range(2,39)

for chrom in chroms:
	point_hap_file = "/home/andrew/point/phasing/pointer_phased/chr%d.phased.haps" % chrom
	point_sample_file = "/home/andrew/point/phasing/pointer_phased/chr%d.phased.sample" % chrom
	non_hap_file = "/home/andrew/point/phasing/non_phased/chr%d.phased.haps" % chrom
	non_sample_file = "/home/andrew/point/phasing/non_phased/chr%d.phased.sample" % chrom

	# Grab haplotype info
	pointers, point_snps = parse_shapeit_files(point_sample_file, point_hap_file)
	non_pointers, non_snps = parse_shapeit_files(non_sample_file, non_hap_file)

	# Store SNP allele info to write out files for XPEHH
	snp_alleles = {}
	for marker in zip(point_snps, non_snps):
		# Sanity check
		if marker[0][0] == marker[1][0]:
			marker_snps = list(set([marker[0][1], marker[0][2], marker[1][1], marker[1][2]]))
			if marker[0][0] not in snp_alleles:
				snp_alleles[marker[0][0]] = {}
			if len(marker_snps) == 1:
				snp_alleles[marker[0][0]][marker_snps[0]] = "0" 
			elif len(marker_snps) == 2:
				snp_alleles[marker[0][0]][marker_snps[0]] = "0" 
				snp_alleles[marker[0][0]][marker_snps[1]] = "1" 
			else:
				print "Something is fucking up"

	# Write out map file
	with open("%s%d.ihsmap" % (MAP_DIR, chrom) ,"w") as map_file:
		for snp in point_snps:
			# Get positional info
			snp_name = snp[0]
			pos_info = [x for x in map_data if x[0] == snp_name][0]
			pos = float(pos_info[2])
			cm_pos = (pos / 1000000) * 0.97
			# Figure out how markers are coded
			alleles = snp_alleles[snp_name].keys()
			if len(alleles) == 1:
				zero_allele = alleles[0]
				map_file.write(" ".join([str(x) for x in [snp_name, pos, cm_pos, zero_allele, zero_allele]]) + "\n")
			else:
				zero_allele = [allele for allele, code in snp_alleles[snp_name].items() if code == "0"][0]
				one_allele = [allele for allele, code in snp_alleles[snp_name].items() if code == "1"][0]
				map_file.write(" ".join([str(x) for x in [snp_name, pos, cm_pos, zero_allele, one_allele]]) + "\n")

	# Write out pointer file
	with open("%s%d.ihshap" % (POINT_DIR, chrom) ,"w") as hap_file:
		with open("%s%d.ident" % (POINT_DIR, chrom) ,"w") as ind_file:
			for dog in pointers:
				hap_1 = list(pointers[dog][0])
				hap_2 = list(pointers[dog][1])

				# Sub in correct 0s and 1s based on what was written to map file
				for i in range(len(hap_1)):

					try:
						hap_1[i] = snp_alleles[point_snps[i][0]][hap_1[i]]	
						hap_2[i] = snp_alleles[point_snps[i][0]][hap_2[i]]
						
					except:
						print "Allele", hap_1[i], 'or', hap_2[i], "was not found in"
						print point_snps[i][0]
						raw_input()
						print snp_alleles[point_snps[i][0]]
						raw_input()

				# Write out haplotypes and identifiers
				hap_file.write(" ".join(hap_1) + "\n" + " ".join(hap_2) + "\n")
				ind_file.write(dog + "\n" + dog + "\n")
				


	# Write out non-pointer file
	with open("%s%d.ihshap" % (NON_DIR, chrom) ,"w") as hap_file:
		with open("%s%d.ident" % (NON_DIR, chrom) ,"w") as ind_file:
			for dog in non_pointers:
				hap_1 = list(non_pointers[dog][0])
				hap_2 = list(non_pointers[dog][1])

				# Sub in correct 0s and 1s based on what was written to map file
				for i in range(len(hap_1)):

					try:
						hap_1[i] = snp_alleles[point_snps[i][0]][hap_1[i]]	
						hap_2[i] = snp_alleles[point_snps[i][0]][hap_2[i]]
						
					except:
						print "Allele", hap_1[i], 'or', hap_2[i], "was not found in"
						print point_snps[i][0]
						raw_input()
						print snp_alleles[point_snps[i][0]]
						raw_input()

				# Write out haplotypes and identifiers
				hap_file.write(" ".join(hap_1) + "\n" + " ".join(hap_2) + "\n")
				ind_file.write(dog + "\n" + dog + "\n")

