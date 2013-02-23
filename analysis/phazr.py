#!/usr/bin/env python
# Andrew Borgman
# 2/19/2012
# Exaining haplotype significance using iHS and XP-EHH statistics.
# Calculate these genome-wide to determine empirical significance
# First need to phase haplotypes for whole genome
# for chr in $(seq 1 38) ; do plink --dog --allow-no-sex --missing-genotype - --exclude /share/ultimateExclusionList.txt --no-fid --noweb --no-parents --no-sex --nonfounders --ped all_dogs.ped --map /share/canFam3.map --chr $chr --geno 0.1 --recode --out all_chroms/chr$chr\.unphased ; sed -i 's/-/0/g' all_chroms/chr$chr\.unphased.ped ; done
# Need to phase haplotypes within a given breed
import os, sys
from itertools import chain
from operator import itemgetter
from point_functions import *
sys.path.append("/share/gwas_runner/scripts")
from utils import *
sys.path.append('/grill/')
os.environ['DJANGO_SETTINGS_MODULE'] = 'bbq.settings'
from pork.models import *
from chicken.models import *

# Store all pariwise comparisons that have 
# already been made so as to not duplicate
comparisons = set()

# Grab all breed data
strong_pointer_breeds = file_to_list('../breed_lists/strong_pointers')
versatile_dog_breeds = file_to_list('../breed_lists/versatile_dogs')
non_pointing_breeds = file_to_list('../breed_lists/non_pointing_breeds')
all_breeds = list(chain(strong_pointer_breeds, versatile_dog_breeds, non_pointing_breeds))

# Grab all of the aliases
all_aliases = file_to_list('../phasing/all_aliases')

# Get in our directory
os.chdir("/home/andrew/point/phasing/")

# Loop through breeds, phase haplotypes within breeds
for b in all_breeds:

	# Make directory to store phased haplotypes
	breed_str = b.replace("(","_").replace(")","_").replace(" ","_")
	dir_loc = "./" + "by_breed/" + breed_str + "/"
	os.system('mkdir %s' % dir_loc)

	# Grab aliases 
	breed_aliases = get_breed_aliases(b)

	# Write out inclusion list for shapeit
	to_include = [x for x in breed_aliases if x in all_aliases]
	with open("inclusion_list", "w") as outfile:
		for z in to_include:
			outfile.write("%s\n" % z)

	# Phase haplotypes for all chromosomes within breed
	for i in range(1,39):
		run_string = "shapeit --input-ped all_chroms/chr%d.unphased.ped all_chroms/chr%d.unphased.map --output-max %schr%d.phased.haps  %schr%d.phased.sample --include-ind inclusion_list --thread 10" % (i, i, dir_loc, i, dir_loc, i)
		os.system(run_string)
		os.system("rm shapeit*")

