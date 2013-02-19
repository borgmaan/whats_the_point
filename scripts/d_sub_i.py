#!/usr/bin/env python
# Andrew Borgman
# Calculating d(i) values for pointers vs. non-pointers and pointers vs. pointers
import os
from itertools import chain
from operator import itemgetter
from point_functions import *
sys.path.append("/share/gwas_runner/scripts")
from utils import *
sys.path.append('/grill/')
os.environ['DJANGO_SETTINGS_MODULE'] = 'bbq.settings'
from pork.models import *
from chicken.models import *

# Folders
FST_DIR = "/home/andrew/point/pairwise_fst/"
WITHIN_DIR = "/home/andrew/point/pointer_vs_pointer/"
BETWEEN_DIR = "/home/andrew/point/pointer_vs_non/"

# Grab all breed data
strong_pointer_breeds = file_to_list('../breed_lists/strong_pointers')
versatile_dog_breeds = file_to_list('../breed_lists/versatile_dogs')
non_pointing_breeds = file_to_list('../breed_lists/non_pointing_breeds')
all_pointers = list(chain(strong_pointer_breeds, versatile_dog_breeds))

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


for pointer in all_pointers:
	print "aggregating values for",pointer

	between_vals = d_sub_i(pointer, non_pointing_breeds)
	within_vals = d_sub_i(pointer, all_pointers)

	outfile_name = pointer.replace("(","_").replace(")","_").replace(" ","_")
	with open(WITHIN_DIR + outfile_name, "w") as outfile:
		for marker in map_data:
			if marker[0] in within_vals:
				outfile.write("\t".join([str(x) for x in [marker[0], marker[1], marker[2], within_vals[marker[0]]]]) + "\n")
			else:
				outfile.write("\t".join([str(x) for x in [marker[0], marker[1], marker[2], "NA"]]) + "\n")

	with open(BETWEEN_DIR + outfile_name, "w") as outfile:
		for marker in map_data:
			if marker[0] in between_vals:
				outfile.write("\t".join([str(x) for x in [marker[0], marker[1], marker[2], between_vals[marker[0]]]]) + "\n")
			else:
				outfile.write("\t".join([str(x) for x in [marker[0], marker[1], marker[2], "NA"]]) + "\n")
