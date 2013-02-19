#!/usr/bin/env python
# Andrew Borgman
# Pairwise Fst comparisons between all breeds
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

DATA_DIR = '/share/bar_data/'
PED_DIR = '/share/ped_files/'
TEMP_DIR = "/home/andrew/point/tmp/"
FST_DIR = "/home/andrew/point/pairwise_fst/"

# Store all pariwise comparisons that have 
# already been made so as to not duplicate
comparisons = set()

# Grab all breed data
strong_pointer_breeds = file_to_list('../breed_lists/strong_pointers')
versatile_dog_breeds = file_to_list('../breed_lists/versatile_dogs')
non_pointing_breeds = file_to_list('../breed_lists/non_pointing_breeds')
all_pointers = list(chain(strong_pointer_breeds, versatile_dog_breeds))
all_breeds = list(chain(strong_pointer_breeds, versatile_dog_breeds, non_pointing_breeds))

# Grab map data to write out markers in sorted order
map_data = []
with open("/share/canFam3.map") as map_file:
	for line in map_file:
		spl = line.replace("\n", "").split()
		if spl[0] == 'X' or spl[0] == 'Y':
			spl[0] = 39
		tup = [spl[1], int(spl[0]), int(spl[3])]
		map_data.append(tup)

map_data = sorted(map_data, key=itemgetter(1,2))

# Loop through all breeds, make some pairwise Fst calculations
with open('../qc_data/ped.log', 'w') as log_file:
	log_file.write("Comparison_Hash\tBreed One\tBreed Two\tBreed One Count\tBreed Two Count\n")
	for breed_uno in all_pointers:
		print comparisons
		print 'Doing pariwise comparisons for', breed_uno
		for breed_dos in all_breeds:
			if breed_uno != breed_dos:
				
				# Make unique strings for comparisons to track which ones we've done
				breed_str = breed_uno + "_" + breed_dos
				breed_str_rev = breed_dos + "_" + breed_uno

				# Make sure we haven't already made the comparison
				if breed_str not in comparisons and breed_str_rev not in comparisons:

					# Add the comparison
					comparisons.add(breed_str)

					# Grab all of the aliases for the two groups
					breed_uno_aliases = get_breed_aliases(breed_uno)
					breed_dos_aliases = get_breed_aliases(breed_dos)

					# Make sure we have at least 6 dogs in each breed
					if len(breed_uno_aliases) > 5 and len(breed_dos_aliases) > 5:

						# After edit, make sure we don't duplicate any analyses
						outfile_name = breed_uno.replace("(","_").replace(")","_").replace(" ","_") + "-" + breed_dos.replace("(","_").replace(")","_").replace(" ","_")
						outfile_name_rev = breed_dos.replace("(","_").replace(")","_").replace(" ","_") + "-" + breed_uno.replace("(","_").replace(")","_").replace(" ","_")

						# Check to see if we've run the Fst on this pair, if so, skip it
						if not os.path.exists(FST_DIR + outfile_name) and not os.path.exists(FST_DIR + outfile_name_rev):

							# Run a pairwise Fst comparison between the two groups, returns a dict keyed on marker names
							# with Fst vals as values. Also returns a hash of the ped used for logging
							pairwise_fst_vals, comparison_hash = pairwise_comp(breed_uno_aliases, breed_dos_aliases, TEMP_DIR)

							# Log which hash goes with which dog combo
							log_file.write("\t".join([comparison_hash, breed_uno, breed_dos, str(len(breed_uno_aliases)), str(len(breed_dos_aliases))]) + "\n")

							# Write to file in chromosomal order
							outfile_name = breed_uno.replace("(","_").replace(")","_").replace(" ","_") + "-" + breed_dos.replace("(","_").replace(")","_").replace(" ","_")
							with open(FST_DIR + outfile_name, "w") as outfile:
								for marker in map_data:
									if marker[0] in pairwise_fst_vals:
										outfile.write("\t".join([str(x) for x in [marker[0], marker[1], marker[2], pairwise_fst_vals[marker[0]]]]) + "\n")
									else:
										outfile.write("\t".join([str(x) for x in [marker[0], marker[1], marker[2], "NA"]]) + "\n")

						else:
							print "Already ran comparison for", breed_dos, 'and', breed_uno
