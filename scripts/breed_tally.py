#!/usr/bin/env python
# Andrew Borgman
# QC data for pointing study
import os, sys
from itertools import chain
from point_functions import *
sys.path.append("/share/gwas_runner/scripts")
from utils import *
sys.path.append('/grill/')
os.environ['DJANGO_SETTINGS_MODULE'] = 'bbq.settings'
from pork.models import *
from chicken.models import *


DATA_DIR = '/share/bar_data/'

bar_files = os.listdir(DATA_DIR)

strong_pointer_breeds = file_to_list('../breed_lists/strong_pointers')
versatile_dog_breeds = file_to_list('../breed_lists/versatile_dogs')

with open("../qc_data/pointer_tallies.tsv", "w") as pointer_tally:
	pointer_tally.write("Breed\tDogs Genotyped\tN_Passing_Call_Rate\n")

	with open("../qc_data/low_call_rate_pointers.tsv" ,"w") as bad_file:
		bad_file.write("Breed\tID\tCall_Rate\n")

		for pointer in chain(strong_pointer_breeds, versatile_dog_breeds):
			if check_breed_existance(pointer):
				genotyped_breed_aliases = [z for z in set(list(chain(*[x.alias_csv.split(",") for x in Dog.objects.filter(breed__name = pointer).distinct() if x != ""]))) if z in bar_files]
				genotyped_count = len(genotyped_breed_aliases)
				bad_count = 0
				for alias in genotyped_breed_aliases:
					call_rate = 1 - check_call_rate(alias)
					if call_rate < 0.90:
						print "Found low call rate dog",alias,'with call rate of', call_rate
						bad_file.write("\t".join([str(x) for x in [pointer, alias, call_rate]]) + "\n")
						bad_count += 1

				pointer_tally.write("\t".join([str(x) for x in [pointer, genotyped_count, (genotyped_count - bad_count)]]) + "\n")

with open("../breed_lists/non_pointing_breeds", "w") as non_pointers:
	with open("../qc_data/non_pointer_tallies.tsv", "w") as non_pointer_tally:
		non_pointer_tally.write("Breed\tDogs Genotyped\tN_Passing_Call_Rate\n")

		with open("../qc_data/low_call_rate_non_pointers.tsv" ,"w") as bad_file:
			bad_file.write("Breed\tID\tCall_Rate\n")

			all_breeds = Breed.objects.all()
			for b in all_breeds:
				
				# Exclude out those pointers 
				if b.name not in strong_pointer_breeds and b.name not in versatile_dog_breeds:
					genotyped_breed_aliases = [z for z in set(list(chain(*[x.alias_csv.split(",") for x in Dog.objects.filter(breed = b).distinct() if x != ""]))) if z in bar_files]

					genotyped_count = len(genotyped_breed_aliases)
					bad_count = 0
					for alias in genotyped_breed_aliases:
						call_rate = 1 - check_call_rate(alias)
						if call_rate < 0.90:
							print "Found low call rate dog",alias,'with call rate of', call_rate
							bad_file.write("\t".join([str(x) for x in [b.name, alias, call_rate]]) + "\n")
							bad_count += 1

					# Only use breeds with at least 6 good dogs
					if (genotyped_count - bad_count) > 5:
						non_pointer_tally.write("\t".join([str(x) for x in [b.name, genotyped_count, (genotyped_count - bad_count)]]) + "\n")
						non_pointers.write(b.name + "\n")

