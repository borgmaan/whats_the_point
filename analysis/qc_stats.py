#!/usr/bin/env python
# Andrew Borgman
# QC stats for RIP
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
import numpy as np
from scipy import stats
import math

DATA_DIR = '/share/bar_data/'
PED_DIR = '/share/ped_files/'
TEMP_DIR = "/home/andrew/point/tmp/"
FST_DIR = "/home/andrew/point/pairwise_fst/"

# Grab all breed data
strong_pointer_breeds = file_to_list('../breed_lists/strong_pointers')
versatile_dog_breeds = file_to_list('../breed_lists/versatile_dogs')
non_pointing_breeds = file_to_list('../breed_lists/non_pointing_breeds')
all_pointers = list(chain(strong_pointer_breeds, versatile_dog_breeds))
all_breeds = list(chain(strong_pointer_breeds, versatile_dog_breeds, non_pointing_breeds))


with open("/home/andrew/point/qc_data/breed_call_rates.csv","w") as outfile:
	outfile.write("Breed,Mean_Call,SD_Call,SE_Call,Type\n")
	for b in all_breeds:
		print b
		aliases = get_breed_aliases(b)
		if len(aliases) > 5:
			call_rates = []
			for a in aliases:
				try:
					call_rate = 1 - check_call_rate(a.replace("_", "\ "))
					call_rates.append(call_rate)
				except:
					pass

			call_rates = np.array(call_rates)
			call_mean = stats.nanmean(call_rates)
			call_sd = stats.nanstd(call_rates)
			call_se = call_sd / math.sqrt(len(aliases))

			if b in all_pointers:
				outfile.write(",".join([str(x) for x in [b, call_mean, call_sd, call_se, "Pointer"]]) + "\n")
			else:
				outfile.write(",".join([str(x) for x in [b, call_mean, call_sd, call_se, "Non-Pointer"]]) + "\n")