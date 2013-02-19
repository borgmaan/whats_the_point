#!/usr/bin/env python
# Andrew Borgman
# Running simulations to get estimates of variance in 
# Fst for small sample sizes
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
import numpy as np
from scipy import stats

# Directories
DATA_DIR = '/share/bar_data/'
TEMP_DIR = "/home/andrew/point/tmp/"
SIM_DIR = "/home/andrew/point/simulations/sim_peds/"

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

# Take some breeds with a bunch of dogs //Doberman Pinscher
breeds = ["Golden Retriever",  "Great Pyrenees"]

# Store deviances from true vals
diffs = {}

# Grab all of the aliases for the two groups
breed_uno_aliases = get_breed_aliases("Golden Retriever")
breed_dos_aliases = get_breed_aliases("Great Pyrenees")

# Get pairwise vals for all 
true_vals, _ = pairwise_comp(breed_uno_aliases, breed_dos_aliases, TEMP_DIR)

# Pull random samples of increasing size. Calculate differences from true values
# at each run. Calculate mean an SD of differences for each sample size.
# Store those and write to file.
for samp in range(6, 96, 6):
	
	# Store temp values
	temp_vals = {}
	
	# Take 50 random samples and calculate their deviance from the true values
	for i in range(50):
		
		# Grab some random samples
		uno_samps = []
		dos_samps = []
		uno_sel = noDups(samp, 0, len(breed_uno_aliases))
		dos_sel = noDups(samp, 0, len(breed_dos_aliases))

		for k in range(len(uno_sel)):
			uno_samps.append(breed_uno_aliases[uno_sel[k]])
			dos_samps.append(breed_dos_aliases[dos_sel[k]])

		# Calculate fst vals for the sample
		samp_fst, _ = pairwise_comp(uno_samps, dos_samps, SIM_DIR)

		# Calculate and store differences
		for marker in samp_fst.keys():
			if marker not in temp_vals:
				temp_vals[marker] = np.array([], dtype=np.float)

			if samp_fst[marker] == "NA":
				temp_vals[marker] = np.append(temp_vals[marker], np.nan)
			else:
				diff = true_vals[marker] - samp_fst[marker]
				temp_vals[marker] = np.append(temp_vals[marker], diff)

	# Get averages and SDs for differences and store them 
	diffs[samp] = {}
	for marker in temp_vals.keys():
		mean_diff = stats.nanmean(temp_vals[marker])
		sd_diff = stats.nanstd(temp_vals[marker])

		diffs[samp][marker] = [mean_diff, sd_diff]

with open("/home/andrew/point/simulations/sim_results.tsv", "w") as outfile:
	header_str = "SNP,CHR,BP," + ",".join([str(x) + "_mean" + "," + str(x) + "_sd" for x in range(6,96,6)]) + "\n"
	outfile.write(header_str)
	for marker in map_data:
		outfile.write(",".join([str(x) for x in [marker[0], marker[1], marker[2]]]) + ",")
		outstr = ""
		for samp in range(6, 96, 6):
			if marker[0] in diffs[samp]:
				outstr += ",".join([str(x) for x in [diffs[samp][marker[0]][0], diffs[samp][marker[0]][1]]]) + ','
			else:
				outstr += ",".join([str(x) for x in ["NA", "NA"]]) + ","

		outstr = outstr[:-1] + "\n"
		outfile.write(outstr)


