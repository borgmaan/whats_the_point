#!/usr/bin/env python
# Andrew Borgman
# Useful functions for pointing analysis

import os, sys
from itertools import chain
sys.path.append('/grill/')
os.environ['DJANGO_SETTINGS_MODULE'] = 'bbq.settings'
from pork.models import *
from chicken.models import *
sys.path.append("/share/gwas_runner/scripts")
from utils import *
import numpy as np
from scipy import stats

def get_breed_aliases(breed_name = ""):
	"""
		Returns all aliases for genotyped dogs with good call rates
	"""
	# Grab which dogs failed the call rate scan
	low_call_aliases = file_to_list("/home/andrew/point/qc_data/all_low_call_aliases")

	# Grab all bar files
	bar_files = os.listdir("/share/bar_data/")

	# Get all aliases for genotyped dogs with good call rates
	breed_aliases = [t for t in [z.replace(" ", "_") for z in set(list(chain(*[x.alias_csv.split(",") for x in Dog.objects.filter(breed__name = breed_name).distinct() if x != ""]))) if z in bar_files] if t not in low_call_aliases]

	return breed_aliases

def check_breed_existance(breed_name = ""):
	"""
		Checks to see if we have a breed in our DB with this name
	"""

	name_parts = breed_name.split()

	try:
		breed = Breed.objects.get(name = breed_name)
		
	except:
		print "Couldn't find matching breed for", breed_name
		potential_matches = []
		for part in name_parts:
			potential_matches.extend(Breed.objects.filter(name__icontains = part))

		print "Is it one of these", potential_matches

	if breed:
		return True
	else:
		return False

def check_call_rate(alias_string = ""):
	missing_markers = 0

	with open("/share/bar_data/" + alias_string) as infile:
		for line in infile:
			spl = line.split()
			if spl[2] == "-":
				missing_markers += 1

	return float(missing_markers) / 173662.0


def fst(case_vals, con_vals):
	"""
		Calculates Fst statistic between two sets of alleles 
		at a given marker
	"""

	# Remove missing alleles
	case_vals = case_vals[~np.isnan(case_vals)]
	con_vals = con_vals[~np.isnan(con_vals)]

	# Make sure we have data left for both markers
	if len(case_vals) != 0 and len(con_vals) != 0:

		# Calculate major/minor allele frequencies
		freqMinor1 = case_vals.sum() / (2 * len(case_vals))
		freqMajor1 = 1 - freqMinor1
		freqMinor2 = con_vals.sum() / (2 * len(con_vals))
		freqMajor2 = 1 - freqMinor2

		#Local expected heterozygosity
		hetExp1 = 1 - (freqMinor1**2 + freqMajor1**2)
		hetExp2 = 1 - (freqMinor2**2 + freqMajor2**2)

		#Allele frequencies over whole poplulation
		minorBar = (case_vals.sum() + con_vals.sum()) / ((2 * len(case_vals)) + (2 * len(con_vals))) 
		majorBar = (((2 * len(case_vals)) - case_vals.sum()) + ((2 * len(con_vals)) - con_vals.sum())) / ((2 * len(case_vals)) + ((2 * len(con_vals))))

		#Global heterozygosity index for subpopulations
		Hs = ((hetExp1 * len(case_vals)) + (hetExp2 * len(con_vals))) / (len(case_vals) + len(con_vals))

		#Global heterozygosity index for total population
		Ht = 2 * minorBar * majorBar

		#Calculating Fst
		if Ht > 0:
			Fst = (Ht - Hs) / Ht
		else:
			Fst = 0

	# Return NA if we didn't have sufficient data
	else:
		Fst = "NA"

	return Fst



def pairwise_comp(case_files = [], con_files = [], storage_dir = ""):
	"""
		Takes a list of case aliases, control aliases, and a directory to store ped files
	"""

	# Make a ped file for these bad boiz
	ped_hash, ped_loc = make_ped_file(case_files, con_files, storage_dir)

	# Recode it to additive coding
	ad_file = "%s%s.raw" % (storage_dir, ped_hash)
	if not os.path.exists(ad_file):
		os.system("plink --dog --map /share/canFam3.map --ped %s --allow-no-sex --missing-genotype - --exclude /share/ultimateExclusionList.txt --no-fid --noweb --no-parents --no-sex --nonfounders --recodeA --out %s%s" % (ped_loc, storage_dir, ped_hash))
	
	#Reading ped files into numpy arrays
	ped = np.loadtxt(ad_file,str)
	
	#Grabbing SNP list
	snps = ped[0, : ]
	snps = snps[6:len(snps)].tolist()

	#Getting rid of first row
	ped = ped[1:len(ped[ :, 0]),:]
	
	#Splitting 
	case_sel = np.where(ped[ :, 5] == '2')
	con_sel = np.where(ped[ :, 5] == '1')
	
	case_ped = ped[case_sel, : ]
	con_ped = ped[con_sel, : ]
	
	#Grabbing only SNPs
	case_ped = case_ped[:,:,6:case_ped.shape[2]]
	con_ped = con_ped[:,:,6:con_ped.shape[2]]
	
	# Converting NA values to NaN and changing type to float
	case_nan_sel = np.where(case_ped == "NA")
	con_nan_sel = np.where(con_ped == "NA")
	case_ped[case_nan_sel] = "NaN"
	con_ped[con_nan_sel] = "NaN"
	case_ped = case_ped.astype(np.float)
	con_ped = con_ped.astype(np.float)

	#Starting Fst Calculations
	fst_vals = {}
	for i in range(0,case_ped.shape[2]):
		
		# Grab marker name and calculate Fst
		snp_name = snps[i][:-2]
		marker_score = fst(case_ped[:,:,i],con_ped[:,:,i])

		# Store those bad bitches
		fst_vals[snp_name] = marker_score

	return fst_vals, ped_hash

def d_sub_i(breed_name = "", to_compare = []):
	"""
		Take a breed and a list of breeds to compare to
		outputs a dictionary of snps and d(i) values
	"""

	# Source and Fst files
	FST_DIR = "/home/andrew/point/pairwise_fst/"

	# Grab string format that file was stored in 
	breed_string = breed_name.replace("(","_").replace(")","_").replace(" ","_")

	# Get all pairwise fst files for breed of interest
	pairwise_files = [x for x in os.listdir(FST_DIR) if breed_string in x]

	# Storage: results[breed][snp] = dj(i)
	results = {}
	snps = set()

	# Grab string format that file was stored in for all other dogs
	file_strings = [x.replace("(","_").replace(")","_").replace(" ","_") for x in to_compare]

	for f in file_strings:
		if f != breed_string:

			# Check which way the file was named
			if breed_string + "-" + f in pairwise_files:
				infile_string = FST_DIR + breed_string + "-" + f
			elif f + "-" + breed_string in pairwise_files:
				infile_string = FST_DIR + f + "-" + breed_string
			else:
				print 'no data for',breed_string + "vs." + f 
				continue


			# Init dict		
			results[f] = {}

			# Load up Fst file
			ped = np.loadtxt(infile_string,str)

			# NA value handling
			nan_sel = np.where(ped == "NA")
			ped[nan_sel] = "NaN"
			fst_vals = ped[: , 3].astype(np.float)

			# Calculating E[Fst] and SD[Fst]
			comp_mean = stats.nanmean(fst_vals)
			comp_sd = stats.nanstd(fst_vals)

			# Normalizing Fst vals
			norm_fst = (fst_vals - comp_mean) / comp_sd

			# Store results
			for i in range(len(norm_fst)):
				results[f][ped[i, 0]] = norm_fst[i]
				snps.add(ped[i, 0])

	# Aggregate results from all comparisons into d(i) stat
	d_vals = {}
	for snp in snps:
		d_stat = 0
		for breed in results.keys():
			if snp in results[breed]:
				if ~np.isnan(results[breed][snp]):
					d_stat += results[breed][snp]

			else:
				print snp, 'not found in', breed_name, "-", breed
				raw_input()
		d_vals[snp] = d_stat

	return d_vals


