#!/usr/bin/env python
# Andrew Borgman
# Aggregate iHs and XP-EHH data
import os,sys
from operator import itemgetter

# Grab all data
xpehh = {}
with open("/home/andrew/projects/whats_the_point/data/xpehh_data/all_results.txt") as infile:
	for line in infile:
		spl = line.replace("\n","").split()
		xpehh[spl[0]] = [spl[2], spl[3], spl[4]]


# Grabbing SNP order
map_data = []
with open("/home/andrew/canFam3.map") as map_file:
	for line in map_file:
		spl = line.replace("\n", "").split()
		if spl[0] == 'X' or spl[0] == 'Y':
			spl[0] = 39
		tup = [spl[1], int(spl[0]), int(spl[3])]
		map_data.append(tup)

map_data = sorted(map_data, key=itemgetter(1,2))

with open("/home/andrew/projects/whats_the_point/data/xpehh_results.tsv", "w") as outfile:
	outfile.write("SNP\tCHR\tBP\tiHs Pop 1\tiHs Pop 2\tXP-EHH\n")
	
	for marker in map_data:
		if marker[0] in xpehh:
			outfile.write("\t".join([str(x) for x in [marker[0], marker[1], marker[2], xpehh[marker[0]][0], xpehh[marker[0]][1], xpehh[marker[0]][2]]]) + "\n")	