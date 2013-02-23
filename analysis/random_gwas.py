#!/usr/bin/env python
# Andrew Borgman
# 2/19/2012
# Running conventional GWAS and EMMAX GWAS on pointer/non-pointer data with 
# randomly sampled dogs from each breed to have 6 dogs max per breed
#The study will include 17 pointing labs and 103 non-pointing labs

import os, sys
from itertools import chain
from operator import itemgetter
sys.path.append("/home/andrew/point/src/")
from point_functions import *
sys.path.append("/share/gwas_runner/scripts")
from utils import *
from gwas_functions import *
sys.path.append('/grill/')
os.environ['DJANGO_SETTINGS_MODULE'] = 'bbq.settings'
from pork.models import *
from chicken.models import *

# Get flash pointer aliases
flash_aliases = file_to_list("/home/andrew/point/flash_point/flash_aliases")

# Get genotyped labs
all_labs = get_breed_aliases('Labrador Retriever')

# Split into cases and controls
cases = []
controls = []
for lab in all_labs:
	if lab in flash_aliases:
		cases.append(lab)
	else:
		controls.append(lab)


# Tell me how many dogs I have
print "The study will include", len(cases), 'pointing labs and', len(controls),'non-pointing labs'
ped_hash, ped_loc = make_ped_file(cases, controls, "/home/andrew/point/flash_point/")
p2p = ped_loc 
os.system('plink --dog --map /share/clean.map --ped ' + p2p + ' --allow-no-sex --missing-genotype - --exclude /share/ultimateExclusionList.txt --no-fid --noweb --no-parents --no-sex --nonfounders --mind 0.05 --geno 0.1 --maf 0.05 --assoc --adjust')
os.system('bash /home/andrew/point/flash_point/run_emmax.sh /home/andrew/point/flash_point/%s.ped' % ped_hash)
