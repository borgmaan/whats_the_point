#!/usr/bin/env python
# Andrew Borgman
# Filtering SNPs and indels per the Jindo criteria 
import sys,os
from collections import OrderedDict
import math

""" SNP Filtering:
     > no indel within a 5-bp flanking region
     > estimated copy number of flanking sequences (<2) 
     > minimum distance between any two SNPs ( >= 5 bp) 
     > overall depth (=<100) at a given position in the reference. 
    
    INDEL Filtering:
     > Multiple indels occurring within a 20-bp window were filtered out
"""

# Allows for unbuffered output
class Unbuffered:
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

sys.stdout=Unbuffered(sys.stdout)

# Variant storage
snps = OrderedDict()
indels = OrderedDict()

#with open(sys.argv[1]) as infile:
while True:
    line = sys.stdin.readline()

    if not line:
        break    

    # Write out comment lines, add indels and snps to separate ordered dictionaries for filtering
    if not line.startswith("#"):
        spl = line.split()
        if len(spl[3]) > 1 or len(spl[4]) > 1:
            indels[int(spl[1])] = line
        else:
            snps[int(spl[1])] = line

# Rip through both dictionaries, delete those that don't pass filter criterion
all_positions = set(snps.keys()) |  set(indels.keys())

for pos in all_positions:

    # If the position is in both dicts, delete them both
    if (pos in snps) and (pos in indels):
        del snps[pos]
        del indels[pos]

    # SNP filtering
    elif pos in snps:
        for x in range(-5, 5):
            if x != 0:
                if (pos + x) in all_positions:
                    del snps[pos]
                    break

    # Indel filtering
    elif pos in indels:
        for x in range(-20, 20):
            if x != 0:
                if (pos + x) in indels:
                    del indels[pos]
                    del indels[pos + x]
                    break


# Print these out... doesn't matter if they aren't ordered
for indel in indels.keys():
    print indels[indel],

for snp in snps.keys():
    print snps[snp],
