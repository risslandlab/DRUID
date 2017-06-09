#!/usr/bin/env python

# argv[1] intersect input

from __future__ import division

import sys
import csv
import re
import numpy as np

class LengthError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

features = {}
info = {}
with open(sys.argv[1], 'rb') as tsvin:
    tsvin = csv.reader(tsvin, delimiter = '\t')
    for row in tsvin:
	# extract chromosome and organism
	m = re.search('^(chr[0-9A-Z]+)([a-z])', row[0])
	chromosome = m.group(1)
	organism = m.group(2)
	
	# extract gene and transcript id and exon number
	m = re.search('gene_id "([^"]+)[a-z]"; transcript_id "[^"]+"; exon_number \d+; exon_id "([^.]+)\.(\d+)[a-z]"; gene_name "[^"]+"', row[8])
	gene = m.group(1)
	transcript = m.group(2)
	exon = m.group(3)

	feature = gene + "_" + organism + "_" + exon

	start = int(row[3]) - 1 # gtf is 1 based and inclusive on both ends

        if (feature) not in features:
            features[feature] = np.zeros(int(row[4]) - start, dtype = np.int32)
            info[feature] = [chromosome, organism, gene, transcript, exon]
        cov_start = max(start, int(row[10])) - start
        cov_end = min(int(row[4]), int(row[11])) - start
        if (cov_end - cov_start) != int(row[13]):
            raise LengthError(row)
        f = features[feature] 
        try:
            f[cov_start:cov_end] += int(row[12])
        except ValueError as e:
            print [row, f.size]

print ",".join(["chromosome", "organism", "gene", "transcript", "exon", "mean", "stdev", "size", "coverage"])
for f in features:
    x = features[f]
    print ",".join(info[f]+[str(x) for x in [np.mean(x), np.std(x), x.size, np.sum( x!=0 ) / x.size]])

