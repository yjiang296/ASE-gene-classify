#!/usr/bin/env python3
import sys
import os
infile = sys.argv[1]
out = sys.argv[2]

with open(infile,'r') as in_1, open(out,'w') as out_1:
    for in_2 in in_1:
        in_3 = in_2.split()
        valid_reads = in_3[-2].count('.')+in_3[-2].count(',')+\
            in_3[-2].count('a')+in_3[-2].count('t')+in_3[-2].count('c')+in_3[-2].count('g')+\
            in_3[-2].count('A')+in_3[-2].count('T')+in_3[-2].count('C')+in_3[-2].count('G')
        if valid_reads >= 10:
            out_1.write(in_2)

##usage example: python ./extract_snp_with_10more_valid_reads.py graph.BM-11DAP.q30.sort.pileup.intersect graph.BM-11DAP.q30.sort.pileup.intersect.with10more.valid.reads

