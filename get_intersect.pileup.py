#!/usr/bin/env python3
import sys
import os
snp = sys.argv[1]
pileup = sys.argv[2]
intersect = sys.argv[3]

loc = {}
#f1 = open('../Mo17_based_on_B73_nucmer.SNPs', 'r')
f1 = open(snp, 'r')
for i in f1:
        i = i.split()
        loc[i[0]+'\t'+i[1]] = i[2]+'\t'+i[3]+'\t'+i[4]
#f2 = open('../MB9DAP.sort.rmdup.pileup', 'r')
f2 = open(pileup, 'r')
for i in f2:
        i = i.split()
        a = i[0]+'\t'+i[1]
        if a in loc.keys():
            with open(intersect,'a') as f3:
                f3.writelines(a+'\t'+loc[a]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\n')

#usage python get_intersect.pileup.py <SNPs> <pileup> <pileup.intersect>
#python ../get_intersect.pileup.py ../contig_hap2_ragtag_hap1_dnadiff.SNPs endosperm.graph.wo0.pileup endosperm.graph.wo0.pileup.intersect
