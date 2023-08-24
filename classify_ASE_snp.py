#!/usr/bin/env python3
from decimal import Decimal
import sys
import os
PATH_to_with10more_valid_reads = sys.argv[1]
PATH_to_with10more_valid_reads_count_csv = sys.argv[2]
with open(PATH_to_with10more_valid_reads,'r') as inputf,open(PATH_to_with10more_valid_reads_count_csv,'w') as outf:
    for i in inputf:
        i=i.split()
        B_count = i[-2].count('.') + i[-2].count(',')
        M_count = i[-2].count('a') + i[-2].count('A') + i[-2].count('c') + i[-2].count('C') + i[-2].count('g') + i[-2].count('G') + i[-2].count('t') + i[-2].count('T')
        valid_count = B_count + M_count
        B_ratio = B_count/valid_count
        M_ratio = M_count/valid_count

        #if B_ratio >= 0.9 or M_ratio >= 0.9:
            #outf.write(i[0]+'\t'+i[1]+'\t'+i[3]+'\t'+i[4]+'\t'+i[-3]+'\t'+repr(valid_count)+'\t'+repr(B_count)+'\t'+repr(M_count)+'\t'+str(round(B_ratio,2))+'\t'+str(round(M_ratio,2))+'\n')
        outf.write(i[0]+'\t'+i[1]+'\t'+i[3]+'\t'+i[4]+'\t'+i[-3]+'\t'+repr(valid_count)+'\t'+repr(B_count)+'\t'+repr(M_count)+'\t'+str(round(B_ratio,2))+'\t'+str(round(M_ratio,2))+'\n')
#usage example: python ./classify_ASE_snp.py b73.BM-10DAP.q30.sort.pileup.intersect.exon.with10more.valid.reads b73.BM-10DAP.exon.with10more.valid.reads.count.csv

