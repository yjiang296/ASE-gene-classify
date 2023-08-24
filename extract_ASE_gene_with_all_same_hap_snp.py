#!/usr/bin/env python3
import sys
import os
import csv
import json
PATH_to_genes_allele_specific_snp_wao1 = sys.argv[1]
PATH_to_ase_gene_hap_snp_count_csv = sys.argv[2]

with open(PATH_to_genes_allele_specific_snp_wao1,'r') as inputf, open(PATH_to_ase_gene_hap_snp_count_csv,'w') as out:
    dic = {}
    for i in inputf:
        i = i.split()
        
        if i[3] not in dic.keys():
            dic[i[3]] = []
        dic[i[3]].append(i[7])
    

    for key,value in dic.items():
        out.write("{}\t{}\t{}".format(key, value.count('hap1_type'), value.count('hap2_type'))+'\n')

# {
# 'evm.TU.chr1.5':['hap2_type'],
# 'evm.TU.chr1.6':['hap1_type','hap1_type'],
#
# }

# geneid    hap1_count  hap2_count
# evm.TU.chr1.5 1   0
# evm.TU.chr1.6 4   0
# evm.TU.chr1.7 10   1


##usage example: python ./extract_ASE_gene_with_all_same_hap_snp.py graph.hap1.genes.allele-specific-snp.wao1 graph.ase.gene.hap.snp.count.csv
