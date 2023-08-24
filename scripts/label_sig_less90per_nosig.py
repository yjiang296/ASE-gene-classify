#!/usr/bin/env python3
import sys
fdr001_reads20 = sys.argv[1]
out = sys.argv[2]

with open(fdr001_reads20) as input, \
    open(out,'w') as out:
    input1 = input.readlines()
    for i in input1:
        i = i.split()
        if i[-1].strip()=='hap1_type':
            if int(i[-5])/(int(i[-4])+int(i[-5])) >=0.9:
                out.writelines(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\t'+i[9]+'\n')
            else:
                out.writelines(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\t'+'NoSig'+'\n')
        if i[-1].strip()=='hap2_type':
            if int(i[-4])/(int(i[-4])+int(i[-5])) >=0.9:
                out.writelines(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\t'+i[9]+'\n')
            else:
                out.writelines(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\t'+'NoSig'+'\n')

        if i[-1].strip()=='NoSig':
            out.writelines(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\t'+i[9]+'\n')

#usage example: 
#python ./label_sig_less90per_nosig.py endosperm.graph.count.3snp.20read.fdr001.label.csv endosperm.graph.count.3snp.20read.fdr001.label.final.csv


#with open('/data21/yjiang22/kqq_rnaseq/endosperm/endosperm.hap1_classification.reads20.fdr001.labeled.csv') as input, \
#    open('/data21/yjiang22/kqq_rnaseq/endosperm/endosperm.hap1_class.reads20.fdr001.labeled.final.csv','a') as out:
#    input1 = input.readlines()
#    for i in input1:
#        i = i.split()
#        if i[-1].strip()=='hap1_type':
#            if int(i[-5])/(int(i[-4])+int(i[-5])) >=0.9:
#                out.writelines(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\n')
#            else:
#                out.writelines(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+'NoSig'+'\n')
#        if i[-1].strip()=='hap2_type':
#            if int(i[-4])/(int(i[-4])+int(i[-5])) >=0.9:
#                out.writelines(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\n')
#            else:
#                out.writelines(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+'NoSig'+'\n')
#
#        if i[-1].strip()=='NoSig':
#            out.writelines(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\n')
