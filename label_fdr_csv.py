import sys
wo0_fdr_csv = sys.argv[1] #path to
wo0_fdr_csv_labeled = sys.argv[2]
fdr = sys.argv[3]

with open(wo0_fdr_csv ,'r') as input, open(wo0_fdr_csv_labeled,'a') as output:
    input1 = input.readlines()
    for i in input1:
        i=i.split(',')
        if 'e' in i[-1]:
            if int(i[-4]) < int(i[-3]):  #hap1_count < hap2_count
                output.writelines('\t'.join(i).strip('\n')+'\t'+'hap2_type'+'\n')
            if int(i[-4]) > int(i[-3]):  #hap1_count > hap2_count
                output.writelines('\t'.join(i).strip('\n')+'\t'+'hap1_type'+'\n')
        #
        else:
            #i[-1] = '%.4f' % i[-1].strip('\n')
            if float(i[-1]) < float(fdr) or float(i[-1]) == float(fdr):
                if int(i[-4]) < int(i[-3]):
                    output.writelines('\t'.join(i).strip('\n')+'\t'+'hap2_type'+'\n')
                if int(i[-4]) > int(i[-3]):
                    output.writelines('\t'.join(i).strip('\n')+'\t'+'hap1_type'+'\n')
            if float(i[-1]) > float(fdr):
                output.writelines('\t'.join(i).strip('\n')+'\t'+'NoSig'+'\n')

##uasge example: 
# python ./label_fdr_csv.py < wo0_fdr_csv >  < wo0_fdr_csv_labeled > <fdr>
# python ./label_fdr_csv.py ./seedling/seed.graph_classification.wo0_fdr.csv ./seedling/seed.graph_classification.wo0.fdr005.labeled.csv 0.05