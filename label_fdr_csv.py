import sys
fdr_csv = sys.argv[1] #path to
fdr_csv_labeled = sys.argv[2]
fdr = sys.argv[3]

with open(fdr_csv ,'r') as input, open(fdr_csv_labeled,'w') as output:
    input1 = input.readlines()
    for i in input1:
        i=i.split(',')
        if 'e' in i[-1]:
            if int(i[-4]) < int(i[-3]):  #primary_count < alternative_count
                output.writelines('\t'.join(i).strip('\n')+'\t'+'hap2_type'+'\n')
            if int(i[-4]) > int(i[-3]):  #primary_count > hap2_count
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
# python ./label_fdr_csv.py < fdr_csv >  < fdr_csv_labeled > <fdr>
# python ./label_fdr_csv.py ./seedling/seedling.graph.count.fdr.csv ./seedling/seedling.graph.count.fdr001.label.csv 0.01