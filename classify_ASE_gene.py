import time
from tqdm.autonotebook import tqdm
from tqdm.auto import trange

import sys
gene_bedfile = sys.argv[1]
intersect = sys.argv[2]
out = sys.argv[3]

#status = ['B_type','M_type','Balanced']

def binary_search(set, num):
    tmp = [r[1] for r in set]
    tmp.append(set[-1][2])
    tmp = list(zip(tmp[0:len(tmp)-1], tmp[1:len(tmp)]))

    left = 0
    right = len(set)
    mid = (0+len(set)-1) // 2

    num = int(num)

    if num < int(tmp[0][0]) or num >= int(tmp[-1][1]):
        return None
    
    # print(tmp)

    while True:
        # print(mid)
        if int(tmp[mid][0]) > num:
            right = mid
            mid = (left+mid) // 2
        elif int(tmp[mid][1]) <= num:
            left = mid
            mid = (right-1+mid) // 2 + 1
        else:
            break
    if num < int(set[mid][2]):
        return mid

if __name__ == '__main__':
    print('Loading files...', end='')
    #with open('MB9DAP_peaks.narrowPeak','r') as f1, open('MB9DAP.pe.q10.sort.rmdup.shift.pileup.intersect','r') as f2:
    with open(gene_bedfile,'r') as f1, open(intersect,'r') as f2:
        #loc = {}
        f2_1 = f2.readlines()
        f1_1 = f1.readlines()
    print('done')
    f1_1 = [line.split() for line in tqdm(f1_1, desc='Splitting f1')]
    f2_1 = [line.split() for line in tqdm(f2_1, desc='Splitting f2')]

    index1 = {}
    index2 = {}
    token = None
    num = 0
    for line in tqdm(f1_1, desc='Creating index for f1'):
        start = line[0]
        if start != token:
            token = start
            index1.update({token:num})
        num += 1
    tmp = index1.values()
    tmp2 = list(tmp)
    tmp2.pop(0)
    tmp2.append(len(f1_1))
    index1 = {key:(start, end) for key, start, end in zip(index1.keys(), tmp, tmp2)}

    token = None
    num = 0
    for line in tqdm(f2_1, desc='Creating index for f2'):
        start = line[0]
        if start != token:
            token = start
            index2.update({token:num})
        num += 1
    tmp = index2.values()
    tmp2 = list(tmp)
    tmp2.pop(0)
    tmp2.append(len(f2_1))
    index2 = {key:(start, end) for key, start, end in zip(index2.keys(), tmp, tmp2)}

    for chr in tqdm(index1.keys(), 'Chr'):
        f1_2 = f1_1[index1[chr][0]:index1[chr][1]]
        f2_2 = f2_1[index2[chr][0]:index2[chr][1]]
    
        join_table = {}
        for line_num, line in tqdm(zip(range(len(f2_2)), f2_2), desc='Join Oprt', total=len(f2_2)):
            res = binary_search(f1_2, line[1])
            if res is not None:
                try:
                    join_table[res].append(line_num)
                except:
                    join_table.update({res:[line_num]})
        
        # print(join_table)

        for r, line_list in tqdm(join_table.items(), desc='Counting primary_count&alternative_count'):
            reads_count = 0
            i = f1_2[r]
            js = [f2_2[j_index] for j_index in line_list]
            for j in js:
                reads_count += int(j[-3])
            if reads_count >= 1:
                primary_count = 0
                alternative_count = 0
                for j in js:
                    primary_count += j[-2].count('.') + j[-2].count(',')
                    alternative_count += j[-2].count('a') + j[-2].count('A') + \
                            j[-2].count('c') + j[-2].count('C') + \
                            j[-2].count('g') + j[-2].count('G') + \
                            j[-2].count('t') + j[-2].count('T')
                #with open('./result.csv','a') as f3:
                with open(out,'a') as f3:
                    if reads_count >= 0:
                            #print(i[3]+'\t'+loc[i[3]]+'\t'+B+'\t'+M+'\t'+status[0])
                            f3.write(i[3]+'\t'+i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+str(len(js))+'\t'+repr(primary_count)+'\t'+repr(alternative_count)+'\n')




####----------------info--------------------####
#author:yjiang yjiang22@dingtalk.com
#usage: python classify_ASE_gene.py <gene_bedfile> <*pileup.intersect> <OUT.csv>
#example : python classify_AES_gene.py final_annotation_gene.bed seedling.graph.pileup.intersect seedling.graph.count.csv

