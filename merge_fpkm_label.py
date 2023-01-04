import pandas as pd
import sys
labeled_file = sys.argv[1] #path to
fpkm_file = sys.argv[2] #path to
out_merged = sys.argv[3]

label = pd.read_csv(labeled_file, sep='\t')

fpkm = pd.read_csv(fpkm_file, sep='\t')
merge = pd.merge(label,fpkm,how='left')

merge.to_csv(out_merged,sep='\t',index=False)

##usage example: python *py <*fdr001.labeled.csv> <gene_fpkm_matrix_tab.csv> <*fdr001.labeled.fpkm.csv>
#python ../merge_fpkm_label.py seed.graph_classification.wo0_fdr001.labeled.csv gene_fpkm_matrix_tab.csv seed.graph_classification.wo0_fdr001.labeled.fpkm.csv