# Parental-specific-expression-gene-classify

#### These scripts are used to handle Hisat2's output bam file to classify parental specific expression genes.
###### (1) we used samtools(1.3.1) mpileup to get the detailed mapping information at each base.
###### example:
###### bedtools bamtobed -i seedling.graph.bam > seedling.graph.bed
###### samtools mpileup -l seedling.graph.bed -f ../ragtag.scaffold.fasta seedling.graph.bam -o seedling.graph.pileup
###### (2) <get_intersect.pileup.py> was used to handle samtools mpileup output file to further extract the mapping information at the SNP sites between primary and alternative assembly.
###### (3) <classify_ASE_gene.py> was used to count reads that support primary or altenative assembly specific expression.
###### example:
###### python classify_ASE_gene.py final_annotation_gene.bed seedling.graph.pileup.intersect seedling.graph_count.csv

###### (4) Based on the reads count, we conducted binomial distribution hypothesis testing and multiple hypothesis testing to classify primary or altenative assembly specific expression genes by <binom_fdr.R>.


![流程图1 drawio](https://user-images.githubusercontent.com/115337217/212926971-0587592c-8395-40f1-9ab9-4e5c7af23972.png)
