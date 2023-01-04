# Parental-specific-expression-gene-classify

#### These scripts are used to handle Hisat2's output bam file to classify parental specific expression genes.
###### Firstly, we used samtools(1.3.1) mpileup to get the detailed mapping information at each base.
###### Then, <get_intersect.pileup.py> was used to handle samtools mpileup output file to further extract the mapping information at the SNP sites between hap1 and hap2 genome.
###### Then, <classify_parent_specific_gene_MODIFIED.py> was used to count reads that support hap1 or hap2 specific expression.
###### Based on the reads count, we conducted binomial distribution hypothesis testing and Multiple hypothesis testing to classify hap1 or hap2 specific expression genes by <binom_fdr.R>.
###### Finally, <plot.R> was used to draw the scatter plots.
