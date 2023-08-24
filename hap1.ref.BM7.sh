#!/bin/bash

hisat2_index[0]=/data21/yjiang22/hybrid_hap_graph_application/parent-specific-acr_analysis/hap1_atac/hisat2_build/hap1
hisat2_index[1]=/data21/yjiang22/new_hybrid_hap_graph_BM/hisat2_index/new_pseudo_hap2_like_hap1/hap1_hap2_like_hap1

scripts_get_intersect_pileup_py=/data21/yjiang22/hybrid_hap_graph_application/scripts/get_intersect.pileup.py
scripts_extract_snp_with_10more_valid_reads_py=/data21/yjiang22/hybrid_hap_graph_application/scripts/extract_snp_with_10more_valid_reads.py
scripts_classify_ASE_gene_py=/data21/yjiang22/hybrid_hap_graph_application/scripts/classify_ASE_gene.py
scripts_ASE_snp_py=/data21/yjiang22/hybrid_hap_graph_application/scripts/classify_ASE_snp.py
scripts_binom_fdr_R=/data21/yjiang22/new_hybrid_hap_graph_BM/scripts/binom_fdr.R
scripts_chisq_BH_R=/data21/yjiang22/new_hybrid_hap_graph_BM/scripts/chisq_BH.R
scripts_label_fdr_csv_py=/data21/yjiang22/hybrid_hap_graph_application/scripts/label_fdr_csv.py
scripts_label_sig_less90per_nosig_py=/data21/yjiang22/hybrid_hap_graph_application/scripts/label_sig_less90per_nosig.py
scripts_filter_eFDR_R=/data21/yjiang22/hybrid_hap_graph_application/scripts/filter_eFDR.R
scripts_extract_ASE_gene_with_all_same_hap_snp_py=/data21/yjiang22/hybrid_hap_graph_application/scripts/extract_ASE_gene_with_all_same_hap_snp.py
scripts_extract_NM_from_bam_py=/data21/yjiang22/new_hybrid_hap_graph_BM/scripts/extract_mismatch_NM_from_bam.py

PATH_to_relative_R_WD=./allele-specific-gene/linear_ref/hap1_ref/BM7  ##need to be modified

DAP=7

rnaseq_em_fp_R1=/data21/yjiang22/hybrid_em_graph_application/rnaseq_meng/BM_0${DAP}_EM_1.fastp.fastq.gz
rnaseq_em_fp_R2=/data21/yjiang22/hybrid_em_graph_application/rnaseq_meng/BM_0${DAP}_EM_2.fastp.fastq.gz

##    source activate rnaseq
##    ##fastp --in1 ../../../rnaseq_meng/BM_07_EM_1.fastq.gz --in2 ../../../rnaseq_meng/BM_07_EM_2.fastq.gz --out1 ../../rnaseq_meng/BM_07_EM_1.fastp.fastq.gz --out2 ../../rnaseq_meng/BM_07_EM_2.fastp.fastq.gz && echo '--fastp done--'
##    hisat2 -p 30 -t -x ${hisat2_index[0]} -1 ${rnaseq_em_fp_R1} -2 ${rnaseq_em_fp_R2} | samtools view -Sbq 30 > hap1.BM-${DAP}DAP.q30.bam && echo "----hisat2 done----"
##    samtools sort -@ 30 -o hap1.BM-${DAP}DAP.q30.sort.bam hap1.BM-${DAP}DAP.q30.bam && echo '--samtools sort done--'
##    samtools index hap1.BM-${DAP}DAP.q30.sort.bam && echo '--samtools index done--'
##    
##    stringtie -p 30 -G /data21/yjiang22/hybrid_hap_graph_application/EVM.all_1new.gtf -e -A hap1.BM-${DAP}DAP.abundance.txt hap1.BM-${DAP}DAP.q30.sort.bam && echo '--stringtie done--'
##    
##    #step1 samtools mpileup
##    bedtools bamtobed -i hap1.BM-${DAP}DAP.q30.sort.bam > hap1.BM-${DAP}DAP.q30.sort.bed && echo '--step1.1 bam2bed done--'
##    
##    source activate base ## base下面重新安装了samtool1.18，旧版本unrecognize mpileup新加的参数
##    samtools mpileup --no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends -l hap1.BM-${DAP}DAP.q30.sort.bed -f /data21/yjiang22/hybrid_hap_graph_application/ASE_analysis/hisat2_build_new/hap1_1to10.fa hap1.BM-${DAP}DAP.q30.sort.bam -o hap1.BM-${DAP}DAP.q30.sort.pileup && echo '--step1.2 mpileup done--'
##    
##    #step2 extract mapping info at snp sites
##    #grep -v 'scaf' b73v5-mo17.snps.all2vcf.vcf.snp > b73v5-mo17.snps.all2vcf.vcf.chr.snp
##    awk '$4!=0 {print $0}' hap1.BM-${DAP}DAP.q30.sort.pileup > hap1.BM-${DAP}DAP.q30.sort.wo0.pileup && echo "--step2.1 hap1.BM-${DAP}DAP awk wo0 done---"
##    mv hap1.BM-${DAP}DAP.q30.sort.wo0.pileup hap1.BM-${DAP}DAP.q30.sort.pileup ##save storage space
##    
##    ##只考虑gene exon上snp处的reads
##    touch hap1.BM-${DAP}DAP.q30.sort.pileup.intersect.exon;rm hap1.BM-${DAP}DAP.q30.sort.pileup.intersect.exon #防止下一步追加写

python ${scripts_get_intersect_pileup_py} /data21/yjiang22/new_hybrid_hap_graph_BM/hisat2_index/snp_filter/03_remove_HIFI_non_heterozygous/hap1_hap2.snps.rm.ajacent.with3more.hifi.heterozygous.in.exon hap1.BM-${DAP}DAP.q30.sort.pileup ./hap1.BM-${DAP}DAP.q30.sort.pileup.intersect.exon && echo '--step2.2 intersect done---'

python ${scripts_extract_snp_with_10more_valid_reads_py} hap1.BM-${DAP}DAP.q30.sort.pileup.intersect.exon hap1.BM-${DAP}DAP.q30.sort.pileup.intersect.exon.with10more.valid.reads

python ${scripts_ASE_snp_py} hap1.BM-${DAP}DAP.q30.sort.pileup.intersect.exon.with10more.valid.reads hap1.BM-${DAP}DAP.exon.with10more.valid.reads.count.csv

##注意，hap1 hap2是反的
sed -i '1i chr\tsite\thap1_geno\thap2_geno\treads_count\tvalid_reads_count\thap1_count\thap2_count\thap1_ratio\thap2_ratio' hap1.BM-${DAP}DAP.exon.with10more.valid.reads.count.csv 

### binom --> FDR

Rscript ${scripts_binom_fdr_R} ${PATH_to_relative_R_WD}/hap1.BM-${DAP}DAP.exon.with10more.valid.reads.count.csv ${PATH_to_relative_R_WD}/hap1.BM-${DAP}DAP.exon.with10more.valid.reads.count.fdr.csv

#final choice
grep -v 'site' hap1.BM-${DAP}DAP.exon.with10more.valid.reads.count.fdr.csv | awk -F ',' '{if($9 >= 0.9 && $12 <= 0.01) print $1,$2,$2,"hap1_type"}'  | tr ' ' '\t' > hap1.fdr.hap1_type.snp
grep -v 'site' hap1.BM-${DAP}DAP.exon.with10more.valid.reads.count.fdr.csv | awk -F ',' '{if($10 >= 0.9 && $14 <= 0.01) print $1,$2,$2,"hap2_type"}'  | tr ' ' '\t' > hap1.fdr.hap2_type.snp

#  wc -l hap1.ref.hap1_type.snp hap1.ref.hap2_type.snp

cat hap1.fdr.hap1_type.snp hap1.fdr.hap2_type.snp | sort -k1.4,1n -k2,2n > hap1.fdr.allele-specific.snp
bedtools intersect -a /data21/yjiang22/hybrid_hap_graph_application/hap1.chr.gene.bed -b hap1.fdr.allele-specific.snp -wao > hap1.fdr.genes.allele-specific-snp.wao
grep -v '\-1' hap1.fdr.genes.allele-specific-snp.wao > hap1.fdr.genes.allele-specific-snp.wao1
python ${scripts_extract_ASE_gene_with_all_same_hap_snp_py} hap1.fdr.genes.allele-specific-snp.wao1 hap1.fdr.ase.gene.hap.snp.count.csv

awk '{ if($2=="0") print $0}' hap1.fdr.ase.gene.hap.snp.count.csv > hap1.fdr.ase.gene.hap2.snp.count
awk '{ if($3=="0") print $0}' hap1.fdr.ase.gene.hap.snp.count.csv > hap1.fdr.ase.gene.hap1.snp.count

awk '{ if($2 >= 3) print $0}' hap1.fdr.ase.gene.hap1.snp.count > hap1.fdr.ase.gene.hap1.snp.count.3more
awk '{ if($3 >= 3) print $0}' hap1.fdr.ase.gene.hap2.snp.count > hap1.fdr.ase.gene.hap2.snp.count.3more

### chisq --> BH

source activate r4.2
#Rscript ${scripts_binom_fdr_R} ${PATH_to_relative_binom_fdr_R_WD}/graph.BM-${DAP}DAP.exon.with10more.valid.reads.count.csv ${PATH_to_relative_binom_fdr_R_WD}/graph.BM-${DAP}DAP.exon.with10more.valid.reads.count.fdr.csv && echo '--step4 binom.R done--' #这里的csv文件目录要相对于binom.R来写
Rscript ${scripts_chisq_BH_R} ${PATH_to_relative_R_WD}/hap1.BM-${DAP}DAP.exon.with10more.valid.reads.count.csv ${PATH_to_relative_R_WD}/hap1.BM-${DAP}DAP.exon.with10more.valid.reads.count.BH.csv && echo '--step4 chisq BH done--' #这里的csv文件目录要相对于chisq BH.R WD来写

source activate base

#final choice
grep -v 'site' hap1.BM-${DAP}DAP.exon.with10more.valid.reads.count.BH.csv | awk -F ',' '{if($9 >= 0.9 && $12 <= 0.01) print $1,$2,$2,"hap1_type"}'  | tr ' ' '\t' > hap1.BH.hap1_type.snp
grep -v 'site' hap1.BM-${DAP}DAP.exon.with10more.valid.reads.count.BH.csv | awk -F ',' '{if($10 >= 0.9 && $12 <= 0.01) print $1,$2,$2,"hap2_type"}'  | tr ' ' '\t' > hap1.BH.hap2_type.snp

#  wc -l graph.ref.hap1_type.snp graph.ref.hap2_type.snp

cat hap1.BH.hap1_type.snp hap1.BH.hap2_type.snp | sort -k1.4,1n -k2,2n > hap1.BH.allele-specific.snp
bedtools intersect -a /data21/yjiang22/hybrid_hap_graph_application/hap1.chr.gene.bed -b hap1.BH.allele-specific.snp -wao > hap1.BH.genes.allele-specific-snp.wao
grep -v '\-1' hap1.BH.genes.allele-specific-snp.wao > hap1.BH.genes.allele-specific-snp.wao1
python ${scripts_extract_ASE_gene_with_all_same_hap_snp_py} hap1.BH.genes.allele-specific-snp.wao1 hap1.BH.ase.gene.hap.snp.count.csv

awk '{ if($2=="0") print $0}' hap1.BH.ase.gene.hap.snp.count.csv > hap1.BH.ase.gene.hap2.snp.count
awk '{ if($3=="0") print $0}' hap1.BH.ase.gene.hap.snp.count.csv > hap1.BH.ase.gene.hap1.snp.count

awk '{ if($2 >= 3) print $0}' hap1.BH.ase.gene.hap1.snp.count > hap1.BH.ase.gene.hap1.snp.count.3more
awk '{ if($3 >= 3) print $0}' hap1.BH.ase.gene.hap2.snp.count > hap1.BH.ase.gene.hap2.snp.count.3more




##从bam文件提取mismatch信息

python ${scripts_extract_NM_from_bam_py} hap1.BM-${DAP}DAP.q30.sort.bam > hap1.BM-${DAP}DAP.q30.sort.bam.NM
sed -i 's/\[//' hap1.BM-${DAP}DAP.q30.sort.bam.NM
sed -i 's/\]//' hap1.BM-${DAP}DAP.q30.sort.bam.NM

cat hap1.BM-${DAP}DAP.q30.sort.bam.NM | tr ', ' '\n' | sed '/^[[:space:]]*$/d' > hap1.BM-${DAP}DAP.q30.sort.bam.NM1;mv hap1.BM-${DAP}DAP.q30.sort.bam.NM1 hap1.BM-${DAP}DAP.q30.sort.bam.NM

