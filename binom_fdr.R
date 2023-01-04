setwd('/Users/yjiang/Documents/script/杂种胚/parent-specific-classification/kqq')
rm(list = ls(all=T))  ##清除变量

#args <- commandArgs(T)
#classification.wo0.csv <- as.character(args[1])
#classification.wo0_fdr.csv <- as.character(args[2])
##seedling ====
count_csv <- read.csv("./seedling/seed.graph_classification.wo0.csv",header=T,sep="\t")

#count_csv <- read.csv(classification.wo0.csv,header=T,sep=",")

hap1_count <- count_csv$hap1_count
hap2_count <- count_csv$hap2_count
totalreads <- hap1_count + hap2_count#总reads

len <- length(hap1_count)

pv <- numeric(len) #设置相同长度的变量

P <- 0.5 

for(i in 1:len){pv[i]=binom.test(hap1_count[i], totalreads[i], P)$p.value} #计算P值

qv <- p.adjust(pv,"fdr") #对p值进行FDR检验
total=cbind(count_csv,pv)
total1 =cbind(total,qv) #将FDR(qv)列加入到原文件后面

write.csv(total1, file = "./seedling/seed.graph_classification.wo0_fdr.csv", row.names=FALSE, quote=FALSE) #输出合并后的文件
#write.csv(total, file = classification.wo0_fdr.csv, row.names=FALSE, quote=FALSE) #输出合并后的文件
##usage:    >Rscript binom_fdr.R <input_classification.wo0.csv> <output_classification.wo0_fdr.csv>

##endosperm ====
setwd('/Users/yjiang/Documents/script/杂种胚/parent-specific-classification/kqq')
rm(list = ls(all=T))  ##清除变量

#args <- commandArgs(T)
#classification.wo0.csv <- as.character(args[1])
#classification.wo0_fdr.csv <- as.character(args[2])
count_csv <- read.csv("./endosperm/endosperm.hap1_classification.reads20.csv",header=T,sep="\t")

#count_csv <- read.csv(classification.wo0.csv,header=T,sep=",")

hap1_count <- count_csv$hap1_count
hap2_count <- count_csv$hap2_count
totalreads <- hap1_count + hap2_count#总reads

len <- length(hap1_count)

pv <- numeric(len) #设置相同长度的变量

P <- 0.667

for(i in 1:len)
  {
  pv[i]=binom.test(max(c(hap1_count[i],hap2_count[i])), totalreads[i], P)$p.value
  } #计算P值

qv <- p.adjust(pv,"fdr") #对p值进行FDR检验
total=cbind(count_csv,pv)
total1 =cbind(total,qv) #将FDR(qv)列加入到原文件后面

write.csv(total1, file = "./endosperm/endosperm.hap1_classification.reads20.fdr.csv", row.names=FALSE, quote=FALSE) #输出合并后的文件
#write.csv(total, file = classification.wo0_fdr.csv, row.names=FALSE, quote=FALSE) #输出合并后的文件
##usage:    >Rscript binom_fdr.R <input_classification.wo0.csv> <output_classification.wo0_fdr.csv>