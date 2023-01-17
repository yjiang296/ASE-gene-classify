setwd('/Users/yjiang/Documents/script/hybrid/parent-specific-classification/kqq')
#rm(list = ls(all=T))  

#args <- commandArgs(T)
#classification.wo0.csv <- as.character(args[1])
#classification.wo0_fdr.csv <- as.character(args[2])
##seedling ====
count_csv <- read.csv("./seedling/seedling.graph.count.csv",header=T,sep="\t")

primary_count <- count_csv$primary_count
alternative_count <- count_csv$alternative_count
totalreads <- primary_count + alternative_count#总reads

len <- length(primary_count)

pv <- numeric(len) #设置相同长度的变量

P <- 0.5 

for(i in 1:len){pv[i]=binom.test(primary_count[i], totalreads[i], P)$p.value} #计算P值

qv <- p.adjust(pv,"fdr") #对p值进行FDR检验
total=cbind(count_csv,pv)
total1 =cbind(total,qv) #将FDR(qv)列加入到原文件后面

write.csv(total1, file = "./seedling/seedling.graph.count.fdr.csv", row.names=FALSE, quote=FALSE) #输出合并后的文件
##usage:    >Rscript binom_fdr.R <input *count.csv> <output *count.fdr.csv>

