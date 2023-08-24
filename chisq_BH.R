setwd('/data21/yjiang22/new_hybrid_hap_graph_BM')

args <- commandArgs(T)
classification.wo0.csv <- as.character(args[1])
#classification.wo0_fdr.csv <- as.character(args[2])
out_BH <- as.character(args[2])
# 加载所需的包
library(tidyverse)

# 读取文件，假设你的文件叫"data.csv"
data <- read.csv(classification.wo0.csv,header=T,sep='\t')

# 执行卡方检验并将结果添加为新列
data <- data %>% 
  mutate(chisq_test = pmap(list(hap1_count, hap2_count), function(hap1, hap2) {
    test <- chisq.test(c(hap1, hap2))
    return(test$p.value)
  }))

# 将 chisq_test 列的内容转化为数字
data$chisq_test <- as.numeric(data$chisq_test)

# 用Benjamini & Hochberg方法对p值进行校正
data <- data %>% 
  mutate(p_adj = p.adjust(chisq_test, method = "BH"))

# 写入新的CSV文件
write_csv(data, out_BH)

##usage:    >Rscript chisq_BH.R <input *count.csv> <output *count.BHcsv>
