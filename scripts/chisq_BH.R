setwd('/data21/yjiang22/new_hybrid_hap_graph_BM')

args <- commandArgs(T)
classification.wo0.csv <- as.character(args[1])
out_BH <- as.character(args[2])

# load package
library(tidyverse)

data <- read.csv(classification.wo0.csv,header=T,sep='\t')

# Perform a chi-square test and add the result as a new column
data <- data %>% 
  mutate(chisq_test = pmap(list(hap1_count, hap2_count), function(hap1, hap2) {
    test <- chisq.test(c(hap1, hap2))
    return(test$p.value)
  }))

# change chisq_test colunm data type to numerical
data$chisq_test <- as.numeric(data$chisq_test)

# use Benjamini & Hochberg method to adjust p value
data <- data %>% 
  mutate(p_adj = p.adjust(chisq_test, method = "BH"))

# write into a new CSV
write_csv(data, out_BH)

##usage:    >Rscript chisq_BH.R <input *count.csv> <output *count.BH.csv>
