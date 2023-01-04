library(ggplot2)
library(ggpubr)
rm(list = ls(all=T))  ##清除变量


## endosperm ====
data <- read.csv('./endosperm/endosperm.hap1_classification.reads20.fdr001.labeled.sig90per.and.nosig.csv',sep=',',header=TRUE)

p1 <- ggplot(data, aes(x=log10(hap1_count), y=log10(hap2_count), group=label)) +
  geom_point(aes(color=label),size=0.8,pch=16) + 
  scale_color_manual(values=c('#F08080','#1E90FF', '#A9A9A9')) + 
  coord_fixed(ratio=1) + 
  ylim(0,5) +
  xlim(0,5) +
  theme_light()
  
p1 +labs(title="endosperm hap1 fdr0.01")

data <- read.csv('./endosperm/endosperm.graph_classification.reads20.fdr001.labeled.sig90per.and.nosig.csv',sep=',',header=TRUE)

p2 <- ggplot(data, aes(x=log10(hap1_count), y=log10(hap2_count), group=label)) +
  geom_point(aes(color=label),size=0.8,pch=16) + 
  scale_color_manual(values=c('#F08080','#1E90FF', '#A9A9A9')) + 
  coord_fixed(ratio=1) + 
  ylim(0,5) +
  xlim(0,5) +
  theme_light()

p2 +labs(title="endosperm graph fdr0.01 ")

ggarrange(p1,p2,
          labels=c("endosperm hap1 fdr0.01",
                   "endosperm graph fdr0.01"),
          ncol=2,nrow=1,
          font.label = list(size = 10, color = "black"))




## embryo ====
data <- read.csv('./embryo/embryo.hap1_classification_fdr001.labeled.reads20.sig90per.and.nosig.csv',sep=',',header=TRUE)

p1 <- ggplot(data, aes(x=log10(hap1_count), y=log10(hap2_count), group=label)) +
  geom_point(aes(color=label),size=0.8,pch=16) + 
  scale_color_manual(values=c('#F08080','#1E90FF', '#A9A9A9')) + 
  coord_fixed(ratio=1) + 
  ylim(0,4.5) +
  xlim(0,4.5) +
  theme_light()

p1 +labs(title="embryo hap1 fdr0.01")

data <- read.csv('./embryo/embryo.graph_classification_fdr001.labeled.reads20.sig90per.and.nosig.csv',sep=',',header=TRUE)

p2 <- ggplot(data, aes(x=log10(hap1_count), y=log10(hap2_count), group=label)) +
  geom_point(aes(color=label),size=0.8,pch=16) + 
  scale_color_manual(values=c('#F08080','#1E90FF', '#A9A9A9')) + 
  coord_fixed(ratio=1) + 
  ylim(0,4.5) +
  xlim(0,4.5) +
  theme_light()

p2 +labs(title="embryo graph fdr0.01 ")

ggarrange(p1,p2,
          labels=c("embryo hap1 fdr0.01",
                   "embryo graph fdr0.01"),
          ncol=2,nrow=1,
          font.label = list(size = 10, color = "black"))


## seedling ====
data <- read.csv('./seedling/seedling.hap1_classification_fdr001.labeled.reads20.sig90per.and.nosig.csv',sep=',',header=TRUE)
#ggplot(data, aes(x=log10(hap1_count), y=log10(hap2_count), color=log10(seed_graph))) +
#       geom_point(pch=16,size=1) +
#       scale_color_gradient(low="blue", high="red")
p1 <- ggplot(data, aes(x=log10(hap1_count), y=log10(hap2_count), group=label)) +
  geom_point(aes(color=label),size=0.8,pch=16) + 
  scale_color_manual(values=c('#F08080','#1E90FF', '#A9A9A9')) + 
  coord_fixed(ratio=1) + 
  ylim(0,4.4) +
  xlim(0,4.4) +
  theme_light()

p1 +labs(title="seedling hap1 fdr0.01")

data <- read.csv('./seedling/seedling.graph_classification_fdr001.labeled.reads20.sig90per.and.nosig.csv',sep=',',header=TRUE)
#ggplot(data, aes(x=log10(hap1_count), y=log10(hap2_count), color=log10(seed_graph))) +
#       geom_point(pch=16,size=1) +
#       scale_color_gradient(low="blue", high="red")
p2 <- ggplot(data, aes(x=log10(hap1_count), y=log10(hap2_count), group=label)) +
  geom_point(aes(color=label),size=0.8,pch=16) + 
  scale_color_manual(values=c('#F08080','#1E90FF', '#A9A9A9')) + 
  coord_fixed(ratio=1) + 
  ylim(0,4.4) +
  xlim(0,4.4) +
  theme_light()

p2 +labs(title="seedling graph fdr0.01 ")

ggarrange(p1,p2,
          labels=c("seedling hap1 fdr0.01",
                   "seedling graph fdr0.01"),
          ncol=2,nrow=1,
          font.label = list(size = 10, color = "black"))





