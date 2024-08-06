#调用R包
library(ggplot2)
library(tidyverse)
library(randomForest)
library(rfUtilities)
library(rfPermute)
#读取数据
setwd("E:/毕业数据/数据处理R/16S降解/随机森林")
mydata <- read.csv("RF_ITS.csv",header = T)
#构建模型查看模型分类准确性
set.seed(123)
treat_rf <- randomForest(rate ~ ., data= mydata,importance=TRUE,proximity=TRUE)
treat_rf
#通过置换检验进一步查看模型的显著性水平
set.seed(123)
treat_perm <- rf.significance(treat_rf, mydata, nperm=99, ntree=500)
treat_perm
#依次检验每个变量的重要性并作图
set.seed(123)
ANPP_rfP<- rfPermute(rate ~ ., data = mydata, ntree = 500,
                     na.action = na.omit, nrep = 100,num.cores = 1)
ANPP_dat <-  importance(ANPP_rfP, sort.by = NULL, decreasing = TRUE)
p <- ANPP_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>%
  ggplot(aes(x = names, y = X.IncMSE))+
  geom_bar(aes(fill = group),stat = "identity")+
  geom_text(aes(y = X.IncMSE + 1,label = label))+
  labs(x = "", y = "%IncMSE")+
  theme(axis.text.x = element_text(color="black",size = 11,face = "plain", angle = 0,vjust = 0.5))+#X轴坐标，size 大小， face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
  theme(axis.text.y = element_text(color="black",size = 11,face = "plain"))+
  coord_flip()
p
ggsave('RF_ITS.pdf', p, width = 9, height = 7)
