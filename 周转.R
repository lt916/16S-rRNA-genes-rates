
#####################################################################################
##############################          ###############################################
#############################  周转  ##############################################
#############################           #################################################
#####################################################################################
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/速率-胞外")
RATE<-read.csv("16S降解速率表_29.csv", header=TRUE, row.names=1)
IDNA<-read.csv("IDNA X COPY.csv", header=TRUE, row.names=1)
EDNA<-read.csv("EDNA.csv", header=TRUE, row.names=1)
data1 <- RATE*EDNA*-1
data2 <- IDNA/data1
data2[sapply( data2, is.infinite )] <- 0  ##Inf变为0
data2[is.na(data2)] = 0  ###NA值变为0
write.csv(data2, file = "16S降解周转.csv")


library(corrplot)
library(matlab)
library(ggplot2)  
library(RColorBrewer)  
library(reshape2) 
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/速率-胞外")
mat<-read.csv("周转分段.csv",header = T,row=1)
mat<- t(mat)
mydata<- melt(as.matrix(mat))
mydata$AbsValue<-abs(mydata$value)
p = ggplot(mydata, aes(x= Var1 , y=Var2)) +
  geom_point(aes(size=AbsValue,fill = value), shape=21, colour="black",stroke=0.5) +
  theme(axis.text.x = element_text(color="black",size = 23,face = "plain", angle = 45,vjust = 0.5))+#X轴坐标，size 大小， face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
  theme(axis.text.y = element_text(color="black",size = 23,face = "plain"))+#Y轴坐标
  xlab("") + theme(axis.title.x = element_text(color="black",size = 20, face = "bold"))+#X轴标题
  ylab("") + theme(axis.title.y = element_text(size = 20, face = "plain"))+#Y轴标题
  scale_fill_gradientn(colours=c(brewer.pal(7,"Set1")[2],"white",brewer.pal(7,"Set1")[1]),na.value=NA)+
  scale_size_area(max_size=15, guide=FALSE) +
  theme(legend.text=element_text(size=16, face = "plain"))+#图例字体
  theme(
    text=element_text(size=15,face="plain",color="black"),
    axis.title=element_text(size=13,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    legend.position="right"
  )
p
ggsave('周转时间分段.pdf', p, width = 17, height = 10)


###################################################################
####################            ###################################
####################周转分类求和###################################
####################            ###################################
###################################################################
setwd("D:\\研究生\\李婷-DNA降解\\16S降解\\180+TP结果\\180+TP结果\\速率-胞外\\分类水平")
RATE<-read.csv("Phylum降解速率表29.csv", header=TRUE, row.names=1)
IDNA<-read.csv("Phylum_IDNA_copy.csv", header=TRUE, row.names=1)
EDNA<-read.csv("Phylum_EDNA.csv", header=TRUE, row.names=1)
data1 <- RATE*EDNA*-1
data2 <- IDNA/data1
data2[sapply( data2, is.infinite )] <- 0  ##Inf变为0
data2[is.na(data2)] = 0  ###NA值变为0
write.csv(data2, file = "phylum降解周转.csv")