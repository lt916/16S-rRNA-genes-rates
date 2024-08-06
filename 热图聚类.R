
setwd("E:\\毕业数据\\数据处理R\\16S降解\\热图聚类\\")
library(pheatmap)
library(ggplot2)
C <-read.csv("OTU降解分段表.csv", header=TRUE,row.names = 1)
C <- C_1*-1
#C <- 0.693/C_2



color_vector = c("#00008B", "#102D9B", "#215AAC", "#3288BD", "#66C2A5",  "#E6F598", "#FFFFBF", "#FED690", "#FDAE61", "#F46D43", "#D53E4F")

p <- pheatmap(mat = C, display_numbers = F, fontsize_number=5,number_color="red",cellwidth =15,cellheight = 15, fontsize_row=12, fontsize = 15,fontsize_col=12, width=30,height=30, angle_col = "45",  #angle_col = "o"调整列命至水平方向
              #fontsize_number=25表示每个块中字符的大小,cellwidth = 55，表示模块的大小。 fontsize_row=20表示X轴上的字体大小。     
              #gaps_row = breaks, ##按照行进行分块
              show_rownames=T, show_colnames=T,legend_breaks=c(),
              color = colorRampPalette(color_vector)(100),
              cluster_cols = F,cluster_rows = T ) 
p
ggsave('OTU降解速率_30.pdf', p, width = 9, height = 5)



#############################################################################
########################### 热图 ############################################
#############################################################################

library(corrplot)
library(matlab)
library(ggplot2)  
library(RColorBrewer)  
library(reshape2) 
setwd("E:/毕业数据/数据处理R/16S降解/降解拟合")
mat<-read.csv("OTU降解分段表_半衰期.csv",header = T,row=1)
mat<- t(mat)
mydata<- melt(as.matrix(mat))
mydata$AbsValue<-abs(mydata$value)
p = ggplot(mydata, aes(x= Var1 , y=Var2)) +
  geom_point(aes(size=AbsValue,fill = value), shape=21, colour="black",stroke=0.5) +
  theme(axis.text.x = element_text(color="black",size = 23,face = "plain", angle = 45,vjust = 0.5))+#X轴坐标，size 大小， face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
  theme(axis.text.y = element_text(color="black",size = 23,face = "plain"))+#Y轴坐标
  xlab("") + theme(axis.title.x = element_text(color="black",size = 20, face = "bold"))+#X轴标题
  ylab("") + theme(axis.title.y = element_text(size = 20, face = "plain"))+#Y轴标题
  scale_fill_gradientn(colours=c(brewer.pal(7,"Set1")[2],"white",brewer.pal(7,"Set1")[1]),na.value=NA,limits=c(0,1100))+
  scale_size_area(max_size=15, guide=FALSE) +
  theme(legend.text=element_text(size=16, face = "plain"))+#图例字体
  theme(
    text=element_text(size=15,face="plain",color="black"),
    axis.title=element_text(size=13,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    legend.position="right"
  )
p
ggsave('半衰期分段.pdf', p, width = 16, height = 7)






