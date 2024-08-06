###NMDS
setwd("E:/毕业数据/数据处理R/16S降解/NMDS")
group <-read.csv("分组.csv", header=T,row.names=1)
otu<-read.csv("降解16S_OTU_rare.csv", header=TRUE, row.names=1)
otut<-t(otu)
#平方根转化
otut<-otut^0.5
library(vegan)
NMDS=metaMDS(otut,autotransform =FALSE)
#提取物种（OTU）排序坐标
NMDS.species <- data.frame(NMDS$species)
#提取样本排序坐标
NMDS.point <- data.frame(NMDS$point)
#可选择将结果输出并保存在本地，例如将样本坐标输出为 csv 格式
#write.csv(NMDS.point, 'nmds_16s.csv')
#提取应力函数值（stress）
NMDS.stress <- NMDS$stress
#提取样本点坐标（前两轴）
sample_site <- NMDS.point[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
#为样本点坐标添加分组信息
sample_site <- merge(sample_site, group, by = 'row.names', all.x = TRUE)


#在sample_site的基础上加group2,分颜色。
p <- ggplot(sample_site, aes(NMDS1, NMDS2, group = group)) +
  geom_point(aes(color = group), size = 5, alpha =1) + #可在这里修改点的透明度、大小
  #scale_shape_manual(values = c(15)) + #可在这里修改点的形状
  #scale_color_manual(values = c("#80b1d3","#8dd3c7","#B4ADE3","#FFA45B","#F76248","blue"))+ #可在这里修改点的颜色
  scale_color_brewer(palette="Blues", direction = -1)+#修改颜色 ，direction = -1 或者1由深到浅，改变颜色顺序
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) + #去掉背景框
  theme(legend.key = element_rect(fill = 'transparent'), legend.title = element_blank()) + #去掉图例标题及标签背景
  labs(x = 'NMDS 1', y = 'NMDS 2', title = paste('Stress =', round(NMDS$stress,3))) +
  theme(plot.title = element_text(hjust = 0.9,size = 20,face = "plain",vjust = -7))+  #标题居中
  theme(axis.text.x = element_text(size = 20,face = "plain"))+#X轴坐标，size 大小， face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
  theme(axis.text.y = element_text(size = 20,face = "plain"))+#Y轴坐标
  theme(axis.title.x = element_text(size = 20, face = "plain"))+#X轴标题
  theme(axis.title.y = element_text(size = 20, face = "plain"))+#Y轴标题
  theme(legend.text=element_text(size=20, face = "plain"))+#图例字体
  theme(panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+#边框加粗
  theme(axis.ticks.x=element_line(color="black",size=1.1,lineend = 1))+#x轴坐标刻度
  theme(axis.ticks.y=element_line(color="black",size=1.1,lineend = 1))#y轴坐标刻度
p 
ggsave('NMDS_16S_降解.pdf',p, width = 6.5, height = 5)





##################################################################################
#######################################################################################
######################################################################################
######渐变颜色
#########作图
library(ggplot2)
#创建颜色梯度
library(RColorBrewer)
color1 <- colorRampPalette(c("#EB1D9C","#EBABD3"))(6)
color2 <- colorRampPalette(c("#B430D1","#DDB5E6"))(6)
color3 <- colorRampPalette(c("#006CF0","#A4C9F5"))(6)
color4 <- colorRampPalette(c("#0BC5D6","#BEE5E8"))(6)
color5 <- colorRampPalette(c("#0BB394","#BDF0E7"))(6)
color6 <- colorRampPalette(c("#1AC963","#AAF7CA"))(6)
color7 <- colorRampPalette(c("#09B312","#A3D9A5"))(6)

color8 <- colorRampPalette(c("#EB1D9C","#EBABD3"))(6)
color9 <- colorRampPalette(c("#B430D1","#DDB5E6"))(6)
color10 <- colorRampPalette(c("#006CF0","#A4C9F5"))(6)
color11 <- colorRampPalette(c("#0BC5D6","#BEE5E8"))(6)

color12 <- colorRampPalette(c("#EB2F2F","#FCC0C0"))(6)
color13 <- colorRampPalette(c("#B430D1","#EFC0FA"))(6)
color14 <- colorRampPalette(c("#00A1F7","#B3E2FC"))(6)
color15 <- colorRampPalette(c("#05A88A","#B3E6DD"))(6)
color16 <- colorRampPalette(c("#D6D600","#E3E3AA"))(6)
color17 <- colorRampPalette(c("#DB2194","#F0AFD7"))(6)
color18 <- colorRampPalette(c("#F07800","#F7D1AA"))(6)
color19 <- colorRampPalette(c("#0036D6","#D1C1F7"))(6)
color20 <- colorRampPalette(c("#02ADBD","#BCF0F5"))(6)
color21 <- colorRampPalette(c("#00DE5D","#C6F2D9"))(6)

color22 <- colorRampPalette(c("#EB1D9C","#EBABD3"))(6)
color23 <- colorRampPalette(c("#B430D1","#DDB5E6"))(6)
color24 <- colorRampPalette(c("#006CF0","#A4C9F5"))(6)
color25 <- colorRampPalette(c("#0BC5D6","#BEE5E8"))(6)
color26 <- colorRampPalette(c("#0BB394","#BDF0E7"))(6)
color27 <- colorRampPalette(c("#1AC963","#AAF7CA"))(6)
color28 <- colorRampPalette(c("#09B312","#A3D9A5"))(6)
color29 <- colorRampPalette(c("#E0E316","#F6F7B7"))(6)
color30 <- colorRampPalette(c("#EB710E","#F2D5A7"))(6)
#颜色分组按照生态系统划分

#提取样本点分组坐标
sample_site<-read.csv("sample_site.csv", header=TRUE, row.names=1)

#dev.new()
#ggplot2 作图
#library(ggnewscale)
nmds_plot_16s <- ggplot(sample_site, aes(NMDS1, NMDS2, group = Group1)) +
  geom_text(aes(label = sample_site$names)) +##加标签
  #geom_point(aes(color = factor(group_2)),size = 5, alpha = 0.6) + #可在这里修改点的透明度、大小
  scale_color_manual(values = c(color1,color2,color3,color4,color5,color6))+
  #scale_shape_manual(values = c(17,18,19,15)) + #可在这里修改点的形状
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black',fill = 'transparent')) + #去掉背景框
  theme(legend.key = element_rect(fill = 'transparent'), 
        legend.position = "none") + #去除图例，position = "bottom"表示图例在底部
  theme(plot.title = element_text(size=16,hjust = 1,vjust=-6),
        axis.title=element_text(size=16,face="bold"),
        axis.text.x = element_text(size = 14,color="black"),
        axis.text.y = element_text(size = 14,color="black"))
  #guides(col=guide_legend(nrow = 3,  byrow = TRUE, direction = "horizontal")) 
nmds_plot_16s


ggsave('NMDS_16S(1).pdf', nmds_plot_16s, width = 9, height = 8)
