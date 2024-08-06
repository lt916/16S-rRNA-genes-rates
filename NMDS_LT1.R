###NMDS
group <-read.csv("group_1.csv", header=F, row.names = 1)
otu<-read.csv("16S_OTU_rare_1.csv", header=TRUE, row.names=1)
otut<-t(otu)
#平方根转化
otut<-otut^0.5
library(ggplot2)
#创建颜色梯度
library(RColorBrewer)
color1 <- colorRampPalette(c("#FF0000","#F08080"))(6)
color2 <- colorRampPalette(c("#003153","#6495ED"))(6)
color3 <- colorRampPalette(c("#00FF00","#CCFF00"))(6)
color4 <- colorRampPalette(c("#136C87","#82A3AD"))(6)
color5 <- colorRampPalette(c("#0FC4B8","#86D4CF"))(6)
color6 <- colorRampPalette(c("#DE1D8D","#D6B2C7"))(6)


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

#提取样本点分组坐标
group<- sample_site$names
group<- as.matrix(group)
dev.new()
#ggplot2 作图
nmds_plot_16s <- ggplot(sample_site, aes(NMDS1, NMDS2, group = group)) +
  geom_point(aes(color = group), size = 4, alpha = 0.8) + #可在这里修改点的透明度、大小
  #, shape = group scale_shape_manual(values = c(17)) + #可在这里修改点的形状
  scale_color_manual(values = c(color1, color2,color3,color4,color5,color6)) + #可在这里修改点的颜色
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black',fill = 'transparent')) + #去掉背景框
  theme(legend.key = element_rect(fill = 'transparent'), 
        legend.title = element_blank(),#去掉图例标题及标签背背景
        legend.text=element_text(size=rel(0.8)),
        legend.spacing.x=unit(0,'cm'),
        legend.spacing.y=unit(0,'cm'),#设置图例大小间隔
        legend.position = "bottom") + 
  labs(x = '', y = '', title = paste('Stress =', round(NMDS.stress, 4))) +
  theme(plot.title = element_text(size=16,hjust = 1,vjust=-6),
        axis.title=element_text(size=16,face="bold"),
        axis.text.x = element_text(size = 14,color="black"),
        axis.text.y = element_text(size = 14,color="black"))+
  guides(col=guide_legend( nrow = 3,  byrow = TRUE, direction = "horizontal")) 
nmds_plot_16s

ggsave('NMDS_16S.pdf', nmds_plot_16s, width = 9, height = 6)
