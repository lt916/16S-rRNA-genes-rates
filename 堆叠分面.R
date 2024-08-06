library(reshape2)
library(ggplot2)
setwd("E:\\毕业数据\\数据处理R\\16S降解\\堆叠图\\")
phylum_top10 <- read.csv('16S降解-堆叠.csv', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
#整理成 ggplot2 作图格式
phylum_top10$Taxonomy <- factor(rownames(phylum_top10), levels = rev(rownames(phylum_top10)))

phylum_top10 <- melt(phylum_top10, id = 'Taxonomy')
#添加分组，这次我们根据样本分组绘制分面
group <- read.csv('group.csv', row.names = 1, stringsAsFactors = FALSE)
names(group)[1] <- 'variable'
phylum_top10 <- merge(phylum_top10, group, by = 'variable')
x = paste0('Site',1:30)#定义X 的顺序
phylum_top10$group = factor(phylum_top10$group,levels = x)#group按照定义的X进行排序


#绘制带分面的柱状图
p <- ggplot(phylum_top10, aes(variable, 100 * value, fill = Taxonomy)) +
  geom_col(position = 'stack', width = 0.6) +
  facet_wrap(~group, scales = 'free_x', ncol = 6) +# 按group组，X轴，分6面
  scale_fill_manual(values =  rev(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray'))) +
  theme_bw() +theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(strip.background = element_rect( size=0.5))+###site1 大小
  theme(strip.text.x = element_text(size = 10, colour = "black")) + # 设置分面的字字体大小
  theme(axis.text = element_blank())+   ## 删去所有刻度标签
  #theme(axis.text.x = element_text(color="black",size = 10,face = "plain"))+#X轴坐标，size 大小， face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
  theme(axis.text.y = element_text(color="black",size = 10,face = "plain"))+#Y轴坐标
  xlab("") + theme(axis.title.x = element_text(size = 8, face = "plain"))+#X轴标题
  ylab("Relative abundance(%)") + theme(axis.title.y = element_text(size = 15, face = "plain"))+#Y轴标题
  theme(legend.text=element_text(size=12, face = "plain"))+#图例字体
  #theme(legend.position="bottom")+#左边left,右边 right, 底部bottom
  guides(fill=guide_legend(title=NULL))+#移除图例标题
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 8))+ 
  theme(panel.border = element_rect(color="black", size=0.1))+#边框加粗
  theme(axis.ticks.x=element_line(color="black",size=0.3,lineend = 1))+#x轴坐标刻度
  theme(axis.ticks.y=element_line(color="black",size=0.3,lineend = 1))#y轴坐标刻度
  #guides(fill=FALSE)#移除图例
  p
  
ggsave('DD.pdf', p, width = 10, height = 7.5)
  
