setwd("E:\\毕业数据\\PMA_QPCR\\")
####α多样性\
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggsci)
species <-read.csv("16S_PMA_QPCR箱线.csv", header=TRUE, row.names = 1)

scaleFUN<-  function(x) sprintf ("%.1f",x)
my_comparisons<-list(c("TDNA","IDNA"))

#dev.new()
p4 = ggplot(species,aes(x=factor(type,level=c("TDNA","IDNA")), y=copies,color = type))+
  stat_boxplot(geom="errorbar",width=0.15)+
  geom_boxplot(alpha=0.9,aes(fill=type))+
  geom_line(aes(group=group) ,colour="grey70",lwd=0.2)+
  geom_point()+
  scale_y_continuous(name= expression(ARGs ~ copies ~(10^9 ~ per~ g ~dry ~soil)),labels=scaleFUN)+
  scale_x_discrete(name=" ")+
  scale_color_manual(values = c("#2958a2","#c6472e"))+  ##这个表示箱线图上点的颜色
  scale_fill_manual(values = c("#2958a2","#c6472e"))+ ##这个表示箱线图的颜色
  theme_bw()+
  theme(legend.position='none')+
  theme(panel.background = element_rect(fill = "white", colour = "black"), #外框设置黑色
        panel.grid = element_blank(),
        panel.border = element_blank(),
        text = element_text(colour="black",size=20),
        axis.title = element_text(colour="black",color="black"),
        axis.text=element_text(colour="black",size=20),#坐标轴刻度
        axis.title.x = element_blank(),#坐标轴标题
        axis.title.y = element_text(colour="black", size = 20))+#设置y轴标题大小
        ylab("16S rRNA gene copies (lg g-1 soil)") +
  theme(axis.text.x = element_text(size = 20,color="black"))+ ##设置X轴刻度上的字体大小
  theme(axis.text.y  = element_text(size = 20,color="black"))+##设置y轴刻度上的字体大小
  theme(axis.title.x  = element_text(size = 20, color="black"))+#设置x轴标题大小
  #theme(axis.title.y  = element_text(size = 18, color="black"))+#设置y轴标题大小
  ylim(9.2,12.3)+ #设置y轴范围
  stat_compare_means(comparisons=my_comparisons,label ="p.signif",size=8, step_increase = 0.1,map_signif_level = T, method = "t.test", paired = TRUE)+ # Add pairwise 
  stat_compare_means(label.y = 12,label.x = 0.6,size=6,method = "t.test",paired = TRUE) # Add global p-value,在这里添加method选择需要进行参数检验的方法，并显示显著性标记到图形上。
#label.y = 50,和label.x = 1.5 是控制图例在图上的位置
#method默认的参数检验方法为"wilcox.test"
p4

#p4<-p4+scale_fill_lancet()

ggsave('16S_PMA_QPCR.pdf', p4, width = 6, height = 6)

