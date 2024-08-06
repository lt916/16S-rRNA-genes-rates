
library(vegan)
library(ggplot2)
#读入数据文件
setwd("E:\\毕业数据\\数据处理R\\16S降解\\柱状图\\") 
databar=read.csv(file='ITS降解速率.csv',header = T,stringsAsFactors = F)
x = paste0('S',1:30)#定义X 的顺序
databar$Sample = factor(databar$Sample,levels = x)#group按照定义的X进行排序



q1<-ggplot(data=databar, mapping=aes(x = Sample, y = TIME,fill=group))+
  geom_bar(stat="identity",position=position_dodge(0.8),width=0.8)+#组内间隔position_dodge(0.7),矩形条的宽度width=0.5达到的效果
  #scale_fill_brewer(palette="Set3")#修改颜色
  theme_bw()+#去背景
  theme(panel.grid=element_blank())+#去网格线
  scale_fill_manual(values =c( "#00AFBB" ))+#颜色的十六进制代码，或直接用red、blue、green等也可,"fc744c"红, "#00AFBB"蓝
  theme(axis.text.x = element_text(color="black",size = 10,face = "plain",angle = 45))+#X轴坐标，size 大小， face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
  theme(axis.text.y = element_text(color="black",size = 14,face = "plain"))+#Y轴坐标
  xlab("Site") + theme(axis.title.x = element_text(size = 18, face = "plain"))+#X轴标题
  ylab("The turnover time for extracellular DNA (lg)") + theme(axis.title.y = element_text(size = 18, face = "plain"))+#Y轴标题
  theme(legend.text=element_text(size=14, face = "plain"))+#图例字体
  guides(fill=guide_legend(title=NULL))#移除图例标题
  #theme(legend.position="bottom")+#左边left,右边 right, 底部bottom
q1
ggsave('ITS半衰期.pdf', q1, width = 8, height = 6)
