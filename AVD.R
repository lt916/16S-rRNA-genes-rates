
setwd("E:/毕业数据/数据处理R/16S降解/群落稳定性")
#读取丰度数据
otu <- read.csv("降解16S_OTU_rare.csv",row.names = 1,header = T)
otu <- otu/52826
#计算变异度
ai <- abs(otu-apply(otu, 1, mean))/apply(otu, 1, sd)
ai1 <- na.omit(ai)#删除空值
#AVD <- apply(ai,2,sum)/(1*nrow(otu))
AVD <- colSums(ai1)/(1*nrow(otu))
#head(AVD)

 #读取分组依据
Group <- read.csv("group.csv",row.names = 1,header = T)
#合并
Group$AVD <- AVD
write.csv(Group,"AVD.csv")
#Group



#绘制箱线图
library(ggplot2)
#设置新罗马字体
#windowsFonts(A=windowsFont("Times New Roman"),
             #B=windowsFont("Arial"))
# Default plot

# Change fill color by group (dose))
p <- ggplot(Group, aes(x = Group, y = AVD)) + 
  geom_point(size = 0.5)+
  geom_boxplot(aes(fill = Group),alpha=0.5) +
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#F39B7FFF","#91D1C2FF", "#8491B4FF"))+
  geom_jitter(aes(colour=Group),
          position = position_jitter(0.25),size=2,alpha=0.8)+
  scale_colour_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#F39B7FFF","#91D1C2FF", "#8491B4FF"))+
  theme_bw()+scale_y_continuous(limits = c(0,0.7),expand = c(0,0))

p
print(p)
ggsave('AVD1.pdf', p, width = 7, height = 5)


###差异显著性检验
Group1 <- read.csv("AVD.csv",row.names = 1,header = T)
Group1$Group<- as.factor(Group1$Group)

###正态分布检验
shapiro <- tapply(Group1$Group, Group1$AVD, shapiro.test)
shapiro
#差性验方差齐性检验
bartlett.test(Group1$AVD ~ Group1$Group, data=Group1)
####重复测量方差分析
oneway<-aov(Group1$AVD ~ Group1$Group,data = Group1)
anova(oneway)


























#绘制箱线图

library(ggplot2)
library(ggpubr)
library(ggsci)
#设置新罗马字体
windowsFonts(A=windowsFont("Times New Roman"),
             B=windowsFont("Arial"))
# Default plot
e <- ggplot(Group, aes(x = Group, y = 1-AVD))
# Change fill color by group (dose))
P <- ggplot(Group, aes(x = Group, y = 1-AVD)) +
  geom_boxplot(aes(fill = Group),alpha=0.5) +
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#F39B7FFF","#91D1C2FF", "#8491B4FF"))+#定义填充颜色
  geom_line(aes(group=group) ,colour="#9C9C9C",lwd=0.2)+
  geom_point(size = 0.5)+
  scale_y_continuous(name="Brary-cris",labels=scaleFUN)+###y轴标题
  scale_x_discrete(name="type")+
  theme_bw()+
  theme(legend.position='none')+
  theme(panel.background = element_rect(fill = "white", colour = "black"), #外框设置黑色
        panel.grid = element_blank(),
        panel.border = element_blank(),
        text = element_text(colour="black",size=14),
        axis.title = element_text(colour="black",face="bold"),
        axis.text=element_text(colour="black",size=14),#坐标轴刻度
        axis.title.x = element_blank(),#坐标轴标题
        axis.title.y = element_text(colour="black", size = 12))+
  xlab("")+ylab("Bray-Curtis similarity of prokaryotic community")+ylim(0.1,0.8)+
  
  scale_x_discrete(labels = c("D1","D3","D6","D12","D24","D48"))+ ##重新设置X轴刻度上的名字
  #stat_compare_means(comparisons=my_comparisons,label ="p.signif",step_increase = 0.1,map_signif_level = T,test = wilcox.test,paired = TRUE)+ # Add pairwise 
  stat_compare_means(label.y = 0.75,label.x = 3,size=6) 

  
  
  
  
  geom_jitter(aes(colour=Group),
              position = position_jitter(0.25),size=5,alpha=0.8)+
  scale_colour_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#F39B7FFF","#91D1C2FF", "#8491B4FF"))+
  theme_bw()+scale_y_continuous(limits = c(-1,1),expand = c(0,0))+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )+
  theme(text=element_text(family="A",size=20))

P




ggboxplot(group, x = 'group', y = 'AVD', fill = 'Group1', 
          color = 'gray30', width = 0.6, size = 1, legend = 'right') +
  scale_fill_manual(values = c('#E7B800', '#00AFBB')) +
  labs(x = '', y = 'AVD', fill = '')
