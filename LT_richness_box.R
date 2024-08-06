#richness
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggsci)
setwd("E:\\毕业数据\\数据处理R\\箱线图\\")
species<-read.csv("降解16S_richness.csv",header=T,row.names = 1)
species$type<-as.factor(species$type)
species$ID<-as.factor(species$ID)

###正态分布检验
shapiro <- tapply(species$richness, species$type, shapiro.test)
shapiro
#差性验方差齐性检验
bartlett.test(richness ~ type, data=species)
####重复测量方差分析
library(rstatix)
res.aov<-anova_test(data = species, dv = richness,   wid= ID, within =  type)
get_anova_table(res.aov)
######################################################################
#####################################################################
scaleFUN<-  function(x) sprintf ("%.1f",x)
#colnames(species)<-c("Row.names","species","type")
species$group<-c(1:30,1:30,1:30,1:30,1:30,1:30)###根据每组样品数调整1:30

p4 = ggplot(species,aes(fill = type,x=factor(type,level=c("D1","D2","D3","D4","D5","D6")), y=richness/1000))+###abundance： y=species
  stat_boxplot(geom="errorbar",width=0.15)+
  geom_boxplot(alpha=0.7,aes(fill=type))+
  geom_line(aes(group=group) ,colour="#9C9C9C",lwd=0.2)+
  geom_point(size = 0.5)+
  scale_y_continuous(name="Prokaryotic richness(103)",labels=scaleFUN)+###y轴标题
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
        axis.title.y = element_text(colour="black", size = 14))+
        xlab("")+ylab("Prokaryotic richness(103)")+ylim(2,12)+
   scale_x_discrete(labels = c("D1","D3","D6","D12","D24","D48"))+ ##重新设置X轴刻度上的名字
  #stat_compare_means(comparisons=my_comparisons,label ="p.signif",step_increase = 0.1,map_signif_level = T,test = wilcox.test,paired = TRUE)+ # Add pairwise 
   stat_compare_means(label.y = 12,label.x = 1.5,size=6) # Add global p-value

p4<-p4+scale_fill_brewer(palette = "Blues",direction = -1) ##palette = "Blues"设置蓝色渐变色，且由深边浅；direction = -1表示颜色换顺序
p4
ggsave('16S_richness_箱线图.pdf', p4, width = 7, height = 5)




##################################################################################
###拷贝数
species<-read.csv("ITS_拷贝数.csv",header=T,row.names = 1)
species$type<-as.factor(species$type)
species$ID<-as.factor(species$ID)
###正态分布检验
shapiro <- tapply(species$species, species$type, shapiro.test)
shapiro
#差性验方差齐性检验
bartlett.test(species ~ type, data=species)
####重复测量方差分析
library(rstatix)
res.aov<-anova_test(data = species, dv = species,   wid= ID, within =  type)
get_anova_table(res.aov)


scaleFUN<-  function(x) sprintf ("%.1f",x)
#colnames(species)<-c("Row.names","species","type")
species$group<-c(1:30,1:30,1:30,1:30,1:30,1:30)###根据每组样品数调整1:30

p4 = ggplot(species,aes(fill = type,x=factor(type,level=c("D1","D2","D3","D4","D5","D6")), y=species))+###abundance： y=species
  stat_boxplot(geom="errorbar",width=0.15)+
  geom_boxplot(alpha=0.7,aes(fill=type))+
  geom_line(aes(group=group) ,colour="#9C9C9C",lwd=0.2)+
  geom_point(size = 0.5)+
  scale_y_continuous(name="Fungal copies(lg)",labels=scaleFUN)+###y轴标题
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
        axis.title.y = element_text(colour="black", size = 14))+
  xlab("")+ylab("Fungal copies(lg)")+ylim(4.5,10)+
  scale_x_discrete(labels = c("D1","D3","D6","D12","D24","D48"))+ ##重新设置X轴刻度上的名字
  #stat_compare_means(comparisons=my_comparisons,label ="p.signif",step_increase = 0.1,map_signif_level = T,test = wilcox.test,paired = TRUE)+ # Add pairwise 
  stat_compare_means(label.y = 9.5,label.x = 4,size=6) # Add global p-value

p4<-p4+scale_fill_brewer(palette = "Blues",direction = -1) ##palette = "Blues"设置蓝色渐变色，且由深边浅；direction = -1表示颜色换顺序
p4
ggsave('ITS_copy_箱线图1.pdf', p4, width = 6, height = 6)


############################################################################
#BC距离
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggsci)
setwd("E:\\毕业数据\\数据处理R\\16S降解\\NMDS\\BC距离\\")
species<-read.csv("相似性.csv",header=T,row.names = 1)
species$type<-as.factor(species$type)
species$ID<-as.factor(species$ID)

###正态分布检验
shapiro <- tapply(species$value, species$type, shapiro.test)
shapiro
#差性验方差齐性检验
bartlett.test(value ~ type, data=species)
####重复测量方差分析
library(rstatix)
res.aov<-anova_test(data =species , dv =value ,   wid= ID, within =  type)
get_anova_table(res.aov)
######################################################################
#####################################################################
scaleFUN<-  function(x) sprintf ("%.1f",x)
#colnames(species)<-c("Row.names","species","type")
species$group<-c(1:30,1:30,1:30,1:30,1:30)###根据每组样品数调整1:30

p4 = ggplot(species,aes(fill = type,x=factor(type,level=c("group1","group2","group3","group4","group5")), y=value))+###abundance： y=species
  stat_boxplot(geom="errorbar",width=0.15)+
  geom_boxplot(alpha=0.7,aes(fill=type))+
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
  
  scale_x_discrete(labels = c("D1/D3","D1/D6","D1/D12","D1/D24","D1/D48"))+ ##重新设置X轴刻度上的名字
  #stat_compare_means(comparisons=my_comparisons,label ="p.signif",step_increase = 0.1,map_signif_level = T,test = wilcox.test,paired = TRUE)+ # Add pairwise 
  stat_compare_means(label.y = 0.75,label.x = 3,size=6) # Add global p-value

p4<-p4+scale_fill_brewer(palette = "Blues",direction = -1) ##palette = "Blues"设置蓝色渐变色，且由深边浅；direction = -1表示颜色换顺序
p4
ggsave('BC距离16S.pdf', p4, width = 7, height = 5)
