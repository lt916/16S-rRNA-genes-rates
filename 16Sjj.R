###稀疏化表格
library(permute)
library(vegan)
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果")
Bfn1<-read.csv("JJ_OTU16s.csv", header=TRUE, row.names=1)
Bfn<-t(Bfn1)
S1 <- specnumber(Bfn)
(raremax1 <- min(rowSums(Bfn)))
newOTU <- rrarefy(Bfn, raremax1)
newOTU1 = t(newOTU)
write.csv(newOTU1,file='JJ_16S_rare.csv')

############################################################################
############################### NMDS分析  #####################################
###########################################################################
group <-read.csv("16S分组.csv", header=T, row.names = 1)
otu<-read.csv("JJ_16S_rare.csv", header=TRUE, row.names=1)
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
write.csv(sample_site,file='sample_site_16S.csv')

library(ggplot2)
p <- ggplot(sample_site, aes(NMDS1, NMDS2, group = group)) +
  geom_point(aes(color = group ), size = 5, alpha =1) + #可在这里修改点的透明度、大小
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

################################################################################################
###############################################################################################

####alpha diversity
library(permute)
library(vegan)
rare1<-read.csv("JJ_16S_rare.csv", header=TRUE, row.names=1)
rare<-t(rare1)
S <- specnumber(rare)
write.csv(S,file='16SJJ_richness.csv')

shannon <- diversity(rare, index = "shannon")
write.csv(shannon,file='16SJJ_shannon_lt.csv')

invsimpson <- diversity(rare, index = "invsimpson")
write.csv(invsimpson,file='16SJJ_invsimpson_lt.csv')

alpha = data.frame(S,shannon,invsimpson,check.names = T)
write.csv(alpha,file='16SJJ_alpha_lt.csv')

###Chao1 ACE
X=estimateR(rare)
write.csv(X,file='Chao_lt.csv')

###Pielou's evenness (J
J <- shannon/log(S)
write.csv(J,file='J_lt.csv')


#####################################################################################
###################################richness箱线图  #########################################
######################################################################################
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggsci)
species<-read.csv("JJRrichness_box.csv",header=T,row.names = 1)
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
ggsave('16S_JJrichness_箱线图.pdf', p4, width = 7, height = 6)


##################################################################################
############################  分类求和  ##########################################
#################################################################################
library(permute)
library(vegan)
OTU_TaxITSCP <-read.csv("JJ_16S_rare+taxa.csv", header=TRUE)

Domain=aggregate(x=OTU_TaxITSCP[,2:181]/53251,by=list(OTU_TaxITSCP$Kingdom),FUN='sum')
write.csv(Domain, file = " Domain_16SJJ.csv")

Phylum=aggregate(x=OTU_TaxITSCP[,2:181]/53251,by=list(OTU_TaxITSCP$Phylum),FUN='sum')
write.csv(Phylum, file = " Phylum_16SJJ.csv")

Class=aggregate(x=OTU_TaxITSCP[,2:181]/53251,by=list(OTU_TaxITSCP$Class),FUN='sum')
write.csv(Class, file = " Class_16SJJ.csv")

Order=aggregate(x=OTU_TaxITSCP[,2:181]/53251,by=list(OTU_TaxITSCP$Order),FUN='sum')
write.csv(Order, file = " Order_16SJJ.csv")

Family=aggregate(x=OTU_TaxITSCP[,2:181]/53251,by=list(OTU_TaxITSCP$Family),FUN='sum')
write.csv(Family, file = " Family_16SJJ.csv")

Genus=aggregate(x=OTU_TaxITSCP[,2:181]/53251,by=list(OTU_TaxITSCP$Gene),FUN='sum')
write.csv(Genus, file = " Genus_16SJJ.csv")

##############################################
###相似性
library(permute)
library(vegan)
library(reshape)

ado1<-read.csv("JJ_16S_rare.csv", header=TRUE, row.names=1)
ado=(t(ado1)/53251)^0.5

prok.dist <- vegdist(ado)

prok.matr<- melt(as.matrix(prok.dist)) #Converting the dist object to matrix using melt
p    <- t(apply(prok.matr[,c(1,2)],1,FUN=sort))

#p = prok.matr[order(prok.matr[,1],prok.matr[,2]),另一种排序方法

rmv1 <- which(p[,1] == p[,2])

p.new <- paste(p[,1],p[,2],sep="|")

rmv2 <- which(duplicated(p.new))

prok.matr1   <- prok.matr[-c(rmv1,rmv2),] 

write.csv(prok.matr1, file = "相似性1.csv")#根据相似性表格筛选出不同的时间点的相似性


library(ggplot2)
library(ggpubr)
library(ggpmisc)
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))
df_1 <-read.csv("时间衰减.csv", header=TRUE, row.names = 1)
##绘制基本的相关散点图
##geom_point()：用于创建散点图。关键参数：color，size和shape，用于设置散点的边缘颜色，大小和形状。
##geom_smooth()：用于添加回归曲线和曲线拟合的置信区间。关键参数：
#color、size和linetype：指定回归曲线的颜色、粗细和线型。
#fill：指定置信区间的填充颜色。
#method：用于指定曲线拟合的方法。
##一个组的相关分析
# Scatter plot with regression line
b1 = ggplot(df_1, aes(x = group, y =similarity))+ geom_point(color = "black",shape = 16,size = 3)+
  geom_smooth(method = "lm", color = "black", fill = "lightgray",size = 1)+
  stat_cor(method = "pearson", label.x = 30, label.y =0.7,size=5)+##添加显著性水平，这里可以改变相关性的方法
  theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ #添加黑色外框
  xlab("inclubation time (day)")+
  ylab("Bray-Curtis similarity")+
  theme(axis.text.x = element_text(size = 15,color="black"))+ ##设置X轴刻度上的字体大小
  theme(axis.text.y  = element_text(size = 15,color="black"))+##设置y轴刻度上的字体大小
  theme(axis.title.y = element_text(size = 15, color="black"))+#设置Y轴标题大小
  theme(axis.title.x  = element_text(size = 15, color="black"))+#设置x轴标题大小
  theme(legend.text=element_text(size=18, color="black"))
#theme_grey()
b1
ggsave('距离衰减-相似性1.pdf', b1, width = 6, height = 6)


library(RColorBrewer)
colourCount <-  length(unique(df_1$ID))
b2 = ggplot(df_1, aes(x = group, y =similarity,color=factor(ID))) + geom_point(aes(color = factor(ID)),size = 3.5)+
  geom_smooth(aes(color = factor(ID), fill = factor(ID)), method = "lm",size = 1.5) +
  geom_rug(aes(color = factor(ID))) +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(colourCount))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(colourCount))+
  ggpubr::stat_cor(method = "pearson",aes(color = ID), label.x = 30,label.y = 0.85,size = 5)+#method = "spearman",
  xlab("inclubation time (day)")+
  ylab("Bray-Curtis similarity")+
  theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ #添加黑色外框
  theme(axis.text.x = element_text(size = 19,color="black"))+ ##设置X轴刻度上的字体大小
  theme(axis.text.y  = element_text(size = 19,color="black"))+##设置y轴刻度上的字体大小
  theme(axis.title.y = element_text(size = 19, color="black"))+#设置Y轴标题大小
  theme(axis.title.x  = element_text(size = 19, color="black"))+#设置x轴标题大小
  theme(legend.text=element_text(size=18, color="black"))+#图例字体大小
  theme(legend.position = "none")

#ylim(0,1)
b2

ggsave('距离衰减-相似性1.2.pdf', b2, width = 6, height = 6)




#####################################################################################
##############################          ###############################################
#############################  降解速率  ##############################################
#############################           #################################################
#####################################################################################
#######和拷贝数相乘
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/降解拟合")
otu1<-read.csv("JJ_16S_rare.csv", header=TRUE, row.names=1)
copy<-read.csv("copy.csv", header=TRUE, row.names=1)
otu2 <- otu1[,1:180]/53251
le <- c(1:180)
data <- c()
for (i in le) {
  copy1 <- as.numeric(otu2[,i])*as.numeric(copy[i,])
  data <- as.data.frame(cbind(data,copy1))
}

write.csv(data, file = "降解_copy.csv") ###手动修改表头
########################################################
######log转化
table<-read.csv("降解_copy.csv", header=TRUE, row.names=1)
table1<-log(table)
table1[sapply(table1, is.infinite )] <- NA  ##Inf变为NA
table1[is.na(table1)] = 0  ###NA值变为0
head(table1)
write.csv(table1, file = "降解_copy_log.csv")###文件生成后，手动修改标签和乱码


library(basicTrendline)
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/降解拟合")
ado1<-read.csv("降解_copy_log.csv", header=TRUE, row.names=1)
ado2<- ado1[,c(1:6)]
x_arry <- c(1,3,6,12,24,48)
data4 <- c()
j <- seq(1,180,6)
data5 <- c()
for (k in 1:30){
  m = j[k]
  ado2 <- ado1[,c(m,m+1,m+2,m+3,m+4,m+5)]
  for (i in 1:89321){###根据行数进行调整
    y_arry <- as.numeric(ado2[i,])
    data1=data.frame(x=x_arry,y=y_arry) 
    lm.data1<-lm(y~ x,data=data1)
    lm.data2 <- lm.data1$coefficients
    lm.data3 <- as.data.frame(lm.data2)[2,]
    data4 <- rbind(data4,lm.data3)
  }
  data5 <- cbind(data5,data4)
  data4 <- c()
}  
 write.csv(data5, file = "16S降解速率表.csv") ##### 文件生成后，手动修改行名和列名
 data6 <- 0.693/data5*-1
 data6[sapply( data6, is.infinite )] <- 0  ##Inf变为0
 #table1[is.na(table1)] = 0  ###NA值变为0
 write.csv(data6, file = "16S降解半衰期.csv")
 
 
 #######################################################################################################
 ###################################  半衰期分段热图  ###################################################
######################################################################################################

 library(corrplot)
 library(matlab)
 library(ggplot2)  
 library(RColorBrewer)  
 library(reshape2) 
 setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/降解拟合")
 mat<-read.csv("分段表.csv",header = T,row=1)
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
 ggsave('半衰期分段.pdf', p, width = 17, height = 9)
 
 
 ###################################################################################
 #R值P值
 library(reshape2)
 library(ggplot2)
 setwd("E:\\R\\DNA\\heatmap\\")
 data<-read.csv("R值.csv",header = T,row=1)###根据出图调整表格内容
 data$ID <- rownames(data)
 data_m <- melt(data, id.vars=c("ID"))
 
 head(data_m)
 p <- ggplot(data_m, aes(x=variable,y=ID)) + 
   xlab("") + theme_bw() + theme(panel.grid.major = element_blank()) + 
   theme(legend.key=element_blank())  + 
   theme(axis.text.x=element_text(size = 30,face = "bold",angle=45,hjust=1, vjust=1)) + 
   theme(axis.text.y = element_text(size = 30,face = "bold"))+#Y轴坐标
   theme(legend.position="right") +  geom_tile(aes(fill=value)) + 
   theme(legend.text=element_text(size=15, face = "bold"))+#图例字体
   scale_fill_gradient2(low = "blue", high = "red")+
   theme(legend.title=element_text(face="bold",size=15))+
   ylab(" ")
 #guides(fill=FALSE)#移除图例
 
 
 p
 ggsave('Rzhi_1.pdf', p, width = 4, height = 15)

 
 #########################################################################################
 ##################################               ##########################################
 ##################################   mantel分析  #############################################
 ##################################               #############################################
 ##########################################################################################
 library(ggcor)
 library(ggplot2)
 library(RColorBrewer)
 setwd("E:\\毕业数据\\数据处理R\\16S降解\\相关\\")
 spec1 <- t(read.table( "otu降解速率表全.csv", sep = ",", header = T, row.names = 1, check.names = F))
 spec <- spec1*-1
 env <- read.table( "ENV_mantel.csv", sep = ",", header = T, row.names = 1, check.names = F)
 head(env)
 head(spec)
 #组合图
 df <- mantel_test(spec, env,spec.select =list(OTU=1:82244))#分组
 
 df <- df %>% 
   mutate(lty = cut(r, breaks = c(-Inf, 0, Inf), 
                    labels = c("r <= 0", "r > 0")),
          col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                    labels = c("< 0.01", "0.01-0.05", ">= 0.05"),
                    right = FALSE, include.lowest = TRUE))
 #可视化
 options(ggcor.link.inherit.aes = FALSE)
 
 p <- quickcor(env, type = "upper", show.diag = FALSE) + geom_square() +#show.diag = FALSE 保留自相关，FALSE 去除自相关
   anno_link(aes(colour = col, size = lty), data = df,curvature=-0.15) +#curvature改变线的弧度
   #scale_fill_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + ##rev函数表示把颜色反转
   scale_size_manual(values = c(0.01, 0.7)) +#改变线的粗线，因为只设了两个值，所以两个值就够了
   scale_colour_manual(values = c("#ED7D31", "#4472C4", "grey")) +
   guides(size = guide_legend(title = "Mantel's r",
                              override.aes = list(colour = "grey35"),order = 1),#order是图例顺序
          colour = guide_legend(title = "Mantel's p", 
                                override.aes = list(size = 3),order = 2),
          fill = guide_colorbar(title = "Pearson's r", order = 3)) +
   expand_axis(x = -9)+
   #scale_fill_gradient2n() 
   scale_fill_gradientn(colours = rev(brewer.pal(11, "RdBu")))
 p
 ggsave('16S相关mantel.pdf', p, width = 9, height = 7)
 
