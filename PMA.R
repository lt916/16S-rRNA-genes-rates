###稀疏化表格
library(permute)
library(vegan)
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/PMA")
Bfn1<-read.csv("PMA_OTU.csv", header=TRUE, row.names=1)
Bfn<-t(Bfn1)
S1 <- specnumber(Bfn)
(raremax1 <- min(rowSums(Bfn)))
newOTU <- rrarefy(Bfn, raremax1)
newOTU1 = t(newOTU)
write.csv(newOTU1,file='PMA_rare.csv')

#################################### NMDS图 ###################################
##################################################################################
library(vegan)	
library(ggplot2)	
library(ape)
#读入 OTU 丰度表
otu <- read.csv("PMA_rare.csv", header=TRUE, row.names=1)
#otu <- t(otu)
otu=(t(otu)/26729)^0.5
#排序（显示 4 个排序轴）
nmds1 <- metaMDS(otu, autotransform =FALSE)
#提取应力函数值（stress）
nmds1.stress <- nmds1$stress
#提取样本排序坐标
nmds1.point <- data.frame(nmds1$point)
#提取物种（OTU）排序坐标
nmds1.species <- data.frame(nmds1$species)
##ggplot2作图
###读入分组文件
group <- read.csv('group_PMA.csv', header=TRUE)
#提取样本点坐标（前两轴）
sample_site <- nmds1.point[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
#为样本点坐标添加分组信息
sample_site <- merge(sample_site, group, by = 'names', all.x = TRUE)
#write.csv(sample_site, file = "sample_site.csv")
#使用 ggplot2 绘制 NMDS 排序图
library(ggthemes)
library(RColorBrewer)
myPalette=colorRampPalette(brewer.pal(12,"Paired"))(29)##设置颜色
p <- ggplot(sample_site, aes(NMDS1, NMDS2, group = group)) +
  geom_point(aes(color = factor(Site),shape = group), size = 5, alpha =1) + #可在这里修改点的透明度、大小
  scale_shape_manual(values = c(17,16)) + #可在这里修改点的形状
  #scale_color_manual(values = c("#80b1d3","#8dd3c7","#B4ADE3","#FFA45B","#F76248"))+ #可在这里修改点的颜色
  scale_color_manual(values = myPalette)+
  #scale_fill_brewer(palette = "Set3")+
  #scale_color_tableau()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) + #去掉背景框
  theme(legend.key = element_rect(fill = 'transparent'), legend.title = element_blank()) + #去掉图例标题及标签背景
  labs(x = 'NMDS 1', y = 'NMDS 2', title = paste('Stress =', round(nmds1$stress,3))) +
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
ggsave('NMDS_PMA.pdf',p, width =8, height = 6)


###################################################################################
####################################  火山图  ########################################
######################################################################################
otu_file <- read.csv("PMA_rare.csv", row.names=1)
group_file <- read.csv('group_PMA.csv', row.names = 1)#分组文件

#otu_file1 <- read.csv('otutab_raw.csv', row.names = 1)
#group_file1 <- read.csv('group.csv', row.names = 1)#分组文件

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = otu_file, colData = group_file, design = ~group)#构建 DESeqDataSet 对象  
dds <- DESeq(dds) #差异分析
suppressMessages(dds)

#结果查看
res <- results(dds, contrast=c('group', 'Total', 'Active'))#提取分析结果
res = res[order(res$pvalue),]
res #查看结果
summary(res)  #简要统计结果
table(res$padj<0.05) #查看fdr校正后的P<0.05的个数
#合并结果
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.table(resdata,file= "DESeq2_PMA.txt",sep="\t",quote=F,row.names=F)
write.csv(resdata,file='DESeq2_PMA.csv')
#火山图绘制
resdata[which(resdata$pvalue < 0.05 & resdata$log2FoldChange <= -1),'sig'] <- 'diff'
resdata[which(resdata$pvalue < 0.05 & resdata$log2FoldChange >= 1),'sig'] <- 'diff'
resdata[which(resdata$pvalue >= 0.05 | abs(resdata$log2FoldChange) < 1),'sig'] <- 'no diff'
#for (i in 1:nrow(resdata)) {
#  if (abs(resdata[i,'log2FoldChange']) >= 0.5) resdata[i,'select_change'] <- 'y' else resdata[i,'select_change'] <- 'n'
#  if (resdata[i,'padj'] %in% NA | abs(resdata[i,'padj']) >= 0.05) resdata[i,'select_pvalue'] <- 'n' else resdata[i,'select_pvalue'] <- 'y'
#  resdata[i,'select'] <- paste(resdata[i,'select_change'], resdata[i,'select_pvalue'], sep = '')
#}
library(ggplot2)
##ggplot2 差异火山图
#resdata$select <- factor(resdata$select, levels = c('depleted', 'enriched', 'no diff'), labels = c('p >= 0.05, FC < 2', 'p < 0.05, FC < 2', 'p >= 0.05, FC >= 2', 'p < 0.05, FC >= 2'))

#纵轴为显著性 p 值
p <- ggplot(resdata, aes(x = log2FoldChange, y=-log10 (pvalue),color = sig)) +
  geom_point(aes(color = sig), alpha = 1) +
  
  scale_color_manual(values = c('#e64e2f',  'gray60')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
  labs(x = 'log2FoldChange', y = '-log10 (pvalue)')
p
ggsave('火山图.pdf', p, width = 6, height = 5)


#########################################################################################
################################## 生态位宽度 ############################################
##########################################################################################
install.packages("devtools")
library(devtools)
devtools::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)

###################直接安装就行
install.packages(spaa)
###############################先去掉全是0的OTU——（因为我是将总的OTU表拆成很多个分别分析的，所以要进行这一步）
rm(list = ls())
data1 <- read.csv("PMA_rare.csv", header=TRUE, row.names = 1)

############################一个文件里的不同处理，例如这个是PMA处理
b1 <- data1[,1:29]
k1 <- b1[which(rowSums(b1) > 0),]
#write.csv(k1, file = "PMA_remv0.csv", row.names = T
#################################################### 生态位宽度计算

#rm(list = ls())
library(spaa)
#ado <- read.csv("PMA_remv0.csv", header=TRUE, row.names=1)
ado1 <- t(k1)
niche_width1 <- niche.width(ado1, method = 'levins')
niche_width_PMA <- t(niche_width1)
write.csv(niche_width_PMA, 'niche_width_16sPMA.csv')
############################一个文件里的不同处理，例如这个是No_PMA
b2 <- data1[,30:58]
k2 <- b2[which(rowSums(b2) > 0),]
ado2 <- t(k2)
niche_width2 <- niche.width(ado2, method = 'levins')
niche_width_NPMA <- t(niche_width2)
write.csv(niche_width_NPMA, 'niche_width_16sNPMA.csv')

###########################################################
##########出图
library(ggplot2)
library(Rmisc)
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/PMA")
data<-read.csv("niche_width_ALL.csv",header = T,row=1)
input <- summarySE(data, measurevar="WIDTH", groupvars=c("GROUP"))

g <- ggplot(input, aes(GROUP, WIDTH, fill=GROUP)) + 
  geom_bar(position=position_dodge(), stat="identity",width=.8) +
  geom_errorbar(aes(ymin=WIDTH-se, ymax=WIDTH+se), width=.3, position=position_dodge(.9)) +  
  scale_fill_manual(values = c("#2958a2","#c6472e"))+
  theme_bw() +theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(color="black",size = 25,face = "plain", angle = 0,vjust = 0.5))+#X轴坐标，size 大小， face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
  theme(axis.text.y = element_text(color="black",size = 25,face = "plain"))+#Y轴坐标
  xlab("Group") + theme(axis.title.x = element_text(color="black",size = 25, face = "plain"))+#X轴标题
  ylab("Niche_width") + theme(axis.title.y = element_text(size = 27, face = "plain"))+#Y轴标题
  theme(panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+#边框加粗
  theme(axis.ticks.x=element_line(color="black",size=1.1,lineend = 1))+#x轴坐标刻度
  theme(axis.ticks.y=element_line(color="black",size=1.1,lineend = 1))+#y轴坐标刻度
  guides(fill=FALSE)+ylim(0.0,2.5)#移除图例
#theme(legend.text=element_text(size=13, face = "bold"))+#图例字体
#guides(fill=guide_legend(title=NULL))#移除图例标题
g 
wilcox.test( WIDTH ~ GROUP,data = data)##比较差异#####
ggsave('生态位宽度.pdf', g, width = 6.5, height = 7)

#####################比较差异########################################’




####################################################################
############################### 多样性 ################################
#######################################################################
####alpha diversity
library(permute)
library(vegan)
rare1<-read.csv("PMA_rare.csv", header=TRUE, row.names=1)
rare<-t(rare1)
S <- specnumber(rare)
write.csv(S,file='PMA_richness.csv')
shannon <- diversity(rare, index = "shannon")
write.csv(shannon,file='16S_shannon_lt.csv')
invsimpson <- diversity(rare, index = "invsimpson")
write.csv(invsimpson,file='降解16S_invsimpson_lt.csv')
alpha = data.frame(S,shannon,invsimpson,check.names = T)
write.csv(alpha,file='降解16S_alpha_lt.csv')
###Chao1 ACE
X=estimateR(rare)
write.csv(X,file='Chao_lt.csv')
###Pielou's evenness (J
J <- shannon/log(S)
write.csv(J,file='J_lt.csv')

#######################################################
#箱线图
####α多样性\ 
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggsci)
species <-read.csv("PMA_richness.csv", header=TRUE, row.names = 1)
scaleFUN<-  function(x) sprintf ("%.1f",x)
my_comparisons<-list(c("TDNA","IDNA"))
#dev.new()
p4 = ggplot(species,aes(x=factor(type,level=c("TDNA","IDNA")), y=richness/1000,color = type))+
  stat_boxplot(geom="errorbar",width=0.15)+
  geom_boxplot(alpha=0.7,aes(fill=type))+
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
  ylab("Prokaryotic richness(103)") +
  theme(axis.text.x = element_text(size = 20,color="black"))+ ##设置X轴刻度上的字体大小
  theme(axis.text.y  = element_text(size = 20,color="black"))+##设置y轴刻度上的字体大小
  theme(axis.title.x  = element_text(size = 20, color="black"))+#设置x轴标题大小
  #theme(axis.title.y  = element_text(size = 18, color="black"))+#设置y轴标题大小
  ylim(1.0,8.5)+ #设置y轴范围
  stat_compare_means(comparisons=my_comparisons,label ="p.signif",size=8, step_increase = 0.1,map_signif_level = T, method = "t.test", paired = TRUE)+ # Add pairwise 
  stat_compare_means(label.y = 10.6,label.x = 0.6,size=6,method = "t.test",paired = TRUE) # Add global p-value,在这里添加method选择需要进行参数检验的方法，并显示显著性标记到图形上。
#label.y = 50,和label.x = 1.5 是控制图例在图上的位置
#method默认的参数检验方法为"wilcox.test"
p4
#p4<-p4+scale_fill_lancet()
ggsave('PMA_richness.pdf', p4, width = 6, height = 6)



###############################################################################################
#####################################################################################################3
#############基于丰度相关性的微生物共发生网络 ######################################################
##################################################################################################
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/PMA/相互作用")
##计算微生物丰度间的相关系数（以spearman相关系数进行计算）
library(Hmisc)
otu <- read.csv('IDNA.csv', header=T,row.names=1)
#otu <- read.delim('ado_s1.txt', row.name = 1, check.names = FALSE)##这里输入的表已经是相对丰度总和高于 0.005的OTU了 的
otu <- otu/26729
#可选事先过滤一些低丰度或低频的类群
otu <- otu[which(rowSums(otu)/29 >= 0.0002), ]  
otu_corr <- rcorr(t(otu), type = 'spearman')
#相关系数 r 值和显著性 p 值矩阵
r <- otu_corr$r
p <- otu_corr$P
#write.table(r, '16s_forest_corr_r.txt', sep = '\t', row.names = TRUE)
#阈值筛选
#将 spearman 相关系数低于 0.8 的关系剔除，即 r>=0.6
r[abs(r) < 0.8] <- 0
#选取显著性 p 值小于 0.01的相关系数，即 p<0.05
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0
#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
#head(z)[1:6,1:6]
#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
#write.table(data.frame(z, check.names = FALSE), 'genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
##获得网络
library(igraph)
#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')

#自相关也可以通过该式去除
g <- simplify(g)

#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
#边列表
edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表
edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation
)
#write.table(edge_list, 'network.edge_HM_0.250.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(edge_list,file='network.edge_IDNA1.csv', row.names = F)

#节点属性列表，对应边列表，记录节点属性，例如
node_list <- data.frame(
  nodes_id = V(g)$name,    #节点名称
  normalized_degree = degree(g,v = V(g),normalized = T),    #节点度
  betweenness = betweenness(g,v = V(g)),normalized = T)
write.csv(node_list,file='network.node_IDNA1.csv',row.names = T)

#########################################################################################
#################################### 鲁棒性 #############################################
#########################################################################################
library(igraph)
library(ggplot2)
#install.packages("rsq")
library(rsq)
record <- c()
record <- cbind(record,"N_PMA")
result<-c()
##第一条线---------
a=read.csv("network.edge_TDNA1.csv",header = TRUE) ###读取edge文件
a$source=as.factor(a$source)
a$target=as.factor(a$target)
b=a[,c(1:2)]
b=t(b)
b=t(b)
ig=graph_from_edgelist(b,directed=FALSE)
natcon <- function(ig) {
  N   <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc  <- log(mean(exp(evals)))
  nc / (N - log(N))
}
nc.attack <- function(ig) {
  hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=TRUE)
  sapply(1:round(vcount(ig)*.8), function(i) {
    ind <- hubord[1:i]
    tmp <- delete_vertices(ig, V(ig)$name[ind])
    natcon(tmp)
  }) }
nc<- nc.attack(ig)

write.csv(nc,"鲁棒性TDNA1.csv")
p <- read.csv("鲁棒性TDNA1.csv")
colnames(p)[2] <-"NC"
p$RM <- p$X/vcount(ig)
ncmod <- glm(NC~RM,data = p)
summary(ncmod)
record <- cbind(record,ncmod$coefficients[2])
record <- cbind(record,rsq(ncmod))
p1 <- ggplot() + 
  geom_abline(slope = ncmod$coefficients[2], intercept = ncmod$coefficients[1], 
              color = "#c6472e", linewidth = 2, alpha = 1)+
  geom_point(data = p, mapping = aes(x = RM, y = NC), color = "#c6472e",alpha = 1) +
  theme_classic()
p1
result<-rbind(result,record)
##第二条线---------
record <- c()
record <- cbind(record,"PMA")

a=read.csv("network.edge_IDNA1.csv",header = TRUE) ###读取edge文件
a$source=as.factor(a$source)
a$target=as.factor(a$target)
b=a[,c(1:2)]
b=t(b)
b=t(b)
ig=graph_from_edgelist(b,directed=FALSE)
natcon <- function(ig) {
  N   <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc  <- log(mean(exp(evals)))
  nc / (N - log(N))
}
nc.attack <- function(ig) {
  hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=TRUE)
  sapply(1:round(vcount(ig)*.8), function(i) {
    ind <- hubord[1:i]
    tmp <- delete_vertices(ig, V(ig)$name[ind])
    natcon(tmp)
  }) }
nc<- nc.attack(ig)
write.csv(nc,"鲁棒性_IDNA1.csv")
p <- read.csv("鲁棒性_IDNA1.csv")
colnames(p)[2] <-"NC"
p$RM <- p$X/vcount(ig)
ncmod <- glm(NC~RM,data = p)
summary(ncmod)
record <- cbind(record,ncmod$coefficients[2])
record <- cbind(record,rsq(ncmod))


p3 <- p1 + 
  geom_abline(slope = ncmod$coefficients[2], intercept = ncmod$coefficients[1], 
              color = "#2958a2", linewidth = 2, alpha = 1)+
  geom_point(data = p, mapping = aes(x = RM, y = NC), color = "#2958a2",alpha = 1) +
  ylim (0, 0.035) +
  theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ #添加黑色外框
  theme(axis.text.y  = element_text(size = 20,color="black"))+##设置y轴刻度上的字体大小
  theme(axis.title.y = element_text(size = 20, color="black"))+#设置Y轴标题大小
  theme(axis.title.x = element_text(size = 20, color="black"))+#设置x轴标题大小
  theme(legend.text=element_text(size=20, color="black"))+#图例字体大小
  theme(axis.text.x = element_text(size = 20,color="black"))+ ##设置X轴刻度上的字体大小
  theme(legend.text=element_text(size=13, face = "bold"))+
  ylab("Natural connectivity") + xlab("Proportion of removed nodes")
  p3
result<-rbind(result,record)
ggsave("鲁棒性1.pdf",p3,  width = 6, height = 5)
write.csv(result,file='16S_slope&r2.csv')
library(microeco)

#####################################################################################
##########################################################################################
###########################################################################################3
library(ggcor)
library(ggplot2)
library(RColorBrewer)
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/PMA")
spec <- t(read.table( "PMA_mantel.csv", sep = ",", header = T, row.names = 1, check.names = F))

env <- read.table( "理化_mantel.csv", sep = ",", header = T, row.names = 1, check.names = F)
#head(env)
#head(spec)
#组合图
df <- mantel_test(spec, env,spec.select =list(P=1:89321,
                                              T= 89322:178642
))#分组

df <- df %>% 
  mutate(lty = cut(r, breaks = c(-Inf, 0, Inf), 
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "0.01-0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))

#可视化
options(ggcor.link.inherit.aes = FALSE)

p <- quickcor(env, type = "upper",show.diag = FALSE) + geom_square() +
  anno_link(aes(colour = col, size = lty), data = df,curvature=-0.15) +#curvature改变线的弧度
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
ggsave('16S_mantel_PMA1.pdf', p, width = 9, height = 7)



#########################################################################################
#############################                #################################################
#############################  降解—胞外相关 #################################################
#############################                ####################################################
############################################################################################

####################计算胞外的量
###总DNA × 拷贝数
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/速率-胞外")
OTU1<-read.csv("TDNA.csv", header=TRUE, row.names=1)
copy<-read.csv("TDNA_copy.csv", header=TRUE, row.names=1)
OTU_TDNA <- OTU1[,1:29]/26729
le <- c(1:29)
data1 <- c()
for (i in le) {
  copy1 <- as.numeric(OTU_TDNA[,i])*as.numeric(copy[i,])
  data1 <- as.data.frame(cbind(data1,copy1))
}
write.csv(data1, file = "TDNA X COPY.csv")
############胞内DNA × 拷贝数
OTU2<-read.csv("IDNA.csv", header=TRUE, row.names=1)
copy<-read.csv("IDNA_copy.csv", header=TRUE, row.names=1)
OTU_IDNA <- OTU2[,1:29]/26729
le <- c(1:29)
data2 <- c()
for (i in le) {
  copy1 <- as.numeric(OTU_IDNA[,i])*as.numeric(copy[i,])
  data2 <- as.data.frame(cbind(data2,copy1))
}
write.csv(data2, file = "IDNA X COPY1.csv")
###二者相减=胞外DNA 
data1<-read.csv("TDNA X COPY.csv", header=TRUE, row.names=1)
data2<-read.csv("IDNA X COPY.csv", header=TRUE, row.names=1)
EDNA <- data1-data2
EDNA1<- 1- data2/data1
EDNA2<- data1/data2
EDNA3 <- OTU1/OTU2
write.csv(EDNA,file='EDNA.csv')
write.csv(EDNA1,file='EDNA相对.csv')
write.csv(EDNA2,file='TDNA比DNA.csv')
write.csv(EDNA3,file='TDNA比DNA相对.csv')

#########正态分布
shapiro <- shapiro.test(df_1$RATE*1000)
shapiro
library(vegan)
library(ggpubr)
ggqqplot(df_1$RATE)
ggqqplot(df_1$EDNA)
shapiro.test(df_1$RATE1)

############################################################################
#筛选丰富种，单个地点分开筛选
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/速率-胞外")
otu1<-read.csv("JJ_16S_rare.csv", header=TRUE, row.names=1)
otu2<-(otu1[,103:108])/53251

##丰富种
otu_AAT <- otu2[apply(otu2, 1, function(x) min(x)>0.0002), ]
write.csv(otu_AAT, '丰富种S18.csv')



 #######相关分析
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/速率-胞外/作图")
library(ggplot2)
library(ggpubr)
library(ggpmisc)
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))
df_1 <-read.csv("JJ胞外29.csv", header=TRUE, row.names = 1)
df_2=df_1^0.5
##绘制基本的相关散点图
##geom_point()：用于创建散点图。关键参数：color，size和shape，用于设置散点的边缘颜色，大小和形状。
##geom_smooth()：用于添加回归曲线和曲线拟合的置信区间。关键参数：
#color、size和linetype：指定回归曲线的颜色、粗细和线型。
#fill：指定置信区间的填充颜色。
#method：用于指定曲线拟合的方法。
b <- ggplot(df_2, aes(x = T比I, y = RATE))
##一个组的相关分析
# Scatter plot with regression line
b1 = b + geom_point(color = "#236ee6",shape = 16,size = 3)+
  geom_smooth(method = "lm", color = "#236ee6", fill = "lightgray",size = 1)+
  stat_cor(method = "pearson", label.x = 3, label.y =0.28,size=5)+##添加显著性水平，这里可以改变相关性的方法
  theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ #添加黑色外框
  theme(axis.text.x = element_text(size = 15,color="black"))+ ##设置X轴刻度上的字体大小
  theme(axis.text.y  = element_text(size = 15,color="black"))+##设置y轴刻度上的字体大小
  theme(axis.title.y = element_text(size = 15, color="black"))+#设置Y轴标题大小
  theme(axis.title.x  = element_text(size = 15, color="black"))+#设置x轴标题大小
  theme(legend.text=element_text(size=18, color="black"))+xlim(0.01,10)
b1
ggsave('地点28T比I.pdf', b1, width = 6, height = 6)

# Scatter plot with regression line
p <- ggplot(df_1, aes(x = EDNA1, y = RATE1))
b2 = p + geom_point(color = "#236ee6",shape = 16,size = 3)+
  geom_smooth(method = "lm", color = "#236ee6", fill = "lightgray",size = 1)+
  stat_cor(method = "pearson", label.x = 8, label.y =0.04,size=5)+##添加显著性水平，这里可以改变相关性的方法
  theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ #添加黑色外框
  theme(axis.text.x = element_text(size = 15,color="black"))+ ##设置X轴刻度上的字体大小
  theme(axis.text.y  = element_text(size = 15,color="black"))+##设置y轴刻度上的字体大小
  theme(axis.title.y = element_text(size = 15, color="black"))+#设置Y轴标题大小
  theme(axis.title.x  = element_text(size = 15, color="black"))+#设置x轴标题大小
  theme(legend.text=element_text(size=18, color="black"))+xlim(6.2,10.5)
b2
#
ggsave('地点26EDNA.pdf', b2, width = 6, height = 6)


#########################################################################################
#############################                #################################################
#############################    圈圈图      #################################################
#############################                ####################################################
############################################################################################
#分类求和
setwd("D:/研究生/李婷-DNA降解/16S降解/180+TP结果/180+TP结果/PMA/圈圈图")
library(permute)
library(vegan)
OTU_TaxITSCP <-read.csv("PMA_rare+taxa.csv", header=TRUE)

Domain=aggregate(x=OTU_TaxITSCP[,2:59]/26729,by=list(OTU_TaxITSCP$Kingdom),FUN='sum')
write.csv(Domain, file = "Domain_PMA.csv")

Phylum=aggregate(x=OTU_TaxITSCP[,2:59]/26729,by=list(OTU_TaxITSCP$Phylum),FUN='sum')

write.csv(Phylum, file = "Phylum_PMA.csv")

Class=aggregate(x=OTU_TaxITSCP[,2:59]/26729,by=list(OTU_TaxITSCP$Class),FUN='sum')
write.csv(Class, file = "Class_PMA.csv")

Order=aggregate(x=OTU_TaxITSCP[,2:59]/26729,by=list(OTU_TaxITSCP$Order),FUN='sum')
write.csv(Order, file = "Order_PMA.csv")

Family=aggregate(x=OTU_TaxITSCP[,2:59]/26729,by=list(OTU_TaxITSCP$Family),FUN='sum')
write.csv(Family, file = "Family_PMA.csv")

Genus=aggregate(x=OTU_TaxITSCP[,2:59]/26729,by=list(OTU_TaxITSCP$Gene),FUN='sum')
write.csv(Genus, file = "Genus_PMA.csv")


###批量计算配对t-test的 p/t-value值
library(reshape2)
data <-read.csv("Genus_PMA.csv", header=TRUE, row.names=1)

pvalue=apply(data,1,function(x) wilcox.test(x[1:29],x[30:58],paired = T,alternative="two.sided",conf.level=0.95)$p.value)
tvalue=apply(data,1,function(x) t.test(x[1:29],x[30:58],paired = T,alternative="two.sided",conf.level=0.95)$statistic)
write.csv(pvalue, 'G配对P.csv')
write.csv(tvalue, 'G配对T.csv')
