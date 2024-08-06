#############基于丰度相关性的微生物共发生网络
##计算微生物丰度间的相关系数（以spearman相关系数进行计算）
setwd("E:/毕业数据/数据处理R/16S降解/网络")
library(Hmisc)

otu <- read.csv('第48天.csv', header=T,row.names=1)
#otu <- read.delim('ado_s1.txt', row.name = 1, check.names = FALSE)##这里输入的表已经是相对丰度总和高于 0.005的OTU了 的
otu <- otu/52826
#可选事先过滤一些低丰度或低频的类群
otu <- otu[which(rowSums(otu)/30 >= 0.0001), ]  

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
#head(edge_list)

#write.table(edge_list, 'network.edge_HM_0.250.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(edge_list,file='network.edge_6_0.8.csv', row.names = F)



#节点属性列表，对应边列表，记录节点属性，例如
node_list <- data.frame(
  nodes_id = V(g)$name,    #节点名称
  normalized_degree = degree(g,v = V(g),normalized = T),    #节点度
  betweenness = betweenness(g,v = V(g)),normalized = T)
write.csv(node_list,file='network.node_6_0.8.csv',row.names = T)
 #---------------------------------------------------------------------------------
# network property
# 边数量 The size of the graph (number of edges)
num.edges = length(E(g)) # length(curve_multiple(igraph))
num.edges
# 顶点数量 Order (number of vertices) of a graph
num.vertices = length(V(g))# length(diversity(igraph, weights = NULL, vids = V(igraph)))
num.vertices
# 连接数(connectance) 网络中物种之间实际发生的相互作用数之和（连接数之和）占总的潜在相互作用数（连接数）的比例，可以反映网络的复杂程度
connectance = edge_density(g,loops=FALSE)# 同 graph.density;loops如果为TRUE,允许自身环（self loops即A--A或B--B）的存在
#connectance
# 平均度(Average degree)
average.degree = mean(igraph::degree(g))# 或者为2M/N,其中M 和N 分别表示网络的边数和节点数。
#average.degree
# 平均路径长度(Average path length)
average.path.length = average.path.length(g) # 同mean_distance(igraph) # mean_distance calculates the average path length in a graph
#average.path.length
# 直径(Diameter)
diameter = diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)
#diameter
# 边连通度 edge connectivity / group adhesion
edge.connectivity = edge_connectivity(g)
#edge.connectivity

#节点连通度（nodes connectivity）
nodes_connectivity <- vertex.connectivity(g)
#nodes_connectivity

# 聚集系数(Clustering coefficient)：分局域聚类系数和全局聚集系数，是反映网络中节点的紧密关系的参数，也称为传递性。整个网络的全局聚集系数C表征了整个网络的平均的“成簇性质”。
clustering.coefficient = transitivity(g) 
clustering.coefficient
no.clusters = no.clusters(g)
no.clusters
#图密度（density）
graph_density <- graph.density(g)
graph_density

# 介数中心性(Betweenness centralization)
centralization.betweenness = centralization.betweenness(g)$centralization 
centralization.betweenness
# 度中心性(Degree centralization)
centralization.degree = centralization.degree(g)$centralization
centralization.degree

network_property <-c(num.edges=num.edges,num.vertices=num.vertices,connectance=connectance,
                     average.degree=average.degree,average.path.length=average.path.length,
                     diameter=diameter,edge.connectivity=edge.connectivity,clustering.coefficient=clustering.coefficient,
                     no.clusters=no.clusters,centralization.betweenness=centralization.betweenness,
                     centralization.degree=centralization.degree)
network_property <- data.frame(network_property)
write.csv(network_property, 'network_property_CD_0.025.csv', row.names = T)

density(E(g), loops=TRUE)
