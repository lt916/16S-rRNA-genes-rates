###############################################################
#######和拷贝数相乘
setwd("E:\\毕业数据\\数据处理R\\降解拟合\\")
otu1<-read.csv("OTU筛选降解.csv", header=TRUE, row.names=1)
copy<-read.csv("copy.csv", header=TRUE, row.names=1)
otu2 <- otu1[,1:180]/52826
le <- c(1:180)
data <- c()
for (i in le) {
  copy1 <- as.numeric(otu2[,i])*as.numeric(copy[i,])
  data <- as.data.frame(cbind(data,copy1))
}

write.csv(data, file = "otu筛选降解_copy.csv")
########################################################
######log转化
table<-read.csv("otu筛选降解_copy.csv", header=TRUE, row.names=1)
table1<-log(table)
write.csv(table1, file = "otu筛选降解_copy_log.csv")###文件生成后，手动修改标签和乱码


library(basicTrendline)
setwd("E:\\毕业数据\\数据处理R\\降解拟合\\")
ado1<-read.csv("otu筛选降解_copy_log.csv", header=TRUE, row.names=1)
ado2<- ado1[,c(1:6)]
x_arry <- c(1,3,6,12,24,48)
data4 <- c()
j <- seq(1,180,6)
data5 <- c()
for (k in 1:30){
  m = j[k]
ado2 <- ado1[,c(m,m+1,m+2,m+3,m+4,m+5)]
for (i in 1:3063){###根据行数进行调整
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
write.csv(data5, file = "otu降解速率表.csv") ##### 文件生成后，手动修改行名和列名

##########################################
######################分类水平
#######和拷贝数相乘
setwd("E:\\毕业数据\\数据处理R\\16S降解\\降解拟合\\分类水平\\")
otu1<-read.csv("Genus_16S.csv", header=TRUE, row.names=1)
copy<-read.csv("copy.csv", header=TRUE, row.names=1)
#otu2 <- otu1[,1:180]/52826
le <- c(1:180)
data <- c()
for (i in le) {
  copy1 <- as.numeric(otu1[,i])*as.numeric(copy[i,])
  data <- as.data.frame(cbind(data,copy1))
}
test <- log(data)#####log转化
write.csv(test, file = "Genus_16S_copy_log.csv")###文件生成后，手动修改标签和乱码


library(basicTrendline)
ado1<-read.csv("Genus_16S_copy_log.csv", header=TRUE, row.names=1)
ado2<- ado1[,c(1:6)]
x_arry <- c(1,3,6,12,24,48)
data4 <- c()
j <- seq(1,180,6)
data5 <- c()
for (k in 1:30){
  m = j[k]
  ado2 <- ado1[,c(m,m+1,m+2,m+3,m+4,m+5)]
  for (i in 1:1630){###根据行数进行调整
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

write.csv(data5, file = "Genus降解速率表.csv") ##### 文件生成后，手动修改行名和列名
#拟合线性函数




#拟合线性函数
#p1 <- trendline(x_arry, y_arry, 
#                model="line2P",
#                linecolor='red',
#                CI.color = NA,
#                col='red',
#                #xlim=c(0,10),ylim=c(0,10),
#                xlab='',ylab='')


##第二种方法

data1=data.frame(x=x_arry,y=y_arry) 
lm.data1<-lm(y~ x,data=data1)
summary(lm.data1)        #输出拟合后信息
lm.data2 <- lm.data1$coefficients
lm.data3 <- as.data.frame(lm.data2)[2,]
