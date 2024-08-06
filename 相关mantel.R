
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

