setwd("E:/毕业数据/数据处理R/16S降解/结构方程模型")


library(nlme)
library(lme4)
library(piecewiseSEM)
library(QuantPsyc)


Rt <-read.csv("ENV_MST.csv", header=T,row.names=1)

model1 <- lm(rate ~  MAT + MAP + Altitude, Rt)
coefs(model1, standardize = 'scale')


beta_MAT <-  summary(model1)$coefficients[2, 1]
beta_MAP <- summary(model1)$coefficients[3, 1]
beta_Altitude <-  summary(model1)$coefficients[4, 1]

climate <- beta_MAT  * Rt$MAT  + beta_MAP * Rt$MAP + beta_Altitude * Rt$Altitude
Rt$climate <- climate
summary(lm(rate ~ climate, Rt))
coefs(lm(rate~ climate, Rt))

model2 <- lm(rate ~  pH + NO3 + NH4 + Moisture + TOC + TN + TP +TK, Rt)
coefs(model2, standardize = 'scale')

beta_pH <-  summary(model2)$coefficients[2, 1]
beta_NO3 <- summary(model2)$coefficients[3, 1]
beta_NH4 <-  summary(model2)$coefficients[4, 1]
beta_Moisture <-  summary(model2)$coefficients[5, 1]
beta_TOC <-  summary(model2)$coefficients[6, 1]
beta_TN <-  summary(model2)$coefficients[7, 1]
beta_TP <-  summary(model2)$coefficients[8, 1]
beta_TK <-  summary(model2)$coefficients[9, 1]

soil<- beta_pH  * Rt$pH  + beta_NO3 * Rt$NO3 + beta_NH4 * Rt$NH4+ beta_Moisture * Rt$Moisture+ 
  beta_TOC * Rt$TOC+ beta_TN * Rt$TN+ beta_TP * Rt$TP+ beta_TK * Rt$TK
Rt$soil <- soil
summary(lm(rate ~ soil, Rt))
coefs(lm(rate~ soil, Rt))



model3 <- lm(rate ~  Richness + abundance + MDS1 + MDS2, Rt)
coefs(model3, standardize = 'scale')

beta_Richness <-  summary(model3)$coefficients[2, 1]
beta_abundance<- summary(model3)$coefficients[3, 1]
beta_MDS1 <-  summary(model3)$coefficients[4, 1]
beta_MDS2 <-  summary(model3)$coefficients[5, 1]

Micro<- beta_Richness  * Rt$Richness  + beta_abundance * Rt$abundance + beta_MDS1 * Rt$MDS1+ beta_MDS2 * Rt$MDS2
Rt$Micro <- Micro
summary(lm(rate ~ Micro, Rt))
coefs(lm(rate~ Micro, Rt))


model4 <- lm(rate ~   Grassland + Desert+ Forest  + Cropland , Rt)
summary(model4)$coefficients

beta_Grassland  <-  summary(model4)$coefficients[2, 1]
beta_Desert<- summary(model4)$coefficients[3, 1]
beta_Forest <-  summary(model4)$coefficients[4, 1]
ecosystem <- beta_Grassland * Rt$Grassland + beta_Desert * Rt$Desert +  beta_Forest * Rt$Forest
Rt$ecosystem <- ecosystem
summary(lm(rate ~ ecosystem, Rt))
coefs(lm(rate ~ ecosystem, Rt))


#基于混合效应模型的多元回归
microbe.list <- list(
  lme(Micro ~   climate  + soil + ecosystem , random = ~ 1 | site , na.action = na.omit,
      data = Rt),
  lme(soil ~ climate + ecosystem, random = ~ 1 | site , na.action = na.omit,
      data = Rt),
  lme(rate ~   climate  + soil + Micro + ecosystem, random = ~ 1 | site , na.action = na.omit,
      data = Rt)
)

microbe.psem <- as.psem(microbe.list)
(new.summary <- summary(microbe.psem, .progressBar = F))
plot(microbe.psem)





