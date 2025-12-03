gc()

{
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(limma)
library(tidyverse)
library(dplyr)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
}

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 加载数据集
TCGA<-read.csv("TCGATimeExp.csv",header = T,row.names = 1)
TCGA <- TCGA[TCGA$OS.time >= 30, ]
TCGA$OS.time=TCGA$OS.time/365

mm <- lapply(mm,function(x){
  x[,-c(1:2)] <- scale(x[,-c(1:2)])
  return(x)})

result <- data.frame()

# TCGA作为训练集
est_data <- mm$TCGA
# GEO作为验证集
val_data_list <- mm
pre_var <- colnames(est_data)[-c(1:2)]  #注意保留基因
est_dd <- est_data[, c('OS.time', 'OS', pre_var)]
val_dd_list <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', pre_var)]})

# 设置种子数和节点数，其中节点数可以调整
rf_nodesize <- 5
seed <- 123


{
##################################
#### 1-1.RSF ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result, cc)


##################################
#### 1-2.RSF + CoxBoost ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]), 
                            trace=TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]), 
                      maxstepno = 500, K = 10, type = "verweij",  penalty = pen$penalty)
fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]), 
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1],  newstatus = x[, 2], type = "lp")))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CoxBoost')
result <- rbind(result, cc)

##################################
#### 1-3.RSF + Enet ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

##################################
#### 1-4.RSF + GBM ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)

# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'GBM')
result <- rbind(result, cc)

##################################
#### 1-5.RSF + Lasso ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Lasso')
result <- rbind(result, cc)

##################################
#### 1-6.RSF + plsRcox ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)

fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'plsRcox')
result <- rbind(result, cc)

##################################
#### 1-7.RSF + Ridge ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 0)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Ridge')
result <- rbind(result, cc)

##################################
#### 1-8.RSF + StpCox ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time, OS)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}



##################################
#### 1-10.RSF + survival-SVM ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
fit = survivalsvm(Surv(OS.time, OS)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'survival-SVM')
result <- rbind(result,cc)

##################################
#### 2.Enet ####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = 'link', newx = as.matrix(x[,-c(1,2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

##################################
#### 3.StepCox ####
##################################
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd), direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction)
  rid <- names(coef(fit))
  #这里不用卡P值，迭代的结果就是可以纳入的基因
  est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
  val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                              trace=TRUE, start.penalty = 500, parallel = T)
  cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                        maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
  fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                  stepno = cv.res$optimal.step, penalty = pen$penalty)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime=x[, 1], newstatus=x[,2], type="lp")))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + CoxBoost')
  result <- rbind(result, cc)
  
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha, nfolds = 10)
    rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
    cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')
    result <- rbind(result, cc)
  }
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  # find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
  cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + GBM')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold=10, #例文描述：10-fold cross-validation
                  family = "cox", alpha = 1)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Lasso')
  result <- rbind(result, cc)
  set.seed(seed)
  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time,
                 event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1,2)])))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + plsRcox')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold = 10, #例文描述：10-fold cross-validation
                  family = "cox", alpha = 0)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Ridge')
  result <- rbind(result, cc)
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
               ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + RSF')
  
  result <- rbind(result, cc)
  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
               censoring.status = est_dd2$OS,
               featurenames = colnames(est_dd2)[-c(1,2)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                       n.fold = 10,
                       n.components = 3,
                       min.features = 2,
                       max.features = nrow(data$x),
                       compute.fullcv = TRUE,
                       compute.preval = TRUE)
  rs <- lapply(val_dd_list2, function(w){
    test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2], RS = rr)
    return(rr2)
  })
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + SuperPC')
  result <- rbind(result, cc)
  fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + survival-SVM')
  result <- rbind(result, cc)
}

##################################
#### 4-3.CoxBoost + GBM####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'GBM')
result <- rbind(result, cc)

##################################
#### 4-4.CoxBoost + lasso####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)

cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno=cv.res$optimal.step, penalty=pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Lasso')
result <- rbind(result, cc)

##################################
#### 4-5.CoxBoost + plsRcox####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)

fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type="lp", newdata = x[, -c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'plsRcox')
result <- rbind(result, cc)

##################################
#### 4-6.CoxBoost + Ridge####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K=10, type="verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 0)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Ridge')
result <- rbind(result, cc)

##################################
#### 4-7.CoxBoost + StepCox####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

##################################
#### 4-8.CoxBoost + SuperPC####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
data <- list(x = t(est_dd2[, -c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 2,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval =TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x=t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'SuperPC')
result <- rbind(result, cc)

##################################
#### 4-9.CoxBoost + survival-SVM####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[,-c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
fit = survivalsvm(Surv(OS.time, OS)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'survival-SVM')
result <- rbind(result, cc)

##################################
#### 5.plsRcox####
##################################
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd[,pre_var], time = est_dd$OS.time, status = est_dd$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd[,pre_var], time = est_dd$OS.time, event = est_dd$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result, cc)

##################################
#### 6.superpc####
##################################
data <- list(x = t(est_dd[, -c(1,2)]), y = est_dd$OS.time, censoring.status = est_dd$OS, featurenames = colnames(est_dd)[-c(1, 2)])
set.seed(seed) 
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 2,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result, cc)

##################################
#### 7.GBM####
##################################
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result, cc)

##################################
#### 8.survival-SVM####
##################################
fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd, gamma.mu = 1)
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('survival - SVM')
result <- rbind(result, cc)

##################################
#### 9.Ridge####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
cvfit = cv.glmnet(x1, x2,
                  nfold = 10, #例文描述：10-fold cross-validation
                  family = "cox")

rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = cvfit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Ridge')
result <- rbind(result, cc)

##################################
####10.Lasso####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = 'cox', alpha = 1)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso')
result <- rbind(result, cc)

##################################
#### 10.1.Lasso + CoxBoost####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[,-c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result, cc)

##################################
#### 10.2.Lasso + GBM####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'GBM')
result <- rbind(result, cc)

##################################
#### 10.3.Lasso + plsRcox####
##################################

x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[,-c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'plsRcox')
result <- rbind(result, cc)

##################################
#### 10.4.Lasso + RSF####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid<-rid[-1]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
             ntree = 1000, nodesize = rf_nodesize, ##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso', ' + RSF')
result <- rbind(result, cc)

##################################
#### 10.5.Lasso + stepcox####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

##################################
#### 10.6.Lasso + superPC####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
data <- list(x = t(est_dd2[,-c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 2,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'SuperPC')
result <- rbind(result, cc)



#######################
#### 10.7.Lasso + survival-SVM####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1)
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'survival-SVM')
result <- rbind(result, cc)

}

save(result, file = "model-result.rda") #保存结果




load("model-result.rda")
#将得到的结果赋给result2变量进行操作
result2 <- result

duplicates <- result2 %>%
  group_by(ID, Model) %>%
  filter(n() > 1)

###将结果的长数据转换为宽数据
dd2 <- pivot_wider(result2, names_from = 'ID', values_from = 'Cindex') %>% as.data.frame()
#将C指数定义为数值型
dd2[,-1] <- apply(dd2[,-1], 2, as.numeric)
#求每个模型的C指数在2个数据集的均值
dd2$All <- apply(dd2[,2:4], 1, mean)
###查看每个模型的C指数
head(dd2)

#输出C指数结果
write.table(dd2,"output_C_index.txt", col.names = T, row.names = F, sep = "\t", quote = F)

dd2=read.table('output_C_index.txt', header=T, sep="\t", check.names=F)
dt <- dd2[, 2:4]

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

SimpleHeatmap <- function(Cindex_mat, avg_Cindex, 
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5, 
                          cluster_columns, cluster_rows){
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)
  
  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = T, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)
  
  Heatmap(as.matrix(Cindex_mat), name = "AUC",
          right_annotation = row_ha, 
          top_annotation = col_ha,
          # col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 黄绿配色
          col = c("#4195C1", "#FFFFFF", "#CB5746"), # 红蓝配色
          rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
          cluster_columns = cluster_columns, cluster_rows = cluster_rows, # 不进行聚类，无意义
          show_column_names = FALSE, 
          show_row_names = TRUE,
          row_names_side = "left",
          width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 10))
          }
  )
}


# 根据C指数排序
dd2 <- dd2[order(dd2$All, decreasing = T),]
#仅绘制GEO验证集的C指数热图
dt <- dd2[, 2:4]
rownames(dt) <- dd2$Model

##热图绘制
Cindex_mat=dt
avg_Cindex <- apply(Cindex_mat, 1, mean)     # 计算每种算法在所有队列中平均C-index
avg_Cindex <- sort(avg_Cindex, decreasing = T)     # 对各算法C-index由高到低排序
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]      # 对C-index矩阵排序

avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
if(ncol(Cindex_mat) < 3) { # 如果用于绘图的队列小于3个
  CohortCol <- c("red","blue") # 则给出两个颜色即可（可自行替换颜色）
} else { # 否则通过brewer.pal赋予超过3个队列的颜色
  CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") # 设置队列颜色
}
names(CohortCol) <- colnames(Cindex_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat, # 主矩阵
                    avg_Cindex, # 侧边柱状图
                    CohortCol, "steelblue", # 列标签颜色，右侧柱状图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path( "Cindex_epma.pdf"), width = cellwidth * ncol(Cindex_mat) + 5, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm)
invisible(dev.off())





direction='forward'
fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = direction)
rid <- names(coef(fit))

est_dd2 <- est_data[,c('OS.time','OS',rid)]
val_dd_list2 <- lapply(val_data_list,function(x){x[,c('OS.time','OS',rid)]})

set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(fit$err.rate)
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd2,
             ntree = best,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})


