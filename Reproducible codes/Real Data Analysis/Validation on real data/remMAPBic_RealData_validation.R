rm(list=ls())
library(DescTools)
setwd("U:/BMASTER/Real Data Analysis/Validation on real data")
# install.packages("U:/BMASTER/Simulation artificial/Remmap/remMap_0.2-0.tar.gz", repos = NULL, type = "source")
library(remMap)

set.seed(1)
## Reading data ################################################################
X <- as.matrix(read.csv('X_Train.csv',header=FALSE))
Y <- as.matrix(read.csv('Y_Train.csv',header=FALSE))

X_Test <- as.matrix(read.csv('X_Test.csv',header=FALSE))
Y_Test <- as.matrix(read.csv('Y_Test.csv',header=FALSE))


N <- dim(X)[1]
p <- dim(X)[2]
q <- dim(Y)[2]


####### remMap analysis ######################################################

lamL1.v <- exp(seq(log(1), log(1000), length = 11))
lamL2.v <- seq(0, 1000, length = 11)

set.seed(1)
remMap_model <- remMap.BIC(X=X, Y=Y,lamL1.v, lamL2.v, C.m = NULL)
pick <- which.min(as.vector(t(remMap_model$BIC)))
result0 <- remMap_model$phi[[pick]]


B_est <- result0$phi

Y_Test_predicted <- X_Test %*% B_est

RMSE <- sqrt(mean((Y_Test - Y_Test_predicted)^2))
MAD <- mean(abs(Y_Test - Y_Test_predicted))
RMSE_MAD <- t(round(c(RMSE,MAD), 3))

write.table(RMSE_MAD,"RMSE_MAD_remMapbic.csv", sep = ',', 
            row.names = FALSE, col.names = FALSE)


