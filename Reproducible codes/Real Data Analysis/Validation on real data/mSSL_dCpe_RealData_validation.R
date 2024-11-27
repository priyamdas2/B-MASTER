rm(list=ls())
library(DescTools)
library(devtools)
setwd("U:/BMASTER/Real Data Analysis/Validation on real data")
# install.packages("U:/BMASTER/Simulation artificial/SSLasso/SSLASSO_1.1.tar.gz", repos = NULL, type = "source")
# devtools::install_github(repo = "skdeshpande91/multivariate_SSL/mSSL")
library(SSLASSO)
library(mSSL)
set.seed(1)
## Reading data ################################################################
X <- as.matrix(read.csv('X_Train.csv',header=FALSE))
Y <- as.matrix(read.csv('Y_Train.csv',header=FALSE))

X_Test <- as.matrix(read.csv('X_Test.csv',header=FALSE))
Y_Test <- as.matrix(read.csv('Y_Test.csv',header=FALSE))


N <- dim(X)[1]
p <- dim(X)[2]
q <- dim(Y)[2]


####### mSSL_dCpe analysis ####################################################

fit_mSSL_dCpe <- mSSL_dcpe(X,Y)
B_mSSL_dpe <- fit_mSSL_dCpe$B


B_est <- B_mSSL_dpe
Y_Test_predicted <- X_Test %*% B_est

RMSE <- sqrt(mean((Y_Test - Y_Test_predicted)^2))
MAD <- mean(abs(Y_Test - Y_Test_predicted))
RMSE_MAD <- t(round(c(RMSE,MAD), 3))

write.table(RMSE_MAD,"RMSE_MAD_mSSL_dCpe.csv", sep = ',', 
            row.names = FALSE, col.names = FALSE)


