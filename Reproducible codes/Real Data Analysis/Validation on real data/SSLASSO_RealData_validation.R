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
##############################################################################

####### SSLASSO analysis ######################################################

#fit_sslg <- sep_sslg(X,Y)
lambda1 <-  1
nlambda <- 100
lambda0 <- seq(1, nrow(X), length.out = nlambda)
L <- length(lambda0)
B_SSLASSO <- matrix(NA, p, q)

TP_FN_FP_TN <- array(NA, dim = c(nlambda, 4, q))

for (j in 1:q) {
  fit_sslg <- SSLASSO(X, Y[, j], penalty = "adaptive", variance = "fixed", 
                      lambda1, lambda0, nlambda, theta = 0.5, sigma = 1, a = 1, b = ncol(X),  
                      eps = 0.001, max.iter = 500,  counter = 10, warn = FALSE)
  B_SSLASSO[, j] <- fit_sslg$beta[, L]
}

Y_Test_predicted <- X_Test %*% B_SSLASSO

RMSE <- sqrt(mean((Y_Test - Y_Test_predicted)^2))
MAD <- mean(abs(Y_Test - Y_Test_predicted))
RMSE_MAD <- t(round(c(RMSE,MAD), 3))

write.table(RMSE_MAD,"RMSE_MAD_SSLasso.csv", sep = ',', 
            row.names = FALSE, col.names = FALSE)
