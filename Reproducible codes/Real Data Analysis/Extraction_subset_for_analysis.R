rm(list=ls())
setwd("U:/BMASTER/Real Data Analysis")
library(ggplot2)

# Names of the variables, X = genera, Y = metabolite

X <- read.table("Data/X.csv", header = FALSE, sep = ",")    # p = 287 entries
Y <- read.table("Data/Y.csv", header = FALSE, sep = ",")    # q = 249 entries

X.names <- read.table("Data/X_names.csv", header = FALSE)    # p = 287 entries
Y.names <- read.table("Data/Y_names.csv", header = FALSE)    # q = 249 entries
p <- dim(X.names)[1]
q <- dim(Y.names)[1]

subset1 <- c("C00163_Propionate", "C00246_Butyrate",
                   "C00429_Dihydrouracil", "C00025_Glu", "C00086_Urea",
                   "C00042_Succinate", "C00431_5.Aminovalerate", 
                   "C00803_Valerate", "C00047_Lys", "C00041_Ala")
size.1 <- length(subset1)

# Extracting the important Y variables

Y.names.vec <- unlist(Y.names)
Y.pos.subset1 <- match(subset1, Y.names.vec)

Y.subset1 <- Y[, Y.pos.subset1]
Y.names.subset1 <- subset1

write.table(Y.subset1, "Data/Y.subset1.csv", sep = ",", col.names = FALSE, 
            row.names = FALSE)
write.table(Y.names.subset1, "Data/Y.names.subset1.csv", sep = ",", col.names = FALSE, 
            row.names = FALSE)
write.table(Y.pos.subset1, "Data/Y.indexes.subset1.csv", sep = ",", col.names = FALSE, 
            row.names = FALSE)


subset2 <- c("X_DCA", "C01921_Glycocholate", "C05122_Taurocholate",
             "C08262_Isovalerate", "C00407_Ile", "C00123_Leu", 
             "C00183_Val", "C00079_Phe", "C00082_Tyr", "C00065_Ser",
             "C00037_Gly")
size.2 <- length(subset2)

# Extracting the important Y variables

Y.names.vec <- unlist(Y.names)
Y.pos.subset2 <- match(subset2, Y.names.vec)

Y.subset2 <- Y[, Y.pos.subset2]
Y.names.subset2 <- subset2

write.table(Y.subset2, "Data/Y.subset2.csv", sep = ",", col.names = FALSE, 
            row.names = FALSE)
write.table(Y.names.subset2, "Data/Y.names.subset2.csv", sep = ",", col.names = FALSE, 
            row.names = FALSE)
write.table(Y.pos.subset2, "Data/Y.indexes.subset2.csv", sep = ",", col.names = FALSE, 
            row.names = FALSE)
