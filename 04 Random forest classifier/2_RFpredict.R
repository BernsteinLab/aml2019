# Using the random forest classifier to predict cell types in single cells.
# Volker Hovestadt (hovestadt@broadinstitute.org), Jan 14, 2019
#
# From the methods:
#
# When applying both classifiers to single cells from tumor samples, we first determined from the
# second classifier if the prediction score was highest for a malignant or normal cell type. If a cell
# was classified as malignant, we then used the highest prediction score of the HSC to myeloid cell
# types (HSC, progenitor, GMP, promonocyte, monocyte, cDC) from the first classifier for cell type
# assignment. If a cell was classified as normal, we used the highest prediction score from the first
# classifier.
#
# This applies the first and second random forest classifier (as generated using the RFtrain script)
# to single cells from any input sample. The output text file contains prediction scores for each cell
# type, as well as the highest scoring cell type.


options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(randomForest)


## load classifier
load("RF1.RData", verbose = TRUE)   # as generated using the RFtrain script
load("RF2.RData", verbose = TRUE)


## load data
f <- commandArgs(trailingOnly=TRUE)[1]  # input RData file (as generated using the prepareBackSpinKNN script)
load(f, verbose = TRUE)

X <- log(E+1)


## predict, within seconds
# RF1
X.predict.all <- predict(rf.all.inner, newdata = t(X[rownames(rf.all.inner$importance), ]), type="prob")
X.predict.all.max <- apply(X.predict.all, 1, max)
X.predict.all.class <- colnames(X.predict.all)[apply(X.predict.all, 1, which.max)]
table(X.predict.all.class)

# RF2
X.predict.tumor <- predict(rf.tumor.inner, newdata = t(X[rownames(rf.tumor.inner$importance), ]), type="prob")
X.predict.tumor.max <- apply(X.predict.tumor, 1, max)
X.predict.tumor.class <- colnames(X.predict.tumor)[apply(X.predict.tumor, 1, which.max)]
table(X.predict.tumor.class)


## write output
d.out <- data.frame(D.stats[names(X.predict.tumor.max), ],
                    RF1.class=X.predict.all.class, RF1.score=X.predict.all.max, RF1=X.predict.all, 
                    RF2.class=X.predict.tumor.class, RF2.score=X.predict.tumor.max, RF2=X.predict.tumor)
write.table(d.out, file=paste0(f, ".predict.txt"), quote = FALSE, sep = "\t")

save(list=c("X.predict.all", "X.predict.tumor"), file=paste0(f, ".predict.RData"))



