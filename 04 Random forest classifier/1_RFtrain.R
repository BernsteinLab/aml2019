# Training of the random forest classifier.
# Volker Hovestadt (hovestadt@broadinstitute.org), Jan 14, 2019
#
# From the methods:
#
# The Random forest algorithm is a machine learning approach that uses a large number of binary
# decision trees that are learned from random subsets of a training set. These trees (the forest) can
# then be applied to a given sample to generate a class probability that reflects its similarity to a
# given class of the training set. If a single class prediction is required, the class with the
# highest probability score is used (majority vote). Random forest classifiers are particularly well
# suited if the dataset contains many different classes, many samples and many features. In our case
# samples represent single-cell expression profiles, features represent genes, and classes represent
# different cell types. For our analysis we used the randomForest R package version 4.6-14.
# 
# We used Random forest-based classification for two different purposes: To predict similarity of
# single cells to the 15 different cell types detected in healthy BM (classifier 1), and to predict if
# a single cell from a tumor sample is malignant or normal (classifier 2). To train the first
# classifier, we first performed a feature selection step to select the most informative genes from
# all 14,554 expressed genes in the dataset (average expression > 0.01). Feature selection was
# performed by training an “outer” random forest classifier on all expressed genes. We trained 1,000
# trees, using a random subset of 50 cells from each cell type for each tree. Based on the reported
# overall gene importance in the “outer” classifier, we then selected only the 1,000 most informative
# genes for training of the “inner” classifier. The reported out-of-bag (OOB) error (i.e.
# misclassification error of cells that were not used for learning of a given tree) was 20% lower for
# the “inner” than for the “outer” classifier, justifying the use of an initial feature selection
# step. The “inner” classifier was further evaluated by performing 5-fold cross-validation by
# splitting the training dataset into five equally sized parts. In each iteration of the
# cross-validation, four of these parts were used to generate a classifier that was then used for
# predicting class probabilities of the remaining part. Results of the cross-validation are provided
# in Figure S3A.
# 
# The second classifier is used for determining if a cell for which we did not detect a mutant
# transcript is malignant or normal, based on its similarity to normal and malignant cells (i.e. cells
# from healthy BM and HSC to myeloid-like cells from tumor samples for which we detected mutant
# transcripts). We first attempted to use a classifier that distinguishes between just these two
# classes. However, we achieved much better results by using all 15 normal and six malignant cell
# types in a combined training set (21 classes), presumably because a malignant monocyte-like cell is
# more similar to a normal monocyte than to a malignant HSC-like cell. For malignant cells we used
# cell type annotations as predicted by the first classifier, with the following exceptions: to have
# at least 65 HSC-like cells for each malignant class (required for having >50 cells for 5-fold
# cross-validation), we reclassified 23 cells initially classified as progenitor-like with highest
# prediction scores for the HSC cell type as HSC-like cells. We also reclassified 29 cells that were
# initially classified as early Erythroid progenitors as progenitor-like cells, if their prediction
# score for the Progenitor cell type was higher than their prediction score for the late Erythroid
# cell type. The second classifier was then generated using the combined training set of 21 classes
# and the same parameters as for the first classifier. The second classifier reached 95.2% sensitivity
# and 99.7% specificity in distinguishing malignant from normal cells, as measured by 5-fold
# cross-validation. Results of the 5-fold cross-validation are provided in Figure S3E.
#
# This script saves two RData files containing the the first and second random forest classifiers.


options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

library(randomForest)


### 1) Prepare
# load expression data
load("BM_6915cell.RData", verbose = TRUE)  # as generated using the prepareBackSpinKNN script
X <- log(E+1)  # filtered gene set
dim(X)

# load BackSPIN cell type annotation
a <- read.table("BM_6915cells.BackSPIN.txt", row.names=1, header=TRUE)  # cell annotations available on GEO
# head(a)
#                   cluster
# BM4_AAAAAACATACG earlyEry
# BM4_AAAAACTTCTAT earlyEry
# BM4_AAAAATCTGAGA earlyEry
# BM4_AAAACACGAATA  ProMono
# BM4_AAAACGTCGCTN       NK
# BM4_AAAAGACTGTAT earlyEry
# ...

# check if annotations available for all cells
setdiff(rownames(a), colnames(X))
setdiff(colnames(X), rownames(a))

unique(a$cluster)
sort(table(a$cluster))


### 2) Train RF classifier #1 (normal BM data only)

pdf("RF1.pdf", width = 6, height = 6)   # updated for cv heatmap
par(mar=c(4, 4, 4, 4))

## outer classifier (variable selection)
rf.all.expr <- names(which(rowMeans(exp(X)-1) > 0.01))
length(rf.all.expr)

set.seed(123)
rf.all.outer <- randomForest(x = t(X[rf.all.expr, rownames(a)]),
                             y = factor(a$cluster, unique(a$cluster)),
                             sampsize = rep(50, length(unique(a$cluster))),
                             ntree = 1000,
                             do.trace=TRUE)

plot(sort(rf.all.outer$importance[, 1], decreasing = TRUE), type="l")
abline(v=1000)

cols <- apply(colorRamp(c("#FFFFFF", "red"))(seq(0, 1, 0.001)), 1, function(rr) rgb(rr[1], rr[2], rr[3], maxColorValue = 255))
image(as.matrix(seq(0, 1, 0.001)), col=cols, main="scale")

z <- rf.all.outer$confusion[, rownames(rf.all.outer$confusion)] / rowSums(rf.all.outer$confusion[, rownames(rf.all.outer$confusion)])
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM, outer classifier, OOB confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.all.outer$confusion[, rownames(rf.all.outer$confusion)]==0, NA, rf.all.outer$confusion[, rownames(rf.all.outer$confusion)])))

1 - sum(diag(rf.all.outer$confusion[, rownames(rf.all.outer$confusion)])) / sum(rf.all.outer$confusion[, rownames(rf.all.outer$confusion)])

# select 1000 most informative genes (genes chosen most often in the outer classifier)
rf.all.outer.used1k <- names(sort(table(rownames(rf.all.outer$importance)[rf.all.outer$forest$bestvar[rf.all.outer$forest$bestvar!=0]]), decreasing = TRUE)[1:1000])
plot(rf.all.outer$importance[rf.all.outer.used1k, 1])  # correlates to importance measure

## inner classifier
set.seed(123)
rf.all.inner <- randomForest(x = t(X[rf.all.outer.used1k, rownames(a)]),
                             y = factor(a$cluster, unique(a$cluster)),
                             sampsize = rep(50, length(unique(a$cluster))),
                             ntree = 1000,
                             do.trace=TRUE)

z <- rf.all.inner$confusion[, rownames(rf.all.inner$confusion)] / rowSums(rf.all.inner$confusion[, rownames(rf.all.inner$confusion)])
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM, inner classifier, OOB confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.all.inner$confusion[, rownames(rf.all.inner$confusion)]==0, NA, rf.all.inner$confusion[, rownames(rf.all.inner$confusion)])))

1 - sum(diag(rf.all.inner$confusion[, rownames(rf.all.inner$confusion)])) / sum(rf.all.inner$confusion[, rownames(rf.all.inner$confusion)])

## inner classifier, 5-fold cross validation
y <- factor(a$cluster, unique(a$cluster))
names(y) <- rownames(a)
x <- rownames(a)
cv <- split(x, rep(1:5, 1E6)[1:length(y)])
lengths(cv)

rf.all.cv <- lapply(cv, function(y2) {
  set.seed(123)
  randomForest(x = t(X[rf.all.outer.used1k, setdiff(x, y2)]),
               y = y[setdiff(x, y2)],
               sampsize = rep(50, length(unique(y))),
               ntree = 1000,
               do.trace=TRUE)
})

# predict the sets that were not used for training
rf.all.cv.predict.prob <- lapply(rf.all.cv, function(rf) {
  predict(rf, t(X[rownames(rf$importance), setdiff(colnames(X), names(rf$y))]), type = "prob")
})
rf.all.cv.predict <- lapply(rf.all.cv.predict.prob, function(x) {
  y <- factor(colnames(x)[apply(x, 1, which.max)], colnames(x))  # max
  names(y) <- rownames(x)
  y
})

# generate CV confusion matrix
rf.all.cv.predict.matrix <- table(unlist(rf.all.cv.predict), y[unlist(lapply(rf.all.cv.predict, names))])
rf.all.cv.predict.matrix <- rf.all.cv.predict.matrix[, rownames(rf.all.cv.predict.matrix)]
sum(rf.all.cv.predict.matrix)
write.table(rf.all.cv.predict.matrix, sep="\t", quote=FALSE, file="RF1.cv.txt")

z <- t(t(rf.all.cv.predict.matrix) / colSums(rf.all.cv.predict.matrix))
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM, inner classifier, CV confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.all.cv.predict.matrix==0, NA, rf.all.cv.predict.matrix)))

1 - sum(diag(rf.all.cv.predict.matrix)) / sum(rf.all.cv.predict.matrix)

dev.off()

## export RData object for RF1
save(list=c("rf.all.outer", "rf.all.inner", "rf.all.cv"), file="RF1.RData")


### 3) Classify cells from tumor samples with detected mutations
# load tumor data (cell mutation status available on GEO)
load("AMLmut_923cells.RData", verbose = TRUE)    # as generated using the prepareBackSpinKNN script
Y <- log(E+1)

## classify
pdf("RF1.AMLclassify.pdf", width = 6, height = 6)
par(mar=c(4, 4, 4, 4))

rf.all.inner.mut <- predict(rf.all.inner, newdata = t(Y[rownames(rf.all.inner$importance), ]), type = "prob")
rf.all.inner.mut.max <- factor(colnames(rf.all.inner.mut)[apply(rf.all.inner.mut, 1, which.max)], colnames(rf.all.inner.mut))
names(rf.all.inner.mut.max) <- rownames(rf.all.inner.mut)
barplot(table(rf.all.inner.mut.max), las=2)

# reclassify Prog with highest score for HSC to HSC to have 65 cells in total
rf.all.inner.mut.prog2hsc <- names(sort(rf.all.inner.mut[rf.all.inner.mut.max == "Prog", "HSC"], decreasing = TRUE)[1:(65-table(rf.all.inner.mut.max)["HSC"])])
plot(rf.all.inner.mut[, "Prog"], rf.all.inner.mut[, "HSC"], col=ifelse(rownames(rf.all.inner.mut) %in% rf.all.inner.mut.prog2hsc, "red", "black"), xlim=c(0, 0.5), ylim=c(0, 0.5), pch=16)
abline(0, 1)

# reclassify earlyEry with higher score for Prog than for lateEry to Prog
rf.all.inner.mut.ery2prog <- names(which(rf.all.inner.mut[rf.all.inner.mut.max=="earlyEry", "Prog"] > rf.all.inner.mut[rf.all.inner.mut.max=="earlyEry", "lateEry"]))
plot(rf.all.inner.mut[, "Prog"], rf.all.inner.mut[, "earlyEry"], xlim=c(0, 0.5), ylim=c(0, 0.5), col=ifelse(rownames(rf.all.inner.mut) %in% rf.all.inner.mut.ery2prog, "red", "black"), pch=16)
plot(rf.all.inner.mut[, "Prog"], rf.all.inner.mut[, "lateEry"], xlim=c(0, 0.5), ylim=c(0, 0.5), col=ifelse(rownames(rf.all.inner.mut) %in% rf.all.inner.mut.ery2prog, "red", "black"), pch=16)

# create list
a.mut <- split(names(rf.all.inner.mut.max), rf.all.inner.mut.max)[c("HSC", "Prog", "GMP", "ProMono", "Mono", "cDC")]
lengths(a.mut)
a.mut$HSC <- c(a.mut$HSC, rf.all.inner.mut.prog2hsc)
a.mut$Prog <- setdiff(a.mut$Prog, rf.all.inner.mut.prog2hsc)
a.mut$Prog <- c(a.mut$Prog, rf.all.inner.mut.ery2prog)
lengths(a.mut)

dev.off()


### 4) Train RF classifier #2 (normal BM and AML cells with mutations)
pdf("RF2.pdf", width = 6, height = 6)
par(mar=c(4, 4, 4, 4))

a2 <- c(split(rownames(a), a$cluster)[unique(a$cluster)], Tumor=a.mut)
lengths(a2)

XY <- cbind(X, Y[rownames(X), ])[, unlist(a2)]
dim(XY)

## outer classifier (variable selection)
rf.tumor.expr <- names(which(rowMeans(XY[, unlist(a2)]) >= 0.01))
length(rf.tumor.expr)

set.seed(123)
rf.tumor.outer <- randomForest(x = t(XY[rf.tumor.expr, unlist(a2)]),
                               y = factor(rep(names(a2), lengths(a2)), names(a2)),
                               sampsize = rep(50, length(a2)),
                               ntree = 1000,
                               do.trace=TRUE)

z <- rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)] / rowSums(rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)])
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM+mut, outer classifier, OOB confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)]==0, NA, rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)])))

1 - sum(diag(rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)])) / sum(rf.tumor.outer$confusion[, rownames(rf.tumor.outer$confusion)])

rf.tumor.outer.used1k <- names(sort(table(rownames(rf.tumor.outer$importance)[rf.tumor.outer$forest$bestvar[rf.tumor.outer$forest$bestvar!=0]]), decreasing = TRUE)[1:1000])
plot(rf.tumor.outer$importance[rf.tumor.outer.used1k, 1])

## inner classifer
set.seed(123)
rf.tumor.inner <- randomForest(x = t(XY[rf.tumor.outer.used1k, unlist(a2)]),
                               y = factor(rep(names(a2), lengths(a2)), names(a2)),
                               sampsize = rep(50, length(a2)),
                               ntree = 1000,
                               do.trace=TRUE)

z <- rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)] / rowSums(rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)])
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM+mut, inner classifier, OOB confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)]==0, NA, rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)])))

1 - sum(diag(rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)])) / sum(rf.tumor.inner$confusion[, rownames(rf.tumor.inner$confusion)])

## inner, 5-fold cv
y <- factor(rep(names(a2), lengths(a2)), names(a2))
names(y) <- unlist(a2)
x <- unlist(a2)
cv <- split(x, rep(1:5, 1E6)[1:length(y)])

rf.tumor.cv <- lapply(cv, function(y2) {
  set.seed(123)
  randomForest(x = t(XY[rf.tumor.outer.used1k, setdiff(x, y2)]),
               y = y[setdiff(x, y2)],
               sampsize = rep(50, length(unique(y))),
               ntree = 1000,
               do.trace=TRUE)
})

rf.tumor.cv.predict.prob <- lapply(rf.tumor.cv, function(rf) {
  predict(rf, t(XY[rownames(rf$importance), setdiff(colnames(XY), names(rf$y))]), type = "prob")
})
rf.tumor.cv.predict <- lapply(rf.tumor.cv.predict.prob, function(x) {
  y <- factor(colnames(x)[apply(x, 1, which.max)], colnames(x))
  names(y) <- rownames(x)
  y
})

rf.tumor.cv.predict.matrix <- table(unlist(rf.tumor.cv.predict), y[unlist(lapply(rf.tumor.cv.predict, names))])
rf.tumor.cv.predict.matrix <- rf.tumor.cv.predict.matrix[, rownames(rf.tumor.cv.predict.matrix)]
sum(rf.tumor.cv.predict.matrix)
write.table(rf.all.cv.predict.matrix, sep="\t", quote=FALSE, file="RF2.cv.txt")

z <- t(t(rf.tumor.cv.predict.matrix) / colSums(rf.tumor.cv.predict.matrix))
heatmap(z, Rowv=NA, Colv=NA, scale="n", zlim=c(0, 1), col=cols, main="BM+mut, inner classifier, CV confusion matrix",
        add.expr = text(rep(1:ncol(z), each=nrow(z)), rep(1:ncol(z), nrow(z)), ifelse(rf.tumor.cv.predict.matrix==0, NA, rf.tumor.cv.predict.matrix)))

1 - sum(diag(rf.tumor.cv.predict.matrix)) / sum(rf.tumor.cv.predict.matrix)

dev.off()


## export RData object for RF2
save(list=c("rf.tumor.outer", "rf.tumor.inner", "rf.tumor.cv"), file="RF2.RData")



