# Script to produce input data for clustering with BackSPIN and vizualisation with SPRING.
# Volker Hovestadt (hovestadt@broadinstitute.org), Jan 14, 2019
#
# From the methods:
#
# For all downstream analysis we required cells to have at least 1,000 UMIs (gene counts, indicative
# of the number of captured transcripts) mapping to at least 500 unique genes. We additionally
# excluded cells for which more than 20% of the gene counts reflected either mitochondrial genes or
# ribosomal RNAs, as these likely reflected poor quality cells. For downstream analyses, we normalized
# gene counts to a total of 10,000 for each cell.
# 
# BackSPIN employs a bi-clustering algorithm which iteratively splits both cells and genes, until a
# predetermined number of splits is reached. We selected the BackSPIN algorithm because it performs
# well when dealing with a relatively large number of cell populations. This is especially true for
# datasets in which some clusters are demarcated by a large number of differentially expressed genes
# (e.g. between the myeloid and lymphoid lineages), and others by relatively few genes (e.g. different
# populations within the myeloid lineage).
# 
# For clustering, we first determined the most variably expressed genes in the dataset. We performed a
# linear fit of the log-transformed average expression values and the log-transformed coefficients of
# variation (standard deviation divided by the average expression). Variably expressed genes were
# determined as genes associated with a residual larger than two times the standard deviation of all
# residuals. From these genes we excluded a set of genes that were associated with cell cycle (ASPM,
# CENPE, CENPF, DLGAP5, MKI67, NUSAP1, PCLAF, STMN1, TOP2A, TUBB). This yielded in the order of 1,000
# to 2,000 variably expressed genes depending on the set of cells. Expression values were
# log-transformed (after addition of 1) before performing BackSPIN clustering. We used default
# settings and a maximum splitting depth of 5. In the healthy BM data this yielded a final set of 31
# clusters.
# 
# For KNN visualization we calculated pairwise correlation coefficients between single cells (using the
# same set of variable genes as for the BackSPIN clustering). Then we constructed a graph by connected
# each cell to its five most highly correlated neighbors. This graph was visualized using SPRING, an
# interactive tool that uses force-directed graph drawing.
# 
# This script reads one or multiple digital expression matrices and produces a single CEF file used for
# BackSPIN clustering, as well as input files required for SPRING vizualisation (KNN graph).


options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)


### 0) define functions
# message("cutf()")
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# message("loadStarExpr()")
loadStarExpr <- function(files, names, min.umis=NULL, min.genes=NULL, cells=NULL) {
  for(i in seq(length(files))) {
    message(files[i])
    load(files[i])  # contains objects d and d.stats
    colnames(d) <- paste0(names[i], "_", colnames(d))
    
    if(!is.null(min.umis) & !is.null(min.genes)) d <- filterCells(d, min.umis=min.umis, min.genes=min.genes)
    if(length(cells)>0) d <- d[, intersect(colnames(d), cells), drop=FALSE]
    
    rownames(d.stats) <- paste0(names[i], "_", rownames(d.stats))
    d.stats <- d.stats[colnames(d), , drop=FALSE]   # some are not in d if 0 reads mapped to transcriptome
    if(i == 1) {
      D.all <- d
      D.stats.all <- d.stats
    } else {
      D.all <- cbind(D.all, d)
      D.stats.all <- rbind(D.stats.all, d.stats)
    }
    gc()
  }
  attr(D.all, "stats") <- D.stats.all
  D.all
}

# message("filterCells()")
filterCells <- function(D, min.umis=1000, min.genes=500) {
  D.umis <- colSums(D)   # transcripts
  D.genes <- colSums(D>0)  # genes, observed at least once
  plot(D.umis, D.genes, xlim=c(0, max(D.umis)), ylim=c(0, max(D.genes)), pch=".", col=ifelse(D.umis>=min.umis & D.genes>=min.genes, "black", "grey"))
  abline(h=min.genes)
  abline(v=min.umis)
  D[, D.genes >= min.genes & D.umis >= min.umis]
}

# message("filterGenes()")
filterGenes <- function(D, g, D.umis=NULL, max.g=1, xlab="total", ylab="fraction", ...) {
  g <- intersect(g, rownames(D))
  message(length(g), " genes filtered")
  if(is.null(D.umis)) D.umis <- colSums(D)  # provide to use same number in all gene sets
  D.umis <- D.umis[colnames(D)]
  D.g <- colSums(D[g, ])
  plot(D.umis, D.g/D.umis, log="x", ylim=c(0, 1), pch=".", col=ifelse(D.g/D.umis > max.g, "red", "black"), xlab=xlab, ylab=ylab, ...)
  abline(h=max.g)
  message(sum(D.g/D.umis > max.g), " cells filtered")
  D[!is.element(rownames(D), g), D.g/D.umis <= max.g]  # filter genes and cells above threshold
}

# message("variableGenes()")
variableGenes <- function(E, E.min=0.01, E.sd=2) {
  E.mean <- rowMeans(E)  # average expression per gene
  E.cv <- apply(E, 1, sd) / E.mean  # coefficient of variation
  
  E.mean <- sort(E.mean[!is.na(E.cv)])  # remove genes that are not expressed at all
  E.cv <- E.cv[names(E.mean)]
  
  plot(E.mean, E.cv, pch=".", xlab="Mean", ylab="CV")
  
  plot(E.mean, E.cv, log="xy", pch=".", xlab="Mean", ylab="CV")
  abline(v=E.min, col="blue")
  
  E.lm <- lm(y~x, subset = E.mean>=E.min, data.frame(x=log10(E.mean), y=log10(E.cv)))  # linear fit
  E.cv.predict <- predict(E.lm, data.frame(x=log10(E.mean)))
  lines(E.mean, 10^E.cv.predict, col="red")
  
  E.cv.residue <- log10(E.cv) - E.cv.predict
  plot(E.mean, E.cv.residue, log="x", pch=".", xlab="Mean", ylab="Residue")
  abline(v=E.min, col="green")
  abline(h=0, col="red")
  abline(h=sd(E.cv.residue[E.mean>=E.min])*c(1, 2, 3), col="red", lty=2)
  axis(4, sd(E.cv.residue[E.mean>=E.min])*c(1, 2, 3), las=2)
  
  plot(seq(-0.25+0.005, 1-0.005, 0.01), as.vector(table(cut(E.cv.residue[E.mean>=E.min], seq(-0.25, 1, 0.01)))), xlab="Residue", ylab="Frequency", type="b")
  abline(v=0, col="red")
  abline(v=sd(E.cv.residue[E.mean>=E.min])*c(1, 2, 3), col="red", lty=2)
  axis(3, sd(E.cv.residue[E.mean>=E.min])*c(1, 2, 3), las=2)
  
  plot(ecdf(E.cv.residue[E.mean>=E.min]), main="ECDF")
  qqnorm(E.cv.residue[E.mean>=E.min])
  
  m <- sort((E.cv.residue)[E.mean>=E.min], decreasing = TRUE)
  attr(m, which = "sd") <- sd(E.cv.residue[E.mean>=E.min])
  
  plot(E.mean, E.cv, log="xy", pch=".", col=ifelse(names(E.cv) %in% names(which(m > E.sd*attr(m, which = "sd"))), "red", "grey"), xlab="Mean", ylab="CV")
  abline(v=E.min, col="blue")
  lines(E.mean, 10^E.cv.predict, col="red")
  
  m  # return
}

# message("write.cef()")
write.cef <- function(x, file, header=c(), attr.col=data.frame(), attr.row=data.frame(), verbose=TRUE) {
  if(length(dim(x)) != 2) {
    if(is.element("header", names(x))) header <- x$header   # use information from list if input is list
    if(is.element("attr.col", names(x))) attr.col <- x$attr.col
    if(is.element("attr.row", names(x))) attr.row <- x$attr.row
    if(!is.element("data", names(x))) stop("need data in list")
    x <- x$data
  }
  
  if(verbose) message("header count: ", length(header))
  if(verbose) message("row attribute count: ", ncol(attr.row))
  if(verbose) message("column attribute count: ", ncol(attr.col))
  if(verbose) message("row count: ", nrow(x))
  if(verbose) message("column count: ", ncol(x))
  if(verbose) message("Flags value: ", 0)
  
  if(verbose) message("write FIRST LINE")
  write.table(t(c("CEF", length(header), ncol(attr.row), ncol(attr.col), nrow(x), ncol(x), 0, rep("", ncol(attr.row)+1+ncol(x)-7))), file=file, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  if(verbose) message("write HEADER")
  write.table(data.frame(names(header), header, rep(list(""), ncol(attr.row)+1+ncol(x)-2)), file=file, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  if(verbose) message("write COLUMN ATTRIBUTES")
  write.table(data.frame(rep(list(""), ncol(attr.row)), colnames(attr.col), t(attr.col)), file=file, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  if(verbose) message("write ROW ATTRIBUTES & DATA")
  write.table(t(c(names(attr.row), rep("", 1+ncol(x)))), file=file, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(data.frame(attr.row, "", x), file=file, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}


### 1) Load data, and perform initial quality filtering

# example sample sheet (sample name, RData file(s) containing the digital expression matrix, as produced by star2count script):
# BM1	BM1.star.expr.RData
# BM2	BM2.star.expr.RData
# BM3	BM3.star.expr.RData
# BM4	BM4.star.expr.RData
# BM5-34p	BM5-34p.star.expr.RData
# BM5-34p38n	BM5-34p38n.star.expr.RData
# .. add cell IDs below if selecting only a subset of cells

a <- commandArgs(trailingOnly=TRUE)[1]  # sample sheet
b <- read.table(a, fill=TRUE, sep = "\t", col.names = paste0("V", seq(9)))
dim(b)

m <- grepl(".expr.RData$", b$V2)
message("loading data from ", sum(m), " files")

# subset cells
m2 <- b[!m, "V1"]
if(length(m2)) message("only use ", length(m2), " cells listed in sample sheet")

pdf(paste0(a, ".cells.pdf"), width = 6, height = 6)
par(mar=c(4, 4, 4, 4))

D <- loadStarExpr(b[m, "V2"], b[m, "V1"], min.umis=1000, min.genes=500, cells=m2)  # change QC thresholds here
D.stats <- attr(D, "stats")

dev.off()


### 2) Filter genes and cells

pdf(paste0(a, ".filter.pdf"), width = 6, height = 6)
par(mar=c(4, 4, 4, 4))

id2chr <- read.table("hg38.RefSeq_refGene.gtf.transcripts")   # mapping of gene IDs to chromosomes
gene2id <- read.table("hg38.refMrna.gene.chrM.ercc.fa.names", fill = NA)   # mapping of gene names to gene IDs

id2chr <- rbind(id2chr, data.frame(V1="chrM", V2=substr(gene2id[grepl("^MT-", gene2id$V1), "V2"], 1, 12)))    # add chrM
id2chr <- id2chr[nchar(id2chr$V1)<=5, ]   # keep only main chromosomes

gene2chr <- lapply(lapply(split(split(id2chr$V1, id2chr$V2)[cutf(gene2id$V2, 1, "\\.")], gene2id$V1), unlist), unique)
table(lengths(gene2chr))   # number of chromosomes per gene

none.filter <- names(gene2chr[lengths(gene2chr) == 0])  # genes not mapping to the genome (ERCC, etc)
multi.filter <- names(gene2chr[lengths(gene2chr) >= 2])  # genes mapping to different chromosomes (mostly X and Y)
xy.filter <- names(which(sapply(lapply(gene2chr, is.element, c("chrX", "chrY")), any)))   # genes mapping to X and Y

mt <- c("MTRNR2", "SNORA109", rownames(D)[grepl("^MT-", rownames(D))])   # mitochondrial genes, SNORA109 is highly correlated
rr <- c("RNA5S", "RNA18S", "RNA28S")   # ribosomal RNAs
rp <- sort(rownames(D)[grepl("^RPS", rownames(D)) | grepl("^RPL", rownames(D))])   # ribosomal proteins

D.umi <- colSums(D)  # calculate once and use for all filtering steps
# only keep genes mapping to a single chromosome (excluding X and Y)
D <- filterGenes(D, none.filter, D.umi, main="none")
D <- filterGenes(D, multi.filter, D.umi, main="multi")
D <- filterGenes(D, xy.filter, D.umi, main="xy")  # use this when comparing samples from both M and F
# exclude cells with high mitochondrial or rRNA transcripts
D <- filterGenes(D, mt, D.umi, 0.2, main="mitochondrial genes")  # cells with >20% mitochondrial RNAs are filtered. change QC thresholds here
D <- filterGenes(D, rr, D.umi, 0.2, main="ribosomal rnas")  # cells with >20% ribosomal RNAs are filtered
D <- filterGenes(D, rp, D.umi, main="ribosomal proteins")  

dev.off()


### 3) Determine variable genes

# normalize
E <- t(t(D)/colSums(D)) *10000   # UMIs per 10,000

# variable genes
pdf(paste0(a, ".var.pdf"), width = 6, height = 6)
par(mar=c(4, 4, 4, 4))

E.var <- variableGenes(E)  # default parameters, return the most variable genes
m <- names(which(E.var > attr(E.var, "sd")))

dev.off()

# exclude genes associated with cell cycle
sign.rev <- c("TOP2A", "CENPF", "MKI67", "NUSAP1", "CENPE", "DLGAP5", "STMN1", "TUBB", "PCLAF", "ASPM")
m <- setdiff(m, sign.rev)


### 4) Write RData file of full dataset 
save(list=c("D", "D.stats", "E", "m"), file  = paste0(a, ".RData"))


### 5) Write output CEF file used for BackSPIN clustering
dim(E[m, ])

set.seed(11)
s <- sample(colnames(E))   # set seed, ordering affects BackSPIN results
write.cef(x = log(E[m, s]+1),  # log values
          attr.col = data.frame(cellID=s, sample=cutf(s, 1, "_"), D.stats[s, ]),
          attr.row = data.frame(GeneType="mRNA", Gene=m, GeneGroup=1),
          header = c(Genome="hg38"),
          file  = paste0(a, ".cef"))


### 6) Write output for SPRING visualization
library(jsonlite)

## generate KNN graph
# log
X <- log(E[sort(m), ]+1)
dim(X)

# cor
X.cor <- cor(X)
dim(X.cor)

# knn
k <- 5
X.knn.list <- lapply(seq(nrow(X.cor)), function(s) {
  data.frame(source=s, target=order(X.cor[s, ], decreasing=TRUE)[2:(k+1)], distance=1)
})
X.knn <- do.call(rbind, X.knn.list)
dim(X.knn)

# make undirected
X.knn[X.knn$source > X.knn$target, c("source", "target")] <- X.knn[X.knn$source > X.knn$target, c("target", "source")]
dim(X.knn)

X.knn <- unique(X.knn)
dim(X.knn)

## output json file
bm_json <- list()
bm_json$nodes <- data.frame(name=colnames(X), number=seq(ncol(X))-1)
bm_json$links <- X.knn-1

dir.create(sub(".txt$", ".SPRING", a))
write(toJSON(bm_json, pretty=TRUE), file=sub(".txt$", ".SPRING/graph_data.json", a))

## create gene_colors files  [MEAN, STANDARD DEVIATION, MIN, MAX, 99-PERCENTILE]
X.stats <- lapply(apply(X, 1, function(x) list(mean(x), sd(x), min(x), max(x), quantile(x, 0.99, names = FALSE))), unlist)
write(toJSON(X.stats, pretty=TRUE), file=sub(".txt$", ".SPRING/color_stats.json", a))

dir.create(sub(".txt$", ".SPRING/gene_colors", a))
X.rowsplit <- unname(split(rownames(X), cut(seq(nrow(X)), breaks = 51, labels = FALSE)))
for(i in seq(length(X.rowsplit))) {
  message(i)
  write.table(X[X.rowsplit[[i]], ], sep=",", quote=FALSE, col.names=FALSE, file=paste0(sub(".txt$", ".SPRING/gene_colors/", a), "color_data_all_genes-", i-1, ".csv"))
}

