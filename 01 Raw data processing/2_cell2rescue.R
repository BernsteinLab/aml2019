# Script to perform the second step of the raw data processing.
# Volker Hovestadt (hovestadt@broadinstitute.org), Jan 14, 2019
#
# From the methods:
#
# For some of the sequencing runs we noticed a higher number of reads that were excluded because the
# library barcode was not detected accurately. We rescued these reads if their CB matched uniquely to
# one of the libraries that were sequenced together in the respective run.
# 
# This script is run once for all sample per run. It generates lists of CBs for every sample
# annotated by whether reads for this CB should be rescued (from the Undetermined files).


options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)


library(data.table)

cutf <- function(x, f=1, d="/", ...) sapply(strsplit(x, d), function(i) i[f], ...)
na0 <- function(x) ifelse(is.na(x), 0, x)

a <- commandArgs(trailingOnly=TRUE)[1]  # directory in which files containing cell barcodes are located

lf.all <- list.files(a, pattern = "all.txt$", full.names = "TRUE")
lf.good <- list.files(a, pattern = "good.txt$", full.names = "TRUE")
# message(length(lf.all), " *all.txt files found")
# message(length(lf.good), " *good.txt files found")

d.all <- lapply(lf.all, fread)
d.good <- lapply(lf.good, fread)

names(d.all) <- cutf(lf.all, 2, "\\.")
names(d.good) <- cutf(lf.good, 2, "\\.")

for(i in names(d.all)) d.all[[i]] <- rbind(d.all[[i]], d.all[[i]][, sum(V2), by=substr(V1, 1, 11)], use.names=FALSE)
for(i in names(d.good)) d.good[[i]] <- d.good[[i]][, sum(V2), by=sub("N$", "", V1)]

for(i in names(d.all)) d.all[[i]]$sample <- paste0(i, ".all")
for(i in names(d.good)) d.good[[i]]$sample <- paste0(i, ".good")

d.dcast <- as.data.frame(dcast(rbindlist(c(d.good, d.all)), sub~sample, value.var="V1", fill = 0))
rownames(d.dcast) <- d.dcast$sub
d.dcast$sub <- NULL
# dim(d.dcast)

d.dcast <- d.dcast[apply(d.dcast[, grep(".good$", colnames(d.dcast)), drop=FALSE] > 0, 1, any), ]
# dim(d.dcast)


for(i in names(d.good)) {
  # i <- "AML150921A-1"

  pdf(sub(".txt", ".rescue.pdf", lf.good[match(i, cutf(lf.good, 2, "\\."))]), width = 6, height = 6)
  par(mar=c(4, 4, 4, 4))

  y <- d.dcast[, paste0(i, ".good")]
  names(y) <- rownames(d.dcast)
  y <- sort(y[y>0], decreasing = TRUE)
  plot(y, log="xy", xlab=NA, ylab="Reads", main=i, type="l", las=3)
  axis(3, c(500, 1000, 2000, 5000), y[c(500, 1000, 2000, 5000)], las=2)
  abline(v=c(500, 1000, 2000, 5000), col="grey")

  if(length(setdiff(names(d.all), i))) {
    x <- rowSums(d.dcast[names(y), paste0(setdiff(names(d.all), i), ".all"), drop=FALSE])
  } else {
    x <- rep(0, length(y))
  }
  plot(x+1, y, log="xy", xlab="Reads+1 (other samples)", ylab="Reads", main=i,
       pch=16, cex=1/2, col=ifelse(x/y<0.02, "black", "grey"))
  abline(0, 1)
  abline(log10(50), 1, lty=2)  # 1/50 = 2%

  dev.off()

  write.table(data.table(ifelse(nchar(names(y))==11, paste0(names(y), "N"), names(y)), y, x/y<0.02),
              file = sub(".txt", ".rescue.txt", lf.good[match(i, cutf(lf.good, 2, "\\."))]),
              sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}

