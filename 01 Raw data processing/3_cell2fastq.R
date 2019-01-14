# Script to perform the third step of the raw data processing.
# Volker Hovestadt (hovestadt@broadinstitute.org), Jan 14, 2019
#
# From the methods:
#
# The final list of CBs was then used to generate a fastq file containing the Read 2 sequences of the
# remaining cells. The library barcode, the CB, and the UMI were appended to the read identifier. Resulting
# fastq files for each sample were deposited in GEO.
# 
# This script is run for every sample in each run. A single sample can have multiple library barcodes
# in the sample sheet. It generates a single fastq file for each sample.


options(max.print = 500)

options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressMessages(library(data.table))
suppressMessages(library(ShortRead))
suppressMessages(library(stringi))


### Functions
cutf <- function(x, f=1, d="/", ...) sapply(strsplit(x, d), function(i) i[f], ...)

### Arguments
d <- commandArgs(trailingOnly=TRUE)[1]  # one or multiple directories containing fastq files, not searched recursively. add location of Undetermined files if reads should be rescued
l <- commandArgs(trailingOnly=TRUE)[2]  # sample sheet
s <- commandArgs(trailingOnly=TRUE)[3]  # sample name, should be in sample sheet, will be used for output files
b <- commandArgs(trailingOnly=TRUE)[4]  # list of cells, including if it should be rescued
u <- ifelse(length(commandArgs(trailingOnly=TRUE)) >= 5, commandArgs(trailingOnly=TRUE)[5], "")  # one or multiple directories containing undetermined fastq files

# sample sheet
a <- read.table(l)
if(nrow(a) == 0) stop("no entries in ", l)

# regular read files
lf <- unique(sort(unlist(lapply(strsplit(d, ","), list.files, pattern="_R1_", full.names=TRUE, recursive=FALSE))))  # bcl2fastq lib barcode split, R1 only
lf <- lf[grepl(".fastq.gz$", lf)]
if(length(lf) == 0) stop("no fastq.gz files in ", d)

# subset for single sample
lf <- lf[apply(do.call(cbind, lapply(paste0(sub("N700.", "N700-", a[a$V3 == s, "V2"]), "_S"), grepl, x=lf)), 1, any)]  # subset
if(length(lf) == 0) stop("no sample match in ", d)

# cells
x <- read.table(b)
if(grepl("00.txt$", b)) x <- x[x$V2 >= as.numeric(head(tail(strsplit(b, "\\.")[[1]], 2), 1)), ]
if(nrow(x) == 0) stop("no entries in ", b)
x$V1 <- sub("N$", "", x$V1)

# undetermined
lfu <- unique(sort(unlist(lapply(strsplit(u, ","), list.files, pattern="_R1_", full.names=TRUE, recursive=FALSE))))  # undetermined, R1 only
lfu <- lfu[grepl(".fastq.gz$", lfu)]
lf <- c(lf, lfu)


### process fastq files
#dir.create(paste0(sub(".txt$", "", l), ".", s))

n <- list()  # empty list to store number of reads

message("\n### ", Sys.time(), " - reading R1 and R2 sequences, searching for ", nrow(x), " cells")
for(f1 in lf) {
  f2 <- sub("_R1_", "_R2_", f1)
  message("file ", match(f1, lf), "/", length(lf), ": ", f1, " ", appendLF = FALSE)
  n[[f1]] <- c(0, 0)

  strm1 <- FastqStreamer(f1, n=1E7)  # 1M reads by default
  strm2 <- FastqStreamer(f2, n=1E7)
  repeat {
    message("*", appendLF = FALSE)
    fq1 <- yield(strm1)
    if(length(fq1) == 0) break
    fq2 <- yield(strm2)

    # match to expected cell barcodes
    if(is.element(f1, lfu)) {   # rescue (undetermined reads)
      fq1.m <- ifelse(is.element(as.vector(subseq(sread(fq1), 1, 12)), x$V1[x$V3]), 12, 0) + ifelse(is.element(as.vector(subseq(sread(fq1), 1, 11)), x$V1[x$V3]), 11, 0)
    } else {  # regular read files
      fq1.m <- ifelse(is.element(as.vector(subseq(sread(fq1), 1, 12)), x$V1), 12, 0) + ifelse(is.element(as.vector(subseq(sread(fq1), 1, 11)), x$V1), 11, 0)
    }
    n[[f1]] <- n[[f1]] + c(length(fq1.m), sum(fq1.m != 0))

    # filter unmatched
    fq1 <- fq1[fq1.m>0]
    fq2 <- fq2[fq1.m>0]
    fq1.m <- fq1.m[fq1.m>0]

    # extract cell barcode and umi
    fq1.cell <- paste0(as.vector(subseq(sread(fq1), 1, fq1.m)), ifelse(fq1.m==11, "N", ""))
    fq1.umi <- as.vector(subseq(sread(fq1), fq1.m+1, fq1.m+8))

    # shorten R2 and add lib barcode, cell barcode, and umi to read bname
    if(!grepl("PCR", s)) fq2 <- narrow(fq2, 1, 50)
    fq2@id <- BStringSet(paste0(sub(" .:N:0:", "_", as.vector(fq2@id)), "_", fq1.cell, "_", fq1.umi))

    # write combined fastq file
    #for(i in sort(unique(fq1.cell))) {
      # message(i)
      #ii <- fq1.cell == i

      #paste0(sub(".txt$", "", l), ".", s, ".good.txt")

      writeFastq(fq2, file = paste0(sub(".txt$", "", l), ".", s, ".fastq.gz"), mode = "a")
      #write.table(fq1.umi[ii], file = paste0(paste0(sub(".txt$", "", l), ".", s, "/"), i, ".umi"), append = TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
    #}

    invisible(gc())
  }

  close(strm1)
  close(strm2)
  message(" done")
  invisible(gc())
}

nn <- do.call(rbind, n)
nn <- rbind(nn, colSums(nn))
rownames(nn)[nrow(nn)] <- "total"
colnames(nn) <- c("all", s)
write.table(nn, file = paste0(sub(".txt$", "", l), ".", s, ".stats.txt"), sep="\t", quote=FALSE)

invisible(gc())

sessionInfo()


