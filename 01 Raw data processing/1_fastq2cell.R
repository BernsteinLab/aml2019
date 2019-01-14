# Script to perform the first step of the raw data processing.
# Volker Hovestadt (hovestadt@broadinstitute.org), Jan 14, 2019
#
# From the methods:
#
# To analyze our single-cell sequencing data, we employed an approach to annotate sequencing reads by
# CB before sequence alignment and quantification. First, we counted all unique 12 bp CBs for each
# library. We excluded CBs occurring less than 100 times, and filtered barcodes containing stretches
# of eight identical nucleotides. Next, we excluded CBs that were associated with non-random UMIs. For
# all reads associated with a given CB, we checked that the frequency of any nucleotide did not exceed
# 90% at each base of the UMI. The majority of reads filtered this way contained part of the Tn5
# binding sequence, i.e. reflected events in which the transposase integrated within the CB/UMI,
# yielding very short fragments with invariable (non-random) UMI sequences.
# 
# We noticed that a number of CBs (5-20%, depending on the batch of barcoded beads) were associated
# with UMIs that contained a Thymine as the last nucleotide. These sequences often represent CBs in
# which a single nucleotide is missing due to errors in the split-pool synthesis. In this case reads
# start with a 11 bp CB, followed by the 8 bp UMI and the first base of the poly-T sequence that
# hybridizes to the poly-A tail of captured mRNAs. If not corrected, this causes a single cell to
# produce four different single-cell transcriptomes. We corrected these barcodes if in fact four
# different CBs were detected with a similar number of total reads that were variable in their last
# base. The UMIs were also corrected accordingly.
# 
# To filter out CBs that likely resulted from sequencing errors, we ranked all CBs according to their
# number of reads (requiring at least 1,000 reads). We filtered out all CBs that had a higher ranked
# CB that was different in only one position (hamming distance 1).
# 
# This script is run for every sample in each run. A single sample can have multiple library barcodes
# in the sample sheet. It generates two text files containing either all or the filtered list of cell
# barcodes.


options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressMessages(library(data.table))
suppressMessages(library(ShortRead))
suppressMessages(library(stringi))


### Functions
cutf <- function(x, f=1, d="/", ...) sapply(strsplit(x, d), function(i) i[f], ...)

### Arguments
d <- commandArgs(trailingOnly=TRUE)[1]  # one or multiple directories containing fastq files, not searched recursively
l <- commandArgs(trailingOnly=TRUE)[2]  # sample sheet
s <- commandArgs(trailingOnly=TRUE)[3]  # sample name, should be in sample sheet, will be used for output files

# example sample sheet (library barcode, name in filename, sample name):
# CTTGCAGA	N700-BC16	AML150921A-1
# TGCGACCT	N700-BC17	AML150921A-1
# CATTCCGA	N700-BC18	AML150921A-2
# GATAACCT	N700-BC19	AML150921A-2

a <- read.table(l)
if(nrow(a) == 0) stop("no entries in ", l)

lf <- unique(sort(unlist(lapply(strsplit(d, ","), list.files, pattern="_R1_", full.names=TRUE, recursive=FALSE))))  # bcl2fastq lib barcode split, R1 only
lf <- lf[grepl(".fastq.gz$", lf)]
if(length(lf) == 0) stop("no fastq.gz files in ", d)

lf <- lf[apply(do.call(cbind, lapply(paste0(sub("N700.", "N700-", a[a$V3 == s, "V2"]), "_S"), grepl, x=lf)), 1, any)]  # subset
if(length(lf) == 0) stop("no sample match in ", d)


### 1 - load R1 fastq file
message("\n### ", Sys.time(), " - reading R1 sequences")
r1 <- list()  # empty
for(f in lf) {
  message("file ", match(f, lf), "/", length(lf), ": ", f, " ", appendLF = FALSE)
  strm <- FastqStreamer(f)  # 1M reads by default
  repeat {
    message("*", appendLF = FALSE)
    fq <- yield(strm)
    if(length(fq) == 0) break
    r1 <- c(r1, list(data.table(r1=substr(sread(fq), 1, 20))[, .N, by="r1"]))
  }
  r1 <- list(rbindlist(r1)[, list(N=sum(N)), by="r1"])
  close(strm)
  message(" done")
  invisible(gc())
}
r1 <- r1[[1]][order(-N, r1)]
invisible(gc())

# process
r1[, r1.1to12:=stri_sub(r1, 1, 12)]
r1.cell <- r1[, list(N=sum(N)), by="r1.1to12"][order(-N)]
message(format(r1[, sum(N)], nsmall=1, big.mark=","), " reads / ", format(r1[, length(unique(r1.1to12))], nsmall=1, big.mark=","), " unique cell barcodes")

# all output plots in one pdf
pdf(paste0(sub(".txt$", "", l), ".", s, ".qc.pdf"), width = 6, height = 6)
par(mar=c(4, 4, 4, 4))


### 2 - filter rare cell bc
message("\n### ", Sys.time(), " - filter cell barcodes that occur less than 100 times")

r1 <- r1[!is.element(r1.1to12, r1.cell[N < 100, r1.1to12])]
message(format(r1[, sum(N)], nsmall=1, big.mark=","), " reads / ", format(r1[, length(unique(r1.1to12))], nsmall=1, big.mark=","), " unique cell barcodes remaining")
invisible(gc())

write.table(r1.cell[N >= 100], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste0(sub(".txt$", "", l), ".", s, ".all.txt"))


### 3 - filter simple cell barcodes
message("\n### ", Sys.time(), " - filter cell barcodes that contain 8bp repeats, or Ns")

r1.cell.simple <- r1.cell[grepl("AAAAAAAA", r1.1to12) | grepl("CCCCCCCC", r1.1to12) | grepl("GGGGGGGG", r1.1to12) | grepl("TTTTTTTT", r1.1to12) | grepl("N", r1.1to12)]

plot(r1.cell[N>=100, N]/1000, log="xy", type="l", xlab="Cell barcodes", ylab="Reads [K]", main="Unfiltered")
axis(3, match(r1.cell.simple$r1.1to12, r1.cell$r1.1to12), labels = FALSE)

r1 <- r1[!is.element(r1.1to12, r1.cell.simple$r1.1to12)]
message(format(r1[, sum(N)], nsmall=1, big.mark=","), " reads / ", format(r1[, length(unique(r1.1to12))], nsmall=1, big.mark=","), " unique cell barcodes remaining")

rm(list=setdiff(ls(), c("r1", "s", "l")))
invisible(gc())


### 4 - umi filter
message("\n### ", Sys.time(), " - filter cell barcodes with broken umi")

message("processing ", appendLF=FALSE)
r1.umi <- lapply(13:20, function(i) {
  message("*", appendLF = FALSE)
  r1[, b:=factor(stri_sub(r1, i, i), c("A", "C", "G", "T", "N"))]
  r1.b <- r1[, list(N=sum(N)), by=c("r1.1to12", "b")]
  r1.b.dcast <- dcast(r1.b, r1.1to12~b, drop = FALSE, value.var = "N", fill=0)
  # dim(r1.b.dcast)
  r1.b.dcast[, sum:=A+C+G+T+N]
  r1.b.dcast[, i:=i]
  r1.b.dcast
})
r1$b <- NULL
message(" done")

barplot(t(as.matrix(rbindlist(r1.umi)[, list(sum(A), sum(C), sum(G), sum(T), sum(N)), by=i])[, 2:6])/1E6, las=2, col=c("green", "blue", "yellow", "red", "grey"),
        names.arg = 13:20, main="UMI base distribution", ylab="Reads [M]")
axis(4, at = seq(0, r1.umi[[1]][, sum(sum)]/1E6, length.out = 5), labels = seq(0, 1, length.out = 5), las=2)

for(i in 13:20) {
  r1.umi.melt <- melt(r1.umi[[i-12]][sum>=500], id.vars=c("r1.1to12", "sum"), measure.vars=c("A", "C", "G", "T"))
  r1.umi.melt <- r1.umi.melt[order(-sum, -value), ]

  plot(r1.umi.melt$sum/1E3, r1.umi.melt$value/r1.umi.melt$sum,
       pch=".", col=c(A="green", C="blue", G="yellow", T="red")[r1.umi.melt$variable],
       xlab="Reads [K]", ylab="Fraction per base", xlim = c(1, 1E4), ylim=c(0, 1), log="x", main=i)
  abline(h=0.9, lty=2)
}

r1.umi.filter <- sort(rbindlist(r1.umi)[A/sum>0.9 | C/sum>0.9 | G/sum>0.9 | (T/sum>0.9 & i <= 19)][, unique(r1.1to12)])  # unexpected barcode, except T in pos 20
r1.umi.syncandidate <- sort(rbindlist(r1.umi)[(T/sum>0.9 & i == 20)][, unique(r1.1to12)])  # candidates for 1bp synth error fix
r1.umi.syncandidate <- r1.umi.syncandidate[!is.element(stri_sub(r1.umi.syncandidate, 1, 11), stri_sub(r1.umi.filter, 1, 11))]  # cant be in filter list

# tn5 adapter
#"CCTGTCTCTTATACACATCT"   5
# "CTGTCTCTTATACACATCTC"  1
#  "TGTCTCTTATACACATCTCC" 3

# ?
# "TGAGCGAAGCAGTGGTATCA"  2
# "      AAGCAGTGGTATCAACGCAG" 4
# message(format(r1[r1.1to12 %in% r1.umi.filter][grepl("TGTCTCTTAT", r1), sum(N)], nsmall=1, big.mark=","), " contain TGTCTCTTAT (Tn5 adapter sequence)")

r1 <- r1[!is.element(r1.1to12, r1.umi.filter)]
invisible(gc())

message(format(r1[, sum(N)], nsmall=1, big.mark=","), " reads / ", format(r1[, length(unique(r1.1to12))], nsmall=1, big.mark=","), " unique cell barcodes remaining")
message(format(r1[r1.1to12 %in% r1.umi.syncandidate, sum(N)], nsmall=1, big.mark=","), " reads are candidates for 1bp barcode synthesis error fix (T at end of UMI)")


### 4 - fix synthesis error
message("\n### ", Sys.time(), " - fix cell barcode synthesis errors")

r1.umi.syncandidate.all <- data.table(r1.1to11=rep(sort(unique(substr(r1.umi.syncandidate, 1, 11))), each=4), r1.12=c("A", "C", "G", "T"))
r1.umi.syncandidate.all[, r1.1to12:=paste0(r1.1to11, r1.12)]
r1.umi.syncandidate.all <- merge(r1.umi.syncandidate.all, r1[, list(N=sum(N)), by="r1.1to12"], by="r1.1to12", all.x = TRUE)
r1.umi.syncandidate.all[is.na(N), N:=0]

r1.umi.syncandidate.all.meansd <- r1.umi.syncandidate.all[, list(sum=sum(N), sd=sd(N)), by="r1.1to11"]
r1.umi.syn <- r1.umi.syncandidate.all[r1.1to11 %in% r1.umi.syncandidate.all.meansd[sum>=1000 & sd/sum*2<=0.25, sort(unique(r1.1to11))], sort(unique(r1.1to12))]
r1.umi.syn.filter <- setdiff(r1.umi.syncandidate, r1.umi.syn)

plot(r1.umi.syncandidate.all.meansd$sum /1E3, r1.umi.syncandidate.all.meansd$sd/r1.umi.syncandidate.all.meansd$sum*2,
     col=ifelse(r1.umi.syncandidate.all.meansd$r1.1to11 %in% stri_sub(r1.umi.syn, 1, 11), "black", "grey"),
     log="x", pch=".", xlim=c(1, 1E4), ylim=c(0, 1), main="Barcode synthesis error candidates", xlab="Total reads [K]", ylab="SD / Sum * 2")
abline(h=0.25)
abline(v=1)

message(format(r1[r1.1to12 %in% r1.umi.syn.filter, sum(N)], nsmall=1, big.mark=","), " candidate reads filtered, remaining reads fixed")

r1 <- r1[!is.element(r1.1to12, c(r1.umi.syn.filter))]
r1[r1.1to12 %in% r1.umi.syn, r1:=paste0(stri_sub(r1, 1, 11), "N", stri_sub(r1, 12, 19))]  # fix r1
# r1$r1.1to12 <- NULL
r1 <- r1[, list(N=sum(N)), by=r1][order(-N, r1)]  # collapse, so column r1 is unique
r1[, r1.1to12:=stri_sub(r1, 1, 12)]
invisible(gc())

message(format(r1[, sum(N)], nsmall=1, big.mark=","), " reads / ", format(r1[, length(unique(r1.1to12))], nsmall=1, big.mark=","), " unique cell barcodes remaining")

rm(list=setdiff(ls(), c("r1", "s", "l")))
invisible(gc())


### 5 - require at least 1000 reads per cell barcode
message("\n### ", Sys.time(), " - only keep cell barcodes with more than 1000 reads")

r1.cell <- r1[, list(N=sum(N)), by="r1.1to12"][N>=1000][order(-N, r1.1to12)]
r1 <- r1[is.element(r1.1to12, r1.cell[, r1.1to12])]

plot(r1.cell[, N]/1E3, log="xy", type="l", xlab="Cell barcodes", ylab="Reads [K]", las=2, main="Filtered")

message(format(r1[, sum(N)], nsmall=1, big.mark=","), " reads / ", format(r1[, length(unique(r1.1to12))], nsmall=1, big.mark=","), " unique cell barcodes remaining")


### 6 - allow 1mm in cell barcodes
message("\n### ", Sys.time(), " - filter barcodes with distance one to a higher ranked barcode")

# all sequences with 1 mismatch
message("processing ", appendLF=FALSE)
r1.cell.1mm <- rbindlist(lapply(seq(1, 12), function(i) {  #
  message("*", appendLF = FALSE)
  rbindlist(lapply(c("A", "C", "G", "T"), function(b) {
    x <- r1.cell[N>=1000]
    x$r1.1to12.mm1 <- x$r1.1to12
    stri_sub(x$r1.1to12.mm1, i, i) <- b
    x
  }))
}))
r1.cell.1mm <- r1.cell.1mm[r1.1to12 != r1.1to12.mm1]  # replaced by the same base (1/4 of rows)
message(" done")

# remove barcodes that are 1mm to a higher ranked barcode
r1.cell.1mm <- merge(r1.cell.1mm, r1.cell, by.x="r1.1to12.mm1", by.y="r1.1to12", all.x = TRUE)
colnames(r1.cell.1mm)[3:4] <- c("N", "N.mm1")
# r1.cell.1mm

r1.cell.1mm.filter <- r1.cell.1mm[N.mm1 >= N, sort(unique(r1.1to12))]
axis(3, at=match(r1.cell.1mm.filter, r1.cell$r1.1to12), labels = FALSE)

r1 <- r1[!is.element(r1.1to12, r1.cell.1mm.filter)]
message(format(r1[, sum(N)], nsmall=1, big.mark=","), " reads / ", format(r1[, length(unique(r1.1to12))], nsmall=1, big.mark=","), " unique cell barcodes remaining")


### 7 - output
message("\n### ", Sys.time(), " - write final list of cell barcodes\n")

r1.cell <- r1[, list(N=sum(N)), by="r1.1to12"][order(-N, r1.1to12)]
write.table(r1.cell, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste0(sub(".txt$", "", l), ".", s, ".good.txt"))

invisible(dev.off())

sessionInfo()

