# Script to quantify genes after reads have been aligned using STAR.
# Volker Hovestadt (hovestadt@broadinstitute.org), Jan 14, 2019
#
# From the methods:
#
# Sequencing reads were aligned to the human genome (hg38) using STAR version 2.5.2b and default
# parameters. Alignments were guided by using RefSeq gene annotations. Transcripts were quantified
# using the “--quantMode TranscriptomeSAM” option. This resulted in two alignment files, one in which
# reads were aligned to the genome, and one which contained pseudo-alignments to the transcriptome.
# 
# The transcriptome alignments were used to quantify gene expression. For every read all the unique
# gene names of the transcripts the read aligned to were recorded. Some reads aligned to multiple
# genes, which often reflected a primary gene and one or more pseudo-, antisense-, or
# readthrough-genes. We checked if any of the gene names was contained in all the other gene names,
# with a “-” before or after (antisense- and readthrough-genes), or followed by “P” and a digit
# (pseudo-genes). If this was the case, we only kept the primary gene. Reads that still mapped to
# multiple genes were filtered. In a second step, all reads that mapped to the same gene and had an
# identical UMI sequence were collapsed. This yielded a digital expression matrix consisting of UMI
# counts for each cell and gene.
# 
# For all downstream analysis we required cells to have at least 1,000 UMIs (gene counts, indicative
# of the number of captured transcripts) mapping to at least 500 unique genes. We additionally
# excluded cells for which more than 20% of the gene counts reflected either mitochondrial genes or
# ribosomal RNAs, as these likely reflected poor quality cells. Resulting digital expression matrices
# for each sample were deposited in GEO. For downstream analyses, we normalized gene counts to a total
# of 10,000 for each cell.


options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)


message("\n### ", Sys.time(), " - start processing")

suppressMessages(library(data.table))
suppressMessages(library(Rsamtools))
suppressMessages(library(stringi))

cutf <- function(x, f=1, d="/", ...) sapply(strsplit(x, d), function(i) i[f], ...)


### gene annotations
message("\n### ", Sys.time(), " - loading gene annotations")
a <- read.table("hg38.refMrna.gene.chrM.ercc.fa.names", sep = "\t", fill = NA)
a$V2 <- cutf(a$V2, 1, "\\.")

# mitochondrial genes
a <- rbind(a, read.table("hg38.chrM.fa.names", sep = "\t"))
# table(table(a$V2))

# fix ribosomal genes
a$V1[grepl("^MTRNR2", a$V1)] <- "MTRNR2"    # mitochondrial 16S
a$V1[grepl("^RNA5S", a$V1)] <- "RNA5S"
a$V1[grepl("^RNA5-8S", a$V1) | grepl("^RNA5-8N2", a$V1)] <- "RNA5-8S"
a$V1[grepl("^RNA18S", a$V1)] <- "RNA18S"
a$V1[grepl("^RNA28S", a$V1)] <- "RNA28S"
a$V1[grepl("^RNA45S", a$V1)] <- NA          # 28S / 5.8S / 18S transcription unit
a$V1[grepl("LOC109910383", a$V1)] <- NA     # 18S pseudogene
a$V1[grepl("LOC109910384", a$V1)] <- NA     # 5.8S pseudogene

# other manual fixes
a$V1[grepl("TALAM1", a$V1)] <- NA           # antisense of MALAT1
message(paste0(length(unique(a$V2[!is.na(a$V2)])), " transcripts"))
message(paste0(length(unique(a$V1[!is.na(a$V1)])), " genes"))


### read data
message("\n### ", Sys.time(), " - reading alignments to transcriptome")

b <- commandArgs(trailingOnly=TRUE)[1]
#b <- "3_fastq/170728.BM10.star/alignment_3562files_Aligned.toTranscriptome.out.bam"
# file.exists(b)

# read stats
d.stats <- read.table(paste(c(head(strsplit(b, "/")[[1]], -1), "readcount.txt"), collapse="/"))  # all reads
rownames(d.stats) <- cutf(d.stats$V1, nchar(d.stats[1, 1])-nchar(gsub("/", "", d.stats[1, 1]))+1)
d.stats$V1 <- NULL
colnames(d.stats) <- "reads"

foo <- read.table(sub("toTranscriptome.out.bam", "sortedByCoord.out.bam.aligned", b))  # aligned to genome
d.stats$genome <- NA
d.stats[foo$V1, "genome"] <- foo$V2
rm(foo)

# bam file
bamFile <- BamFile(b)

yieldSize(bamFile) <- 1E6
open(bamFile)

x.rbind <- list()

message("processing ", appendLF=FALSE)
repeat {
  message("*", appendLF = FALSE)
  x <- scanBam(bamFile, param = ScanBamParam(what = c("rname", "qname")))
  if(length(x[[1]][[1]]) == 0) break

  # match to annotations
  ax <- match(sub("_dup.", "", levels(x[[1]]$rname)), a$V2)  # what is the _dup. supposed to mean?
  # levels(x[[1]]$rname)[which(is.na(ax))]
  x.gene <- a$V1[ax[as.numeric(x[[1]]$rname)]]

  # merge
  x <- as.data.table(stri_split_fixed(x[[1]]$qname, "_", simplify = TRUE))
  x$V2 <- NULL   # lib bc not needed
  x$V3 <- factor(x$V3)
  x$gene <- x.gene   # dont use factors for gene, because split below
  rm(x.gene)

  # unique
  x <- x[!is.na(x$gene)]  # some genes are removed
  x <- unique(x)

  x.rbind <- c(x.rbind, list(x))

  invisible(gc())
}

# rbind, unique
x <- rbindlist(x.rbind)
rm(x.rbind)
invisible(gc())

x <- unique(x)
invisible(gc())

close(bamFile)
message(" done")


### if multiple names, check if one contained in all others (e.g. readthrough genes, antisense, pseudogenes, ...)
message("\n### ", Sys.time(), " - filtering readthrough, pseudogenes, etc")

x[, N:=.N, by=V1]  # count how many genes a read is aligned to

# stats
foo <- x[, length(unique(V1)), by=V3]  # aligned to transcriptome
d.stats$transcriptome <- 0
d.stats[as.vector(foo$V3), "transcriptome"] <- foo$V1

foo <- x[N==1, length(unique(V1)), by=V3]  # uniquely aligned to transcriptome (before filtering)
d.stats$unique.prefilter <- 0
d.stats[as.vector(foo$V3), "unique.prefilter"] <- foo$V1


# split
xx <- split(x[!is.na(gene) & N>1, gene], x[!is.na(gene) & N>1, V1])

# process only unique gene combinations
xx.unique <- unique(xx)
xx.filter <- lapply(xx.unique, function(xxx) {
  xxxm <- lapply(xxx, regexpr, xxx)  # match each to all
  xxxw <- sapply(lapply(xxxm, `!=`, -1), all)  # maximum 1 name can be contained in all others
  if(any(xxxw)) {
    xxxw0 <- substr(xxx[!xxxw], xxxm[[which(xxxw)]][!xxxw]-1, xxxm[[which(xxxw)]][!xxxw]-1)  # character before
    xxxw1 <- substr(xxx[!xxxw], xxxm[[which(xxxw)]][!xxxw]+nchar(xxx[xxxw]), xxxm[[which(xxxw)]][!xxxw]+nchar(xxx[xxxw]))  # character after
    xxxw2 <- substr(xxx[!xxxw], xxxm[[which(xxxw)]][!xxxw]+nchar(xxx[xxxw])+1, xxxm[[which(xxxw)]][!xxxw]+nchar(xxx[xxxw])+1)  # character two after
    if(all(xxxw0 == "-" | xxxw1 == "-" | (xxxw1 == "P" & is.element(xxxw2, sapply(seq(9), toString))))) {  # either dash before, after, or followed by PX (pseudogene)
      return(xxx[xxxw])
    } else {
      return(xxx)
    }
  } else {
    return(xxx)
  }
})
invisible(gc())

# only continue with the ones that are now unique
xx.unique <- xx.unique[lengths(xx.filter) == 1]
xx.filter <- xx.filter[lengths(xx.filter) == 1]

# match split to unique
xx.match <- match(xx, xx.unique)

# put back together in x
x <- rbind(x[N==1],
           merge(data.table(V1=names(xx)[!is.na(xx.match)], gene=unlist(xx.filter)[xx.match[!is.na(xx.match)]]), x[N>1], by=c("V1", "gene")))  # unique alignments + unique alignments after filtering - remaining

invisible(gc())

rm(list=ls()[grepl("^xx", ls())])
invisible(gc())

# merge
x$N <- NULL  # not needed anymore
x$gene <- factor(x$gene, sort(unique(a$V1)))  # now ok to use factor

foo <- x[, length(unique(V1)), by=V3]  # uniquely aligned to transcriptome (after gene filter)
d.stats$unique <- 0
d.stats[as.vector(foo$V3), "unique"] <- foo$V1


###  digital expression matrix
message("\n### ", Sys.time(), " - generate digital expression matrix")
x$V1 <- NULL   # remove read names and make unique
x <- unique(x)   # much smaller, one row per umi
invisible(gc())

d <- as.data.frame(dcast(x[, .N, by=c("V3", "gene")], gene~V3, fill = 0, value.var = "N", drop = FALSE))
rownames(d) <- d$gene
d$gene <- NULL

gz <- gzfile(paste0(paste(head(strsplit(b, "/")[[1]], -1), collapse="/"), ".expr.txt.gz"), "w")
write.table(d, file = gz, sep="\t", quote=FALSE)
close(gz)

# stats
d.stats$umi <- 0
d.stats$umi.ribo <- 0
d.stats$umi.mito <- 0
d.stats$umi.ercc <- 0

d.stats[colnames(d), "umi"] <- colSums(d)
d.stats[colnames(d), "umi.ribo"] <- colSums(d[is.element(rownames(d), c("RNA5S", "RNA18S", "RNA28S")), ])
d.stats[colnames(d), "umi.mito"] <- colSums(d[is.element(rownames(d), c("MTRNR2")) | grepl("^MT-", rownames(d)), ])
d.stats[colnames(d), "umi.ercc"] <- colSums(d[grepl("^ERCC-", rownames(d)), ])

save(list=c("d", "d.stats"), file = paste0(paste(head(strsplit(b, "/")[[1]], -1), collapse="/"), ".expr.RData"))
write.table(d.stats, paste0(paste(head(strsplit(b, "/")[[1]], -1), collapse="/"), ".stats.txt"), sep="\t", quote=FALSE)

message("\n### ", Sys.time(), " - all done\n")


