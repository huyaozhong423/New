################################################################################
####                                   circRNA                           ####
################################################################################
setwd('/home/share/bowen/circrna/')
library(GenomicRanges)
library(dplyr)
library(permute)
library(readr)
library(rtracklayer)
library(BSgenome)
library(GenomicFeatures)
#dir.create('/home/share/bowen/circrna/')

############################ circRNA
data <- read.csv('/home/share/bowen/circrna/circRNA_206.csv')
data[24,]
data$id <- paste0('circRNA_',1:nrow(data))
table(data$Genome.assembly)
indx <- which(data$Genome.assembly == 'hg38')
data_38 <- data[indx,]
data_19 <- data[-indx,]

circ38 <- GRanges(seqnames = data_38$seqnames,
                IRanges(start = as.numeric(data_38$start),
                        end = as.numeric(data_38$end)),
                strand = data_38$strand,
                id = data_38$id)
  
circ19 <- GRanges(seqnames = data_19$seqnames,
                  IRanges(start = as.numeric(data_19$start),
                          end = as.numeric(data_19$end)),
                  strand = data_19$strand,
                  id = data_19$id) 

chain <- import.chain("/home/share/bowen/circrna/hg19ToHg38.over.chain")
circ19To38 <- unlist(liftOver(circ19,chain))
indx <- which(duplicated(circ19To38$id))
circ19To38 <- circ19To38[-indx]
circ38 <- c(circ38,circ19To38)
rm(circ19To38,indx,data_19,data_38,data,circ19,chain)
saveRDS(circ38,'circRNA_hg38.rds')

##################################### m6A 
data <- readRDS('/home/share/bowen/circm6a/hg38_CL_Tech.rds')
olp <- findOverlaps(circ38,data)
length(unique(queryHits(olp)))
length(unique(subjectHits(olp)))
q <- queryHits(olp)
s <- subjectHits(olp)
m6a_circ_olp <- data[s]
m6a_circ_olp$circRNA <- circ38$id[q]
saveRDS(m6a_circ_olp,'m6a_circ_olp.rds')

## negative
negative_circRNA_hg38 <- readRDS("~/circrna/negative_circRNA_hg38.rds")
circ <- unlist(negative_circRNA_hg38)
olp <- findOverlaps(circ,data)
length(unique(queryHits(olp)))
length(unique(subjectHits(olp)))
q <- queryHits(olp)
s <- subjectHits(olp)
m6a_circ_olp <- data[s]
m6a_circ_olp$circRNA <- circ$circID[q]
saveRDS(m6a_circ_olp,'m6a_circ_olp_neg.rds')








  