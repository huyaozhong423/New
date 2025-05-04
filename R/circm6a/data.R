setwd('/home/share/bowen/circm6a/')
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


######################## 筛选过滤human的protein coding CircRNAs (P data)
######################## data from riboCIRC database: http://www.ribocirc.com  
data <- read.csv('./EVTCR.csv',header = T)
table(data$Genome.assembly)
lst <- strsplit(data$Genome.position,':')
remain <- unlist(lapply(lst,function(x) x[1]))
data$seq <- remain
remain2 <- unlist(lapply(lst,function(x) x[2]))
lst <- strsplit(remain2,'-')
data$start <- unlist(lapply(lst,function(x) x[1]))
data$end <- unlist(lapply(lst,function(x) x[2]))
data[24,]

indx <- which(data$Genome.assembly == 'hg38' | data$Genome.assembly == 'hg19')
data <- data[indx,]
data$id <- paste0('CircRNA_',1:nrow(data))
indx <- which(data$Genome.assembly == 'hg38')
lst <- strsplit(data$Strand,"'")
data$Strand <- unlist(lapply(lst,function(x) x[2]))
saveRDS(data,'./EVTCR_human.rds')

gr_hg38 <- GRanges(seqnames = data$seq[indx],
                   IRanges(start = as.numeric(data$start[indx]),
                           end = as.numeric(data$end[indx])),
                   strand = data$Strand[indx],
                   id = data$id[indx])

gr_hg19 <- GRanges(seqnames = data$seq[-indx],
                   IRanges(start = as.numeric(data$start[-indx]),
                           end = as.numeric(data$end[-indx])),
                   strand = data$Strand[-indx],
                   id = data$id[-indx])

chain <- import.chain("/home/share/bowen/circm6a/circrna_LL/hg19ToHg38.over.chain")
circ19To38 <- unlist(liftOver(gr_hg19,chain))
gr_hg38 <- c(gr_hg38,circ19To38)
saveRDS(gr_hg38,'./circ_hg38.rds')


######################## data from TransCirc database
data <- read.table('/home/share/bowen/circm6a/transCirc_144.tsv',header = T, sep = '\t')
table(data$species)
data$id <- paste0('TransCirc_',1:nrow(data))

gr <- GRanges(seqnames = data$chrom,
              IRanges(start = as.numeric(data$start),
                      end = as.numeric(data$end)),
              strand = data$strand,
              id = data$id)

gr$sequence <- DNAStringSet(Views(Hsapiens,gr))

## combind
gr$sequence <- as.character(gr$sequence)
gr <- c(circ_hg38,gr)
gr$label <- 'positive'
saveRDS(gr,'circrna_positive.rds')


######################## 筛选过滤human的 non protein coding CircRNAs (N data)
######################## data from Circ-Atlas 3.0 database  
data <- read.table('./humanCirc_all_bed_v3.0.txt',sep = '\t',header = T)
#data[900,]

gr <- GRanges(seqnames = data$Chro,
              IRanges(start = as.numeric(data$Start),
                      end = as.numeric(data$End)),
              strand = data$Strand,
              id = data$circAltas_ID)
gr <- keepStandardChromosomes(gr,pruning.mode="coarse")
seqlevels(gr)

circrna_positive <- readRDS("~/circm6a/circrna_positive.rds")
olp <- findOverlaps(gr,circrna_positive,ignore.strand = T)
qhits <- queryHits(olp)
gr <- gr[-unique(qhits)]
gr_n_all <- GRanges()

for (j in 1:10) {
  gr_n <- GRanges()
  
  for (i in seq_along(circrna_positive)) {
    width_pos <- width(circrna_positive[i])
    
    indx <- which(width(gr) == width_pos)
    if(length(indx) != 0){
      indx3 <- sample(indx,1)
    }else{
      indx2 <- which(width(gr) >= (width_pos - 50) & width(gr) <= (width_pos + 50))
      indx3 <- sample(indx2,1)
    }
    
    gr_n <- c(gr_n,gr[indx3])
    gr <- gr[-indx3]
    rm(width_pos,indx,indx2,indx3)
  }
  
  gr_n$sequence <- as.character(DNAStringSet(Views(Hsapiens,gr_n)))
  gr_n$label <- paste0('negative_set_',j)
  gr_n_all <- c(gr_n_all,gr_n)
  
  rm(gr_n)
  print(j)
}

length(unique(gr_n_all))
saveRDS(gr_n_all,'circrna_negative.rds')







