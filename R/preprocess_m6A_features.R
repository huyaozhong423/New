library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(dplyr)
library(reticulate)
np <- import("numpy")
# Clean metadata
metadata <- mcols(m6a_circ_olp_neg_updated)
metadata[] <- lapply(metadata, function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.character(x)) return(ifelse(x == "" | is.na(x), "0", x))
  if (is.numeric(x)) return(ifelse(is.na(x), 0, x))
  return(x)
})

# Min-max normalization
normalize_min_max <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
metadata[c("RBP_Num", "SplicingSite_Num", "ORF_Count")] <- 
  lapply(metadata[c("RBP_Num", "SplicingSite_Num", "ORF_Count")], normalize_min_max)

# Compute distances
dist_start <- as.numeric(metadata$Dist_to_start)
dist_end <- as.numeric(metadata$Dist_to_end)
metadata$Dist_to_ORF <- ifelse(dist_start + dist_end == 0, 0, dist_start / (dist_start + dist_end))

dist_ires_start <- as.numeric(metadata$Dist_to_IRES_start)
dist_ires_end <- as.numeric(metadata$Dist_to_IRES_end)
metadata$Dist_to_IRES <- ifelse(dist_ires_start + dist_ires_end == 0, 0, dist_ires_start / (dist_ires_start + dist_ires_end))

mcols(m6a_circ_olp_neg_updated) <- metadata
# Convert GRanges to data.frame
pos_df <- as.data.frame(pos)
neg_df <- as.data.frame(neg)

# Sample 1103 negative samples
set.seed(42)
neg_sample <- neg_df[sample(1:nrow(neg_df), 1103), ]

# Remove genomic position columns
neg_df <- neg_df %>% select(-seqnames, -start, -end, -width, -strand)

# Expand sequence region by ±20 bp and extract genome sequence
pos_expanded <- resize(pos, width = width(ranges(pos)) + 40, fix = "center")
pos$seq_count <- as.character(getSeq(Hsapiens, pos_expanded))

# Base pairing score encoding
paired <- function(x, y, lambda = 0.8) {
  if ((x == 5 && y == 6) || (x == 6 && y == 5)) return(2)
  else if ((x == 4 && y == 7) || (x == 7 && y == 4)) return(3)
  else if ((x == 4 && y == 6) || (x == 6 && y == 4)) return(lambda)
  else return(0)
}

encode_sequence <- function(seq) {
  sapply(strsplit(seq, "")[[1]], function(base) {
    if (base %in% c("A", "a")) return(5)
    else if (base %in% c("C", "c")) return(7)
    else if (base %in% c("G", "g")) return(4)
    else if (base %in% c("T", "U", "t", "u")) return(6)
    else return(3)
  })
}

gaussian <- function(x) exp(-0.5 * x^2)

create_score_matrix <- function(encoded_seq, base_range = 30, lambda = 0.8) {
  len <- length(encoded_seq)
  score_matrix <- matrix(0, nrow = len, ncol = len)
  for (offset in 1:base_range) {
    for (i in 1:(len - offset)) {
      j <- i + offset
      score <- paired(encoded_seq[i], encoded_seq[j], lambda)
      if (score != 0) {
        g <- gaussian(offset)
        score_matrix[i, j] <- score * g
        score_matrix[j, i] <- score * g
      }
    }
  }
  return(score_matrix)
}

# Load samples
neg_all <- readRDS("neg_new1.rds")
pos <- readRDS("pos_new1.rds")
pos_df <- read.csv("pos_df1.csv")

# Binarize ORF
neg_all$ORF <- ifelse(neg_all$ORF != 0, 1, 0)

# Generate 10 datasets with sampled negatives
for (i in 1:10) {
  set.seed(100 + i)
  neg_sampled <- neg_all[sample(1:length(neg_all), 1103)]
  neg_df <- as.data.frame(mcols(neg_sampled))
  neg_df$seq_count <- NULL
  neg_df$label <- 0
  
  pos_df1 <- pos_df[, !grepl("^seq_count\\.V", colnames(pos_df)) | colnames(pos_df) == "label"]
  pos_df1$label <- 1
  
  combined_df <- rbind(pos_df1, neg_df)
  phenotype_features <- combined_df[, setdiff(colnames(combined_df), "label")]
  phenotype_array <- as.matrix(phenotype_features)
  labels_array <- as.numeric(combined_df$label)
  
  pos_matrices <- lapply(pos$seq_count, function(seq) {
    create_score_matrix(encode_sequence(seq), base_range = 30)
  })
  neg_matrices <- lapply(neg_sampled$seq_count, function(seq) {
    create_score_matrix(encode_sequence(seq), base_range = 30)
  })
  all_matrices <- c(pos_matrices, neg_matrices)
  score_array <- array(unlist(all_matrices), dim = c(41, 41, length(all_matrices)))
  score_array <- aperm(score_array, c(3, 1, 2))  # [n, 41, 41]
  
  np$save(paste0("score_matrices_", i, ".npy"), score_array)
  np$save(paste0("phenotype_features_", i, ".npy"), phenotype_array)
  np$save(paste0("labels_", i, ".npy"), labels_array)
}


