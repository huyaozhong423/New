library(package, help, pos = 2, lib.loc = NULL)
# 将 GRanges 转换为 data.frame
pos_df <- as.data.frame(pos)
neg_df <- as.data.frame(neg)
neg<-readRDS("neg_sampleall.rds")

# 1. 随机抽取 1103 个 neg 样本
set.seed(42)  # 保证可复现
neg_sample <- neg_df[sample(1:nrow(neg_df), 1103), ]

# 2. 合并 pos 和 neg
df <- rbind(pos_df, neg_sample)
library(dplyr)

neg_df <- neg_df %>%
  select(-seqnames, -start, -end, -width, -strand)




library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

# 让 data 的前后各扩展 20bp
data_expanded <- resize(data, width = width(data) + 40, fix = "center")

# 直接使用 Views() 提取序列
seqs <- DNAStringSet(Views(Hsapiens, data_expanded))

# 添加到原 GRanges
mcols(data)$seq_expanded <- as.character(seqs)

# 查看结果
head(data)

neg_new1<-readRDS("neg_new1")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# 1. 确保 Bioconductor 及 BSgenome 安装
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
}
library(BSgenome.Hsapiens.UCSC.hg38)

# 2. 提取 neg 对应的 ranges
pos_ranges <- ranges(pos)  # 只获取 GRanges 的 ranges 信息

# 3. 调整范围，前后各扩展 20bp
pos_expanded <- resize(pos, width = width(pos_ranges) + 40, fix = "center")

# 4. 获取基因组序列
pos_sequences <- getSeq(Hsapiens, pos_expanded)

# 5. 转换为 DNAStringSet
pos_dna_string_set <- DNAStringSet(pos_sequences)

# 6. 查看结果
pos_dna_string_set
# 确保 neg_dna_string_set 是 DNAStringSet 并转换为 character 类型
pos$seq_count <- as.character(pos_dna_string_set)
head(pos)


one_hot_encode <- function(sequence) {
  bases <- c("A", "C", "G", "T")
  seq_chars <- unlist(strsplit(sequence, split = ""))
  one_hot <- t(sapply(seq_chars, function(base) as.numeric(bases == base)))
  return(one_hot)
}

# 对 `seq_count` 中的所有序列进行 one-hot 编码
one_hot_list <- lapply(pos$seq_count, one_hot_encode)


# 将 one-hot 编码展平，每个样本变成一个长向量
flatten_one_hot <- function(one_hot) {
  return(as.vector(t(one_hot)))  # 转换为行向量
}
flatten_one_hot <- function(one_hot) {
  return(as.vector(t(one_hot)))  # 转置后展平，变为行向量
}

# 对所有样本进行转换
svm_input <- do.call(rbind, lapply(one_hot_list, flatten_one_hot))
# 将 svm_input 转换为数据框
svm_input_df <- as.data.frame(svm_input)
dim(svm_input_df)
# 替换 neg 数据框中的 seq_count 列
pos$seq_count <- svm_input_df
pos$seq_count
head(pos)
saveRDS(neg_sample,file="neg_.rds")
saveRDS(pos,file="neg_sampleall.rds")


neg<-readRDS("neg_new1.rds")
neg
