# 1. 定义感兴趣的 RBP 名称
rbp_list <- c('HNRNPC', 'YTHDC1', 'YTHDF1', 'YTHDF2', 'METTL3', 'METTL14', 'METTL16', 'WTAP', 'ALKBH5', 'FTO')
RBP_site_gr_38 <- readRDS("RBP_site_gr_38.rds")
# 2. 初始化新的列，默认值为 FALSE
# 假设 m6a_circ_olp_updated 是一个 GRanges 对象
for (rbp in rbp_list) {
  mcols(m6a_circ_olp_neg_updated)[[rbp]] <- FALSE
}

# 3. 遍历 RBP 列表，计算是否有重叠
for (rbp in rbp_list) {
  # 过滤出当前 RBP 的位点
  rbp_sites <- RBP_site_gr_38[mcols(RBP_site_gr_38)$RBP_name == rbp]
  
  # 找到 m6A 位点和当前 RBP 位点的重叠
  overlaps <- findOverlaps(m6a_circ_olp_updated, rbp_sites)
  
  # 标记重叠的 m6A 位点
  m6a_circ_olp_updated[[rbp]][queryHits(overlaps)] <- TRUE
}

# 4. 查看更新后的数据框
head(m6a_circ_olp_updated)

for (rbp in rbp-list){
  rbp_sites <- RBP_site_gr_38[mcols(RBP_site_gr_38)$RBP_name==rbp]
  overlaps<-findOverlaps(m6a_circ_olp_updated,rbp_sites)
  m6a_circ_olp_updated[[rbp]][queryHits(overlaps)]<-TRUE
}

# 假设 rbp_list 是 RBP 名称的列表
rbp_list <- c('HNRNPC', 'YTHDC1', 'YTHDF1', 'YTHDF2', 'METTL3', 'METTL14', 'METTL16', 'WTAP', 'ALKBH5', 'FTO') # 替换为实际的 RBP 名称

# 遍历每个 RBP
for (rbp in rbp_list) {
  # 过滤出当前 RBP 的位点
  rbp_sites <- RBP_site_gr_38[mcols(RBP_site_gr_38)$RBP_name == rbp]
  
  # 找到 m6A 位点和当前 RBP 位点的重叠
  overlaps <- findOverlaps(m6a_circ_olp_neg_updated, rbp_sites)
  
  # 初始化当前 RBP 的元数据列（如果尚未存在）
  if (!rbp %in% colnames(mcols(m6a_circ_olp_neg_updated))) {
    mcols(m6a_circ_olp_neg_updated)[[rbp]] <- FALSE
  }
  
  # 标记重叠的 m6A 位点
  mcols(m6a_circ_olp_neg_updated)[[rbp]][queryHits(overlaps)] <- TRUE
}


