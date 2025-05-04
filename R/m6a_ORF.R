library(GenomicRanges)
setwd("")
# 使用 findOverlaps 找到 m6a_circ_olp 和 circrna_m6a 之间的重叠关系
overlaps <- findOverlaps(m6a_circ_olp, circrna_data)

# 初始化 m6a_in_ORF 列
m6a_circ_olp$m6a_in_ORF <- 0

# 遍历重叠的 m6a 位点
for (hit in seq_along(overlaps)) {
  query_idx <- queryHits(overlaps)[hit]      # m6a_circ_olp 的索引
  subject_idx <- subjectHits(overlaps)[hit]  # circrna_m6a 的索引
  
  # 获取当前 m6a 位点的起始位置
  m6a_pos <- start(m6a_circ_olp[query_idx])
  
  # 获取当前 circrna_m6a 的 ORF 范围
  orf_range_raw <- as.character(circrna_data$ORF[subject_idx])
  
  # 如果 ORF 为空或者 "0"，跳过
  if (is.na(orf_range_raw) || orf_range_raw == "0") next
  
  # 解析 ORF 范围
  orf_ranges <- strsplit(orf_range_raw, ";")[[1]]
  
  # 遍历所有 ORF，检查 m6a 是否在 ORF 范围内
  for (orf in orf_ranges) {
    orf_start_end <- as.integer(strsplit(orf, "-")[[1]])
    orf_start <- orf_start_end[1]
    orf_end <- orf_start_end[2]
    
    # 判断 m6a 是否在 ORF 范围内
    if (m6a_pos >= orf_start && m6a_pos <= orf_end) {
      m6a_circ_olp$m6a_in_ORF[query_idx] <- 1  # 只要有一个匹配，就标记为 1
      break  # 找到一个匹配后即可跳出循环
    }
  }
}


# 检查结果
table(m6a_circ_olp$m6a_in_ORF)
library(GenomicRanges)

# 使用 findOverlaps 找到 m6a_circ_olp 和 circrna_m6a 之间的重叠关系
overlaps <- findOverlaps(m6a_circ_olp, circrna_data)

# 初始化 ORF 和 ORF_Count 列
m6a_circ_olp$ORF <- NA

m6a_circ_olp$ORF_Count <- 0

# 遍历所有匹配项
for (i in seq_along(overlaps)) {
  if (i %% 100 == 0) {
    cat("正在处理第", i, "行，共", length(m6a_circ_olp), "行\n")
  }
  query_idx <- queryHits(overlaps)[i]      # m6a_circ_olp 的索引
  subject_idx <- subjectHits(overlaps)[i]  # circrna_m6a 的索引
  
  # 获取 circrna_m6a 对应的 ORF 和 ORF_Count
  orf_value <- as.character(circrna_data$ORF[subject_idx])
  orf_count <- circrna_data$ORF_Count[subject_idx]
  
  # 处理 NA 值
  if (is.na(orf_value)) next
  
  # 如果 `m6a_circ_olp` 已经有 ORF 数据，则合并
  if (!is.na(m6a_circ_olp$ORF[query_idx])) {
    m6a_circ_olp$ORF[query_idx] <- paste(m6a_circ_olp$ORF[query_idx], orf_value, sep = ";")
    m6a_circ_olp$ORF_Count[query_idx] <- m6a_circ_olp$ORF_Count[query_idx] + orf_count
  } else {
    m6a_circ_olp$ORF[query_idx] <- orf_value
    m6a_circ_olp$ORF_Count[query_idx] <- orf_count
  }
}

# 检查结果
head(m6a_circ_olp$ORF)
head(m6a_circ_olp$ORF_Count)

# 确保 ORF 列是字符型
m6a_circ_olp$ORF <- as.character(m6a_circ_olp$ORF)

# 遍历 m6a_circ_olp 的每一行
for (i in seq_along(m6a_circ_olp$ORF)) {
  if (!is.na(m6a_circ_olp$ORF[i]) && m6a_circ_olp$ORF[i] != "0") {
    # 拆分 ORF 记录（以 ";" 分割）
    unique_orfs <- unique(strsplit(m6a_circ_olp$ORF[i], ";")[[1]])
    
    # 重新组合 ORF（去重后用 ";" 连接）
    m6a_circ_olp$ORF[i] <- paste(unique_orfs, collapse = ";")
    
    # 更新 ORF_Count（去重后数量）
    m6a_circ_olp$ORF_Count[i] <- length(unique_orfs)
  }
}

# 检查结果
head(m6a_circ_olp$ORF)
head(m6a_circ_olp$ORF_Count)


# 初始化 m6a_in_ORF 列，默认值为 0
m6a_circ_olp$m6a_in_ORF <- 0

# 遍历 m6a_circ_olp 数据的每一行
for (i in seq_along(m6a_circ_olp$ORF)) {
  orf_range_raw <- as.character(m6a_circ_olp$ORF[i])  # 获取 ORF 范围
  m6a_pos <- start(m6a_circ_olp[i])  # 获取 m6A 位点
  
  # 如果 ORF 为空或为 "0"，则跳过
  if (is.na(orf_range_raw) || orf_range_raw == "0") {
    next
  }
  
  # 拆分 ORF 记录（以 ";" 分割多个 ORF）
  orf_ranges <- strsplit(orf_range_raw, ";")[[1]]
  
  # 遍历所有 ORF 范围，检查 m6A 是否落入其中
  for (orf in orf_ranges) {
    orf_start_end <- as.integer(strsplit(orf, "-")[[1]])  # 解析 ORF 起始和终止坐标
    orf_start <- orf_start_end[1]
    orf_end <- orf_start_end[2]
    
    # 判断 m6A 是否在 ORF 范围内
    if (m6a_pos >= orf_start && m6a_pos <= orf_end) {
      m6a_circ_olp$m6a_in_ORF[i] <- 1  # 只要有一个匹配，就设为 1
      break  # 找到一个匹配后即可跳出循环
    }
  }
}

# 查看前几行结果
head(m6a_circ_olp$m6a_in_ORF)




# 初始化 Dist_to_start 和 Dist_to_end 列，默认值为 NA
m6a_circ_olp$Dist_to_start <- NA
m6a_circ_olp$Dist_to_end <- NA

# 遍历 m6a_circ_olp 数据的每一行
for (i in seq_along(m6a_circ_olp$ORF)) {
  if (i %% 100 == 0) {
    cat("正在处理第", i, "行，共", length(m6a_circ_olp), "行\n")
  }
  orf_range_raw <- as.character(m6a_circ_olp$ORF[i])  # 获取 ORF 范围
  m6a_pos <- start(m6a_circ_olp[i])  # 获取 m6A 位点
  
  # 如果 ORF 为空或为 "0"，跳过
  if (is.na(orf_range_raw) || orf_range_raw == "0") {
    next
  }
  
  # 拆分 ORF 记录（以 ";" 分割多个 ORF）
  orf_ranges <- strsplit(orf_range_raw, ";")[[1]]
  
  # 初始化存储距离的向量
  start_distances <- c()
  end_distances <- c()
  
  # 遍历所有 ORF 范围，计算 m6A 到 ORF 两端的距离
  for (orf in orf_ranges) {
    orf_start_end <- as.integer(strsplit(orf, "-")[[1]])  # 解析 ORF 起始和终止坐标
    orf_start <- orf_start_end[1]
    orf_end <- orf_start_end[2]
    
    # 仅当 m6A 在 ORF 内时，计算距离
    if (m6a_pos >= orf_start && m6a_pos <= orf_end) {
      start_distances <- c(start_distances, m6a_pos - orf_start)
      end_distances <- c(end_distances, orf_end - m6a_pos)
    }
  }
  
  # 如果找到 ORF 并计算了距离，则存入数据框
  if (length(start_distances) > 0) {
    m6a_circ_olp$Dist_to_start[i] <- paste(start_distances, collapse = ",")
    m6a_circ_olp$Dist_to_end[i] <- paste(end_distances, collapse = ",")
  }
}

# 查看前几行计算结果
head(m6a_circ_olp[, c("ORF", "m6a_in_ORF", "Dist_to_start", "Dist_to_end")])

m6a_circ_olp$ORF <- ifelse(m6a_circ_olp$ORF != "0" & !is.na(m6a_circ_olp$ORF), 1, 0)