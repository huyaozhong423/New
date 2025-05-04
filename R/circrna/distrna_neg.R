library(GenomicRanges)

# 初始化空向量来保存每条 circrna_data 的 m6a 到 ORF 的距离
dist_to_start_list <- vector("list", length(circrna_m6a ))
dist_to_end_list <- vector("list", length(circrna_m6a ))

# 遍历 circrna_data
for (i in seq_along(circrna_m6a )) {
  # 确保当前 circrna_data 这一行的 seqnames 是字符型
  seqname <- as.character(seqnames(circrna_m6a[i]))
  
  # 获取当前 circrna_data 的 ORF 范围，检查 ORF 列是否为有效数据
  orf_range_raw <- as.character(circrna_m6a$ORF[i])
  
  if (orf_range_raw == "0") {
    # 如果 ORF 列为 "0"，跳过这行数据
    dist_to_start_list[[i]] <- NA
    dist_to_end_list[[i]] <- NA
    next
  }
  
  # 如果 ORF 有值，拆分为多个 ORF 范围
  orf_ranges <- strsplit(orf_range_raw, ";")[[1]]
  orf_ranges <- lapply(orf_ranges, function(x) {
    range_values <- as.integer(strsplit(x, "-")[[1]])
    IRanges(start = range_values[1], end = range_values[2])
  })
  
  # 创建空的向量来存储每个 m6a 到 ORF 的距离
  m6a_to_start <- c()
  m6a_to_end <- c()
  
  # 获取当前 circrna_data 对应的 m6a 数据
  m6a_data <- m6a_circ_olp_neg[seqnames(m6a_circ_olp_neg) == seqname & 
                             start(m6a_circ_olp_neg) >= start(circrna_m6a[i]) & 
                             end(m6a_circ_olp_neg) <= end(circrna_m6a[i])]
  
  # 遍历每个 m6a 数据并计算与 ORF 距离
  for (j in seq_along(m6a_data)) {
    m6a <- m6a_data[j]
    m6a_start <- start(m6a)
    m6a_end <- end(m6a)
    
    # 遍历 ORF 范围
    for (orf_range in orf_ranges) {
      orf_start <- start(orf_range)
      orf_end <- end(orf_range)
      
      # 如果 m6a 在 ORF 范围内，计算距离
      if (m6a_start >= orf_start && m6a_end <= orf_end) {
        # 计算 m6a 到 ORF start 和 end 的距离
        dist_to_start <- m6a_start - orf_start
        dist_to_end <- m6a_end - orf_end
        
        # 将距离添加到向量中
        m6a_to_start <- c(m6a_to_start, dist_to_start)
        m6a_to_end <- c(m6a_to_end, dist_to_end)
      }
    }
  }
  
  # 如果没有计算到任何距离，添加 NA
  if (length(m6a_to_start) == 0) {
    dist_to_start_list[[i]] <- NA
    dist_to_end_list[[i]] <- NA
  } else {
    # 将计算得到的距离添加到列表中，使用逗号分隔
    dist_to_start_list[[i]] <- paste(m6a_to_start, collapse = ",")
    dist_to_end_list[[i]] <- paste(m6a_to_end, collapse = ",")
  }
}

7# 将距离列添加到 circrna_data 中
circrna_m6a$Dist_to_start <- unlist(dist_to_start_list)
circrna_m6a$Dist_to_end <- unlist(dist_to_end_list)

# 查看结果
head(circrna_m6a )
