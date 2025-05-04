human_ires_info <- read.delim("Human_ires_info.txt", header = TRUE, stringsAsFactors = FALSE)

ires_positions <- unlist(human_ires_info[, 3])  # 假设 IRES 信息在第3列

# 查看更新后的 circRNA 数据
head(circrna_data)
# 清理逗号，确保位置是数字
ires_positions_clean <- gsub(",", "", ires_positions)
# 去掉包含 NA 的位置
ires_positions_clean <- ires_positions_clean[!is.na(start_positions) & !is.na(end_positions)]


# 再次尝试解析 IRES 的位置
ires_granges <- GRanges(
  seqnames = sub(":(.*)", "", ires_positions_clean),  # 提取染色体信息
  ranges = IRanges(
    start = as.numeric(sub(".*:(\\d+)-(\\d+)", "\\1", ires_positions_clean)),  # 提取起始位置
    end = as.numeric(sub(".*:(\\d+)-(\\d+)", "\\2", ires_positions_clean))  # 提取结束位置
  ),
  location = paste0("IRES_", seq_along(ires_positions_clean))  # 添加一个标识每个 IRES 的位置
)

# 检查转换后的 GRanges 对象
head(ires_granges)

# 第二步：使用 findOverlaps 来匹配 IRES 和 circRNA 数据
overlaps <- findOverlaps(circrna_data, ires_granges)

# 第三步：提取匹配的 IRES 位置和数量
ires_location <- ires_granges$location[subjectHits(overlaps)]  # 获取匹配的 IRES 位置
ires_count <- rep(1, length(ires_location))  # 如果每个 IRES 匹配计数为 1，可以根据需求调整

# 在 circRNA 数据中添加新的列
circrna_data$IRES_Location <- NA
circrna_data$IRES_Count <- 0
circrna_data$IRES_Location[queryHits(overlaps)] <- ires_location
circrna_data$IRES_Count[queryHits(overlaps)] <- ires_count

overlap_results <- findOverlaps(m6a_circ_olp_updated, ires_granges)

# 初始化新的列来保存匹配结果
m6a_circ_olp_updated$m6a_in_ires <- 0  # 默认为 0，即不在 IRES 区域内
m6a_circ_olp_updated$dist_to_ires_start <- NA  # 默认为 NA
m6a_circ_olp_updated$dist_to_ires_end <- NA  # 默认为 NA

# 如果匹配到重叠区域，更新相关列
m6a_circ_olp_updated$m6a_in_ires[subjectHits(overlap_results)] <- 1

# 计算 m6a 到 IRES 两端的距离
m6a_circ_olp_updated$dist_to_ires_start[subjectHits(overlap_results)] <- 
  start(ires_granges[queryHits(overlap_results)]) - start(m6a_circ_olp_updated[subjectHits(overlap_results)])

m6a_circ_olp_updated$dist_to_ires_end[subjectHits(overlap_results)] <- 
  end(ires_granges[queryHits(overlap_results)]) - end(m6a_circ_olp_updated[subjectHits(overlap_results)])

# 查看更新后的数据
head(m6a_circ_olp_updated)
