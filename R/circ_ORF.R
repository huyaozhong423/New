library(Biostrings)
library(dplyr)
library(stringr)
library(GenomicRanges)

# Step 1: Read ORF fasta
orf_fasta <- readDNAStringSet("longest_orfs.cds")

# Step 2: Parse ORF information
orf_info <- tibble(name = names(orf_fasta)) %>%
  mutate(
    last_field = str_extract(name, "[^ ]+$"),
    id = str_extract(last_field, "^[^:]+"),
    start = as.integer(str_extract(last_field, "(?<=:)\\d+")),
    end = as.integer(str_extract(last_field, "(?<=-)\\d+")),
    strand = str_extract(last_field, "[+-](?=\\))")
  )

# Step 3: Keep positive strand
orf_info_positive <- orf_info %>%
  filter(strand == "+")

# Step 4: Aggregate ORF per circRNA
orf_aggregated <- orf_info_positive %>%
  group_by(id) %>%
  summarise(
    ORF = paste0(start, "-", end, collapse = ";"),
    ORF_Count = n()
  )

# Step 5: Prepare circRNA table
circrna_positive_df <- as.data.frame(circrna_positive) %>%
  dplyr::select(seqnames, start, end, strand, id, sequence, label)

# Step 6: Merge circRNA and ORF information
circrna_data <- circrna_positive_df %>%
  left_join(orf_aggregated, by = "id") %>%
  mutate(
    ORF = ifelse(is.na(ORF), NA, ORF),
    ORF_Count = ifelse(is.na(ORF_Count), 0, ORF_Count)
  )

# Step 7: Convert to GRanges
circrna_data_gr <- GRanges(
  seqnames = circrna_data$seqnames,
  ranges = IRanges(start = circrna_data$start, end = circrna_data$end),
  strand = circrna_data$strand,
  id = circrna_data$id,
  sequence = circrna_data$sequence,
  label = circrna_data$label,
  ORF = circrna_data$ORF,
  ORF_Count = circrna_data$ORF_Count
)

# Step 8: Adjust ORF to genomic coordinates
circrna_data_fixed <- circrna_data_gr

orf_absolute <- sapply(seq_along(circrna_data_fixed), function(i) {
  orf_entry <- circrna_data_fixed$ORF[i]
  if (is.na(orf_entry)) return(NA)
  
  circ_start <- start(circrna_data_fixed[i])
  circ_strand <- as.character(strand(circrna_data_fixed[i]))
  
  orf_ranges <- str_split(orf_entry, ";")[[1]]
  
  orf_genomic_ranges <- sapply(orf_ranges, function(rg) {
    rg_split <- str_split(rg, "-")[[1]]
    orf_start <- as.integer(rg_split[1])
    orf_end <- as.integer(rg_split[2])
    
    if (circ_strand == "+") {
      abs_start <- circ_start + orf_start - 1
      abs_end <- circ_start + orf_end - 1
    } else if (circ_strand == "-") {
      circ_end <- end(circrna_data_fixed[i])
      abs_start <- circ_end - orf_end + 1
      abs_end <- circ_end - orf_start + 1
    } else {
      stop("Unknown strand type")
    }
    
    paste0(abs_start, "-", abs_end)
  })
  
  paste0(orf_genomic_ranges, collapse = ";")
})

mcols(circrna_data_fixed)$ORF <- orf_absolute

# Output
circrna_data_fixed
