## pooja.singh09@gmail 
## re-analysing LD based on FST matched background 
## top 5% FST SNPs used for random background permutations and also for the top rwas snps
## oct 2026 revisions (with >100kb / interchromosomal physical linkage filter for BOTH background and top outliers)


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(svglite)
})

setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ld")

label_map <- c(
  "entantleter_MG_vs_newdegeni_MG"        = "E. antleter vs P. 'new degeni' Mwanza Gulf",
  "ecoprologous_blue_MG_vs_eparopius_MG"  = "E. paropius vs E. coprologous 'blue' Mwanza Gulf",
  "ecinctus_MG_vs_ecoprologous_blue_MG"   = "E. cinctus vs E. coprologous 'blue' Mwanza Gulf",
  "yplumbus_vs_ypyrrocephalus"            = "Y. plumbus vs Y. pyrrocephalus Makobe island",
  "lscraper_makobe_vs_lyellowchin_makobe" = "L. 'scraper' vs L. 'yellow chin' Makobe island",
  "nomni_makobe_vs_nuniscaper_makobe"     = "N. omnicaeruleus vs N. 'uniscraper' Makobe island",
  "ngigas_makobe_vs_paracyaneus_makobe"   = "N. gigas vs P. cyaneus Makobe island",
  "mmbipi_makobe_vs_ppinkanal_makobe"     = "M. mbipi vs P. 'pink anal' Makobe island",
  "pnye_kiss_vs_ppun_kiss"                = "P. 'nyererei-like' vs P. 'pundamilia-like' Kissenda island",
  "pnye_python_vs_ppun_python"            = "P. 'nyererei-like' vs P. 'pundamilia-like' Python island",
  "pnye_ruti_vs_ppun_ruti"                = "P. nyererei vs P. pundamilia Ruti island",
  "pnye_makobe_vs_ppun_makobe"            = "P. nyererei vs P. pundamilia Makobe island"
)

sample_file_map <- c(
  "entantleter_MG_vs_newdegeni_MG"        = "eantleter_pnewdegeni",
  "ecoprologous_blue_MG_vs_eparopius_MG"  = "eparopius_ecoprologous",
  "ecinctus_MG_vs_ecoprologous_blue_MG"   = "ecinctusE_ecorpologousblue2",
  "yplumbus_vs_ypyrrocephalus"            = "ypyrrocephalus_yplumbus",
  "lscraper_makobe_vs_lyellowchin_makobe" = "lscraper_lyellowchin",
  "nomni_makobe_vs_nuniscaper_makobe"     = "nomni_nuni",
  "ngigas_makobe_vs_paracyaneus_makobe"   = "ngigas_pcyaneus",
  "mmbipi_makobe_vs_ppinkanal_makobe"     = "mmbipi_ppinkanal",
  "pnye_kiss_vs_ppun_kiss"                = "pundamilias_kiss",
  "pnye_python_vs_ppun_python"            = "pundamilias_python",
  "pnye_ruti_vs_ppun_ruti"                = "pundamilias_ruti",
  "pnye_makobe_vs_ppun_makobe"            = "pundamilias_mak"
)

target_comparisons <- names(label_map)

cat("Loading combined Top LD file...\n")
master_ld_file <- fread("LD_r2_results_allpairs_combined.txt", header = TRUE, data.table = FALSE)

vcf_path <- "/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ld/all.SNPs.filt2.colour.imp.random1m.recode.vcf"

con <- file(vcf_path, "r")
skip_lines <- 0
while(TRUE) {
  line <- readLines(con, n = 1)
  if (length(line) == 0 || grepl("^#CHROM", line)) break
  skip_lines <- skip_lines + 1
}
close(con)

cat("Loading VCF metadata...\n")
my.data.full <- fread(vcf_path, skip = skip_lines, header = TRUE, data.table = FALSE, check.names = FALSE, comment.char = "")

chrom_col <- if ("#CHROM" %in% colnames(my.data.full)) "#CHROM" else "CHROM"
my.data.full$POS <- as.numeric(my.data.full$POS)
my.data.full$snp_id <- paste(my.data.full[[chrom_col]], my.data.full$POS, sep = "_")

all_fst_files <- list.files(path = ".", pattern = "_20kb\\.windowed\\.weir\\.fst$", full.names = TRUE)
fst_files <- all_fst_files[sapply(all_fst_files, function(f) {
  sub("_20kb\\.windowed\\.weir\\.fst$", "", basename(f)) %in% target_comparisons
})]

n_permutations <- 100
sample_size <- 50  
master_stats_ld <- data.frame()
master_stats_fst <- data.frame()
all_combined_results <- list()
all_fst_data <- list()

get_pair_fst_dt <- function(snp_data, fst_dt, col_name) {
  if (nrow(snp_data) == 0) return(rep(NA, nrow(snp_data)))
  snp_data$row_id <- seq_len(nrow(snp_data))
  
  s1_dt <- setDT(snp_data %>% select(row_id, CHROM = SNP1chr, POS = SNP1pos))
  setkey(s1_dt, CHROM, POS, POS)
  m1 <- fst_dt[s1_dt, on = .(CHROM, BIN_START <= POS, BIN_END >= POS), nomatch = NA]
  
  s2_dt <- setDT(snp_data %>% select(row_id, CHROM = SNP2chr, POS = SNP2pos))
  setkey(s2_dt, CHROM, POS, POS)
  m2 <- fst_dt[s2_dt, on = .(CHROM, BIN_START <= POS, BIN_END >= POS), nomatch = NA]
  
  res_dt <- data.table(row_id = snp_data$row_id)
  res_dt$FST1 <- m1[match(res_dt$row_id, m1$row_id)]$FST_VAL
  res_dt$FST2 <- m2[match(res_dt$row_id, m2$row_id)]$FST_VAL
  
  return((res_dt$FST1 + res_dt$FST2) / 2)
}

for (f in fst_files) {
  raw_comp <- sub("_20kb\\.windowed\\.weir\\.fst$", "", basename(f))
  cat("\nProcessing Unmatched Top-5% Background for comparison:", raw_comp, "\n")
  
  sample_txt <- sample_file_map[raw_comp]
  if (!file.exists(sample_txt)) next
  comparison_samples <- scan(sample_txt, what = character(), quiet = TRUE)
  valid_geno_cols <- intersect(comparison_samples, colnames(my.data.full))
  if (length(valid_geno_cols) == 0) next
  
  fst_df <- fread(f, header = TRUE, data.table = TRUE)
  fst_col <- if ("MEAN_FST" %in% colnames(fst_df)) "MEAN_FST" else colnames(fst_df)[grep("FST", colnames(fst_df), ignore.case = TRUE)[1]]
  setnames(fst_df, fst_col, "FST_VAL")
  fst_df$FST_VAL[fst_df$FST_VAL < 0] <- 0
  setDT(fst_df)
  setkey(fst_df, CHROM, BIN_START, BIN_END)
  
  target_top_string <- paste0(raw_comp, "_20kb.windowed.weir.fst_top")
  top_df <- master_ld_file %>% filter(comparison == target_top_string)
  if (nrow(top_df) == 0) next
  
  top_df <- top_df %>% filter((SNP1chr != SNP2chr) | (SNP1chr == SNP2chr & abs(SNP1pos - SNP2pos) > 100000))
  if (nrow(top_df) == 0) next
  
  top_df$Mean_FST <- get_pair_fst_dt(top_df, fst_df, "FST_VAL")
  top_df <- top_df %>% filter(!is.na(Mean_FST))
  top_df$comparison <- paste0(raw_comp, "_top")
  
  threshold_95 <- quantile(fst_df$FST_VAL, probs = 0.95, na.rm = TRUE)
  top5_windows <- fst_df %>% filter(FST_VAL >= threshold_95)
  if (nrow(top5_windows) == 0) next
  
  filtered_list <- list()
  for (chr in unique(top5_windows$CHROM)) {
    vcf_chr <- my.data.full[my.data.full[[chrom_col]] == chr, ]
    win_chr <- top5_windows[top5_windows$CHROM == chr, ]
    if (nrow(vcf_chr) > 0 && nrow(win_chr) > 0) {
      in_window <- logical(nrow(vcf_chr))
      for (i in seq_len(nrow(win_chr))) {
        in_window <- in_window | (vcf_chr$POS >= win_chr$BIN_START[i] & vcf_chr$POS <= win_chr$BIN_END[i])
      }
      filtered_list[[chr]] <- vcf_chr[in_window, ]
    }
  }
  my.data.filtered <- bind_rows(filtered_list) %>% distinct(snp_id, .keep_all = TRUE)
  if (nrow(my.data.filtered) < sample_size) next
  
  geno_mat <- as.matrix(my.data.filtered[, valid_geno_cols, drop = FALSE])
  suppressWarnings({
    hap_mat <- cbind(
      matrix(as.numeric(substr(geno_mat, 1, 1)), nrow = nrow(geno_mat)),
      matrix(as.numeric(substr(geno_mat, 3, 3)), nrow = nrow(geno_mat))
    )
  })
  
  all_perm_results <- list()
  set.seed(42)
  for (perm in 1:n_permutations) {
    samp_idx <- sample(1:nrow(hap_mat), sample_size)
    samp_haps <- t(hap_mat[samp_idx, ])
    
    d_matrix <- cov(samp_haps, use = "pairwise.complete.obs")
    rsq_matrix <- cor(samp_haps, use = "pairwise.complete.obs")^2
    
    upper_tri_idx <- which(upper.tri(rsq_matrix), arr.ind = TRUE)
    snp1_meta <- my.data.filtered[samp_idx[upper_tri_idx[, 1]], ]
    snp2_meta <- my.data.filtered[samp_idx[upper_tri_idx[, 2]], ]
    
    iter_df <- data.frame(
      SNP1 = snp1_meta$snp_id, SNP2 = snp2_meta$snp_id,
      SNP1chr = snp1_meta[[chrom_col]], SNP1pos = snp1_meta$POS,
      SNP2chr = snp2_meta[[chrom_col]], SNP2pos = snp2_meta$POS,
      D = d_matrix[upper_tri_idx], Rsq = rsq_matrix[upper_tri_idx],
      comparison = paste0(raw_comp, "_bg"), stringsAsFactors = FALSE
    )
    
    iter_df <- iter_df %>% filter((SNP1chr != SNP2chr) | (SNP1chr == SNP2chr & abs(SNP1pos - SNP2pos) > 100000))
    all_perm_results[[perm]] <- iter_df[!is.na(iter_df$Rsq), ]
  }
  
  master_perm_df <- bind_rows(all_perm_results)
  master_perm_df_unique <- master_perm_df[!duplicated(master_perm_df[, c("SNP1", "SNP2")]), ]
  master_perm_df_unique$Mean_FST <- get_pair_fst_dt(master_perm_df_unique, fst_df, "FST_VAL")
  master_perm_df_unique <- master_perm_df_unique %>% filter(!is.na(Mean_FST))
  
  # NO STRICT BIN MATCHING: Use the full unique top-5% background pool directly
  matched_bg_df <- master_perm_df_unique
  
  if (nrow(matched_bg_df) > 0) {
    combined_results <- rbind(
      matched_bg_df %>% select(SNP1, SNP2, SNP1chr, SNP1pos, SNP2chr, SNP2pos, D, Rsq, comparison),
      top_df %>% select(SNP1, SNP2, SNP1chr, SNP1pos, SNP2chr, SNP2pos, D, Rsq, comparison)
    )
    
    fwrite(combined_results, paste0(raw_comp, "_LD_r2_permuted_results.txt"), sep = "\t", quote = FALSE)
    all_combined_results[[raw_comp]] <- combined_results
    
    bg_rsq <- as.numeric(combined_results$Rsq[grepl("_bg$", combined_results$comparison)])
    top_rsq <- as.numeric(combined_results$Rsq[grepl("_top$", combined_results$comparison)])
    
    if (length(bg_rsq) > 1 && length(top_rsq) > 1) {
      stats_ld <- wilcox.test(bg_rsq, top_rsq)
      sig_ld <- ifelse(stats_ld$p.value < 0.001, "***", ifelse(stats_ld$p.value < 0.01, "**", ifelse(stats_ld$p.value < 0.05, "*", "ns")))
      master_stats_ld <- rbind(master_stats_ld, data.frame(Comparison = raw_comp, sig_label = sig_ld))
    }
    
    matched_bg_fst <- matched_bg_df$Mean_FST
    top_fst_vals <- top_df$Mean_FST
    if (length(matched_bg_fst) > 0 && length(top_fst_vals) > 0) {
      comp_fst_df <- rbind(
        data.frame(Raw_Comp = raw_comp, FST = matched_bg_fst, Set = "Genomic background"),
        data.frame(Raw_Comp = raw_comp, FST = top_fst_vals, Set = "Top RWAS & FST outliers")
      )
      all_fst_data[[raw_comp]] <- comp_fst_df
      
      stats_fst <- wilcox.test(matched_bg_fst, top_fst_vals)
      sig_fst <- ifelse(stats_fst$p.value < 0.001, "***", ifelse(stats_fst$p.value < 0.01, "**", ifelse(stats_fst$p.value < 0.05, "*", "ns")))
      master_stats_fst <- rbind(master_stats_fst, data.frame(Raw_Comp = raw_comp, sig_label = sig_fst))
    }
  }
}

# ==============================================================================
# GENERATE FIGURE A & FIGURE B PLOTS
# ==============================================================================
cat("\nGenerating Figure A & Figure B SVGs...\n")

master_ld_df <- bind_rows(all_combined_results)
master_ld_df$Set <- factor(ifelse(grepl("_top$", master_ld_df$comparison), "Top RWAS & FST outliers", "Genomic background"), levels = c("Top RWAS & FST outliers", "Genomic background"))
master_ld_df$Raw_Comp <- gsub("_top$|_bg$", "", master_ld_df$comparison)
master_ld_df$Formatted_Name <- factor(label_map[master_ld_df$Raw_Comp], levels = rev(unname(label_map)))

star_positions_ld <- master_ld_df %>% group_by(Formatted_Name, Raw_Comp) %>% summarise(max_val = max(Rsq, na.rm = TRUE), .groups = "drop") %>% left_join(master_stats_ld, by = c("Raw_Comp" = "Comparison"))

p_ld <- ggplot(master_ld_df, aes(x = Formatted_Name, y = Rsq, fill = Set)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.7, colour = "black", linewidth = 0.4) +
  geom_text(data = star_positions_ld, aes(x = Formatted_Name, y = max_val + 0.05, label = sig_label), inherit.aes = FALSE, size = 4, angle = -90, vjust = 0.5) +
  coord_flip(ylim = c(0, 1.1)) + theme_bw() +
  labs(x = NULL, y = bquote("Linkage Disequilibrium ("~r^2~")"), fill = NULL) +
  theme(
    axis.text.y = element_text(size = 10, colour = "black", face = "italic"),
    axis.text.x = element_text(size = 10, colour = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 12, colour = "black", margin = margin(t = 10)),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.9), legend.background = element_blank(), legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_fill_manual(values = c("Top RWAS & FST outliers" = "darkgrey", "Genomic background" = "white")) +
  annotate("text", x = 12.5, y = -0.05, label = "A", size = 8, fontface = "bold")

ggsave("Figure_A_LD_Comparison_BGtop5.svg", plot = p_ld, width = 9, height = 7)

master_fst_df <- bind_rows(all_fst_data)
master_fst_df$Set <- factor(master_fst_df$Set, levels = c("Top RWAS & FST outliers", "Genomic background"))
master_fst_df$Formatted_Name <- factor(label_map[master_fst_df$Raw_Comp], levels = rev(unname(label_map)))

star_positions_fst <- master_fst_df %>% group_by(Formatted_Name, Raw_Comp) %>% summarise(max_val = max(FST, na.rm = TRUE), .groups = "drop") %>% left_join(master_stats_fst, by = "Raw_Comp")

p_fst <- ggplot(master_fst_df, aes(x = Formatted_Name, y = FST, fill = Set)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.7, colour = "black", linewidth = 0.4) +
  geom_text(data = star_positions_fst, aes(x = Formatted_Name, y = max_val + 0.02, label = sig_label), inherit.aes = FALSE, size = 4, angle = -90, vjust = 0.5) +
  coord_flip() + theme_bw() +
  labs(x = NULL, y = expression(paste("Fixation Index (", F[ST], ")")), fill = NULL) +
  theme(
    axis.text.y = element_text(size = 10, colour = "black", face = "italic"),
    axis.text.x = element_text(size = 10, colour = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 12, colour = "black", margin = margin(t = 10)),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.9), legend.background = element_blank(), legend.key = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_fill_manual(values = c("Top RWAS & FST outliers" = "darkgrey", "Genomic background" = "white")) +
  annotate("text", x = 12.5, y = -0.02, label = "B", size = 8, fontface = "bold")

ggsave("Figure_B_FST_Comparison_BGtop5.svg", plot = p_fst, width = 9, height = 7)
cat("Pipeline complete! Figures saved.\n")




### ### correlation between genome-wide LD and fst interchormosome and 100kb chromsomal - FST MATCHED DIST

library(dplyr)
library(data.table)
library(ggplot2)
library(svglite)

setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ld")

ld_files <- list.files(path = ".", pattern = "_LD_r2_permuted_results\\.txt$", full.names = TRUE)
bg_list <- list()
out_list <- list()

cat("Processing real data for dual-layer correlation plot...\n")

for (ld_f in ld_files) {
  raw_comp <- gsub("_LD_r2_permuted_results\\.txt$", "", basename(ld_f))
  fst_file <- paste0(raw_comp, "_20kb.windowed.weir.fst")
  
  if (!file.exists(fst_file)) next
  
  # 1. Load combined LD results file
  ld_df <- fread(ld_f, header = TRUE, data.table = FALSE)
  
  # 2. Load windowed FST for on-the-fly coordinate mapping
  fst_df <- fread(fst_file, header = TRUE, data.table = TRUE)
  fst_col <- if ("MEAN_FST" %in% colnames(fst_df)) "MEAN_FST" else colnames(fst_df)[grep("FST", colnames(fst_df), ignore.case = TRUE)[1]]
  setnames(fst_df, fst_col, "FST_VAL")
  fst_df$FST_VAL[fst_df$FST_VAL < 0] <- 0
  setDT(fst_df)
  setkey(fst_df, CHROM, BIN_START, BIN_END)
  
  # Helper function to map FST to pairs
  map_fst_to_pairs <- function(sub_df) {
    if (nrow(sub_df) == 0) return(sub_df)
    sub_df$row_id <- seq_len(nrow(sub_df))
    
    s1_dt <- setDT(sub_df %>% select(row_id, CHROM = SNP1chr, POS = SNP1pos))
    setkey(s1_dt, CHROM, POS, POS)
    m1 <- fst_df[s1_dt, on = .(CHROM, BIN_START <= POS, BIN_END >= POS), nomatch = NA]
    
    s2_dt <- setDT(sub_df %>% select(row_id, CHROM = SNP2chr, POS = SNP2pos))
    setkey(s2_dt, CHROM, POS, POS)
    m2 <- fst_df[s2_dt, on = .(CHROM, BIN_START <= POS, BIN_END >= POS), nomatch = NA]
    
    res_dt <- data.table(row_id = sub_df$row_id)
    res_dt$FST1 <- m1[match(res_dt$row_id, m1$row_id)]$FST_VAL
    res_dt$FST2 <- m2[match(res_dt$row_id, m2$row_id)]$FST_VAL
    
    sub_df$geom_mean_fst <- sqrt(res_dt$FST1 * res_dt$FST2)
    return(sub_df)
  }
  
  # Process Background
  bg_sub <- ld_df %>% filter(grepl("_bg$", comparison))
  if (nrow(bg_sub) > 0) {
    bg_sub <- map_fst_to_pairs(bg_sub)
    clean_bg <- bg_sub %>% filter(!is.na(Rsq) & !is.na(geom_mean_fst)) %>% select(geom_mean_fst, ld_r2 = Rsq)
    clean_bg$locus_set <- "Permuted Background"
    bg_list[[raw_comp]] <- clean_bg
  }
  
  # Process Outliers
  out_sub <- ld_df %>% filter(grepl("_top$", comparison))
  if (nrow(out_sub) > 0) {
    out_sub <- map_fst_to_pairs(out_sub)
    clean_out <- out_sub %>% filter(!is.na(Rsq) & !is.na(geom_mean_fst)) %>% select(geom_mean_fst, ld_r2 = Rsq)
    clean_out$locus_set <- "GWAS Outliers"
    out_list[[raw_comp]] <- clean_out
  }
}

bg_df <- bind_rows(bg_list)
out_df <- bind_rows(out_list)
df <- bind_rows(bg_df, out_df)

if (nrow(bg_df) == 0 || nrow(out_df) == 0) {
  stop("Error: Missing background or outlier data. Verify input files.")
}

# 2. Calculate statistics dynamically for the annotation box
rho_bg  <- cor(bg_df$geom_mean_fst, bg_df$ld_r2, method = "spearman", use = "complete.obs")
rho_out <- cor(out_df$geom_mean_fst, out_df$ld_r2, method = "spearman", use = "complete.obs")
fit     <- lm(ld_r2 ~ geom_mean_fst * locus_set, data = df)
p_inter <- summary(fit)$coefficients[grepl("geom_mean_fst:locus_set", rownames(summary(fit)$coefficients)), "Pr(>|t|)"]
if (length(p_inter) == 0) p_inter <- 1.0

stat_label <- sprintf(
  "Background rho = %.3f\nGWAS Outliers rho = %.3f\nSlope interaction p = %.2e",
  rho_bg, rho_out, p_inter
)

# 3. Build the plot
cat("Generating dual-layer correlation plot (SVG)...\n")

p_dual <- ggplot() +
  # 2D binned heatmap for Permuted Background
  stat_bin_2d(
    data = bg_df,
    aes(x = geom_mean_fst, y = ld_r2, fill = after_stat(count)),
    bins = 45
  ) +
  scale_fill_viridis_c(
    option = "mako",
    trans = "log10",
    name = "Background\nPair Count (log10)"
  ) +
  
  # Scatter points for GWAS Outliers
  geom_point(
    data = out_df,
    aes(x = geom_mean_fst, y = ld_r2),
    color = "#E41A1C",
    alpha = 0.45,
    size = 1.6
  ) +
  
  # Linear regression trend lines for both groups
  geom_smooth(
    data = df,
    aes(x = geom_mean_fst, y = ld_r2, color = locus_set, linetype = locus_set),
    method = "lm",
    se = TRUE,
    linewidth = 1.2
  ) +
  scale_color_manual(
    name = "Locus Set",
    values = c("Permuted Background" = "#1E90FF", "GWAS Outliers" = "#E41A1C")
  ) +
  scale_linetype_manual(
    name = "Locus Set",
    values = c("Permuted Background" = "dashed", "GWAS Outliers" = "solid")
  ) +
  
  # Top-left statistical summary box
  annotate(
    "label",
    x = 0.02, y = 0.98,
    hjust = 0, vjust = 1,
    label = stat_label,
    size = 3.6,
    fill = "white",
    color = "black",
    fontface = "bold",
    label.size = 0.3
  ) +
  
  # Titles, axis labels, and theme styling
  labs(
    title = "Unlinked LD vs. FST: GWAS Outliers vs. Genomic Background",
    subtitle = "Interchromosomal & >100 kb SNP pairs across sister species",
    x = expression(Geometric~Mean~Windowed~F[ST]~(sqrt(F[ST1] %*% F[ST2]))),
    y = expression(Linkage~Disequilibrium~(r^2))
  ) +
  coord_cartesian(xlim = c(0, 0.85), ylim = c(0, 1.02)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray25"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

print(p_dual)
ggsave("Dual_Layer_LD_vs_FST_Correlation.svg", plot = p_dual, width = 9, height = 6)
cat("Dual-layer correlation plot successfully saved to Dual_Layer_LD_vs_FST_Correlation.svg!\n")







