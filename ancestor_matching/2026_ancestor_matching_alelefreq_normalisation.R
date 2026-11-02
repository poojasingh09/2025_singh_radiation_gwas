### ancestor matching 2026 oct by pooja singh
##t.haryngalis and h.gracilior  make the nile lineage
## A stappesi and A yaekama make the congo lineage
## null distribution calcualted from allele frequency normalised background
## background allelele freq is normalised by the distritbution of MAF of top RWAS outliers
## not sure how much sense this makese before then obviously you will get less enrichment since most of the
## genomes of lvrs species is from congo and nile, and rwas allales are high freq, so if most background is high freq 
## then it will be enriched for congo/nile both ancestry too
## enrichment as calculated as > 75th and >95th percentile of null dist.

# ============================================================
# ANCESTRY ENRICHMENT ANALYSIS
# MAF-MATCHED RANDOM BACKGROUND
# ============================================================

setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/2026_ancestor_matching/")

# ------------------------------------------------------------
# 1. LIBRARIES
# ------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)
library(stringr)
library(ggpattern)

theme_set(theme_pubr())

# ------------------------------------------------------------
# 2. READ OBSERVED TOP RWAS SNPs
# ------------------------------------------------------------

af <- read.table(
  "../ancestor_matching/top_rwas_all_99.9_and_pip01.allthegoodstuff.txt",
  header = TRUE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 3. READ ANCESTRAL GENOTYPES
# ------------------------------------------------------------

all <- read.table(
  "allchr.anc.rwas.GT.2026.GT.numeric",
  header = TRUE,
  check.names = FALSE
)

all$id <- paste(
  gsub("^chr", "", all$CHROM),
  all$POS,
  sep = "-"
)

all1 <- all[, c(
  "id", "71001", "71014", "81032",
  "64538", "64539", "64767", "64843",
  "64860", "81342", "81343", "81344",
  "81345", "81309"
)]

all_af <- merge(af, all1, by = "id")

all_af1 <- unique(all_af)

all_af1$special <- paste0(
  all_af1$chr, "_", all_af1$start.x
)

all_af2 <- all_af1[
  !duplicated(all_af1$special),
]

rownames(all_af2) <- all_af2$special

# ------------------------------------------------------------
# 4. IDENTIFY MINOR AND MAJOR ALLELES
# ------------------------------------------------------------

all_af1$major_af <- ifelse(
  all_af1$ref_af > all_af1$alt_af,
  all_af1$ref_af,
  ifelse(
    all_af1$ref_af < all_af1$alt_af,
    all_af1$alt_af,
    NA_real_
  )
)

all_af1$major_a <- ifelse(
  all_af1$ref_af > all_af1$alt_af,
  all_af1$ref,
  ifelse(
    all_af1$ref_af < all_af1$alt_af,
    all_af1$alt,
    NA_character_
  )
)

all_af1$minor_af <- ifelse(
  all_af1$ref_af < all_af1$alt_af,
  all_af1$ref_af,
  ifelse(
    all_af1$ref_af > all_af1$alt_af,
    all_af1$alt_af,
    NA_real_
  )
)

all_af1$minor_a <- ifelse(
  all_af1$ref_af < all_af1$alt_af,
  all_af1$ref,
  ifelse(
    all_af1$ref_af > all_af1$alt_af,
    all_af1$alt,
    NA_character_
  )
)

# Split according to which allele is minor

all_af1_a <- all_af1[
  !is.na(all_af1$minor_a) &
    all_af1$minor_a == all_af1$alt,
]

all_af1_r <- all_af1[
  !is.na(all_af1$minor_a) &
    all_af1$minor_a == all_af1$ref,
]

# ------------------------------------------------------------
# 5. ANCESTRAL GENOTYPE CLASSIFICATION
# ------------------------------------------------------------

tp_cols <- c(
  "71001", "71014", "81032",
  "64538", "64539", "64767",
  "64843", "64860"
)

as_cols <- c(
  "81342", "81343", "81344",
  "81345", "81309"
)

# Classify whether the RWAS minor allele is present.
# NA is retained if genotypes are incomplete and no carrier
# is identified.

check_allele <- function(data, columns, present_codes) {
  
  genotypes <- as.matrix(
    data[, columns, drop = FALSE]
  )
  
  apply(genotypes, 1, function(g) {
    
    g <- as.character(g)
    
    if (any(g %in% present_codes, na.rm = TRUE)) {
      return("yes")
    }
    
    if (anyNA(g) || any(g == "", na.rm = TRUE)) {
      return(NA_character_)
    }
    
    return("no")
  })
}

# Minor allele is ALT: genotype codes 1 and 2 carry ALT

all_af1_a$gwas_allele_TP <- check_allele(
  all_af1_a,
  tp_cols,
  c("1", "2")
)

all_af1_a$gwas_allele_AS <- check_allele(
  all_af1_a,
  as_cols,
  c("1", "2")
)

# Minor allele is REF: genotype codes 0 and 1 carry REF

all_af1_r$gwas_allele_TP <- check_allele(
  all_af1_r,
  tp_cols,
  c("0", "1")
)

all_af1_r$gwas_allele_AS <- check_allele(
  all_af1_r,
  as_cols,
  c("0", "1")
)

# Combine both allele groups

all_af2 <- rbind(
  all_af1_a,
  all_af1_r
)

# ------------------------------------------------------------
# 6. SUMMARISE OBSERVED TOP SNPs BY TRAIT
# ------------------------------------------------------------

traits <- c(
  "cs", "ir", "kd", "tm", "pc1", "pc2",
  "vb", "mlb", "dlb", "fb", "yf", "rc",
  "mb", "bb"
)

categories <- c(
  "ancestor_Tpharyngalis",
  "ancestor_Astappersi",
  "both",
  "neither"
)

Res <- list()

for (trait in traits) {
  
  all_af3 <- all_af2[
    all_af2$trait == trait,
  ]
  
  known <- !is.na(all_af3$gwas_allele_TP) &
    !is.na(all_af3$gwas_allele_AS)
  
  all_af3_known <- all_af3[known, ]
  
  tp_only <- sum(
    all_af3_known$gwas_allele_TP == "yes" &
      all_af3_known$gwas_allele_AS == "no"
  )
  
  as_only <- sum(
    all_af3_known$gwas_allele_TP == "no" &
      all_af3_known$gwas_allele_AS == "yes"
  )
  
  both <- sum(
    all_af3_known$gwas_allele_TP == "yes" &
      all_af3_known$gwas_allele_AS == "yes"
  )
  
  neither <- sum(
    all_af3_known$gwas_allele_TP == "no" &
      all_af3_known$gwas_allele_AS == "no"
  )
  
  counts <- c(
    tp_only,
    as_only,
    both,
    neither
  )
  
  total <- sum(counts)
  
  res <- data.frame(
    trait = trait,
    source = c(
      "ancestor_Tpharyngalis",
      "ancestor_Astappersi",
      "both",
      "neither",
      "total"
    ),
    value = c(counts, total)
  )
  
  res$percent <- NA_real_
  
  if (total > 0) {
    res$percent[1:4] <- counts / total * 100
    res$percent[5] <- 100
  }
  
  message(
    trait, ": ",
    sum(!known),
    " SNPs excluded because ancestry status is unknown."
  )
  
  print(res)
  
  Res[[trait]] <- res
}

Res1 <- do.call(rbind, Res)

# ------------------------------------------------------------
# 7. PLOT OBSERVED ANCESTRAL COMPOSITION
# ------------------------------------------------------------

trait_order <- traits

plot_data <- Res1 %>%
  filter(source != "total") %>%
  mutate(
    source = recode(
      source,
      "ancestor_Tpharyngalis" = "T. pharyngalis only",
      "ancestor_Astappersi" = "A. stappersi only",
      "both" = "Both",
      "neither" = "Neither"
    ),
    trait = factor(
      trait,
      levels = trait_order
    )
  ) %>%
  group_by(trait) %>%
  mutate(
    percent_plot = if (sum(value, na.rm = TRUE) > 0) {
      value / sum(value, na.rm = TRUE) * 100
    } else {
      NA_real_
    }
  ) %>%
  ungroup()

plot_data$source <- factor(
  plot_data$source,
  levels = c(
    "Neither",
    "Both",
    "T. pharyngalis only",
    "A. stappersi only"
  )
)

p_observed <- ggplot(
  plot_data,
  aes(
    x = trait,
    y = percent_plot,
    fill = source,
    pattern = source
  )
) +
  geom_col_pattern(
    width = 0.88,
    colour = "grey80",
    linewidth = 0.4,
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.25,
    pattern_spacing = 0.035,
    pattern_size = 0.3
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Neither" = "white",
      "Both" = "lightgrey",
      "T. pharyngalis only" = "white",
      "A. stappersi only" = "black"
    )
  ) +
  scale_pattern_manual(
    values = c(
      "Neither" = "none",
      "Both" = "none",
      "T. pharyngalis only" = "pch",
      "A. stappersi only" = "none"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_discrete(
    labels = function(x) {
      result <- toupper(x)
      result[x == "bb"] <- "BF"
      result
    }
  ) +
  labs(
    x = NULL,
    y = "Proportion of top SNPs (%)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      colour = "black",
      size = 11
    ),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    legend.position = "right"
  )

print(p_observed)



# ==============================================================================
# PIPELINE: ANCESTRAL ENRICHMENT (75th [blue] & 95th [red] Percentile Thresholds) : BG ALT ALLELE
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpattern)
  library(dplyr)
  library(tidyr)
  if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
  library(svglite)
})

# ------------------------------------------------------------------------------
# 1. READ DATA FILES
# ------------------------------------------------------------------------------

cat("Loading background ancestral genotypes...\n")
data1 <- read.table(
  "allchr.anc.100krandom.mafcorrected.GT.2026.GT.numeric",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Loading background mapping flock MAF file...\n")
bg_maf_df <- read.table(
  "random100k_maf_checked_snps2_with_MAF.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(bg_maf_df)[1:3] <- c("CHROM", "POS", "MAF")

# ------------------------------------------------------------------------------
# 2. COMPARE MAF: MAF-CORRECTED BACKGROUND VS. GWAS HITS
# ------------------------------------------------------------------------------

bg_maf_vals <- as.numeric(bg_maf_df$MAF)
bg_maf_vals <- bg_maf_vals[!is.na(bg_maf_vals) & bg_maf_vals >= 0 & bg_maf_vals <= 0.5]

gwas_maf_vals <- as.numeric(af$minor_af)
gwas_maf_vals <- gwas_maf_vals[!is.na(gwas_maf_vals) & gwas_maf_vals >= 0 & gwas_maf_vals <= 0.5]

cat("\n=== Sample Sizes ===\n")
cat("GWAS MAF count:      ", length(gwas_maf_vals), "\n")
cat("Background MAF count:", length(bg_maf_vals), "\n")

# Kolmogorov-Smirnov test
ks_res <- ks.test(gwas_maf_vals, bg_maf_vals)
cat("\n=== Kolmogorov-Smirnov Test ===\n")
print(ks_res)

comp_df <- rbind(
  data.frame(MAF = gwas_maf_vals, Set = "GWAS Top Hits"),
  data.frame(MAF = bg_maf_vals,   Set = "Matched Background (100k)")
)

p_maf_validation <- ggplot(comp_df, aes(x = MAF, fill = Set, colour = Set)) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 0.02,
    position = "identity",
    alpha = 0.45,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_density(linewidth = 0.9, fill = NA, adjust = 1.2) +
  scale_fill_manual(
    name = NULL,
    values = c("GWAS Top Hits" = "#1f77b4", "Matched Background (100k)" = "#ff7f0e")
  ) +
  scale_colour_manual(
    name = NULL,
    values = c("GWAS Top Hits" = "#0d47a1", "Matched Background (100k)" = "#b34700")
  ) +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05), limits = c(0, 0.51)) +
  labs(
    title = "Validation of Allele Frequency Matching in Mapping Flock",
    subtitle = paste0(
      "Kolmogorov-Smirnov D = ", round(ks_res$statistic, 4),
      ", p-value = ", signif(ks_res$p.value, 4)
    ),
    x = "Minor Allele Frequency (MAF)",
    y = "Density"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 11),
    axis.text = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 12, colour = "black")
  )

print(p_maf_validation)

# ------------------------------------------------------------------------------
# 3. CLASSIFY FULL UNTHINNED BACKGROUND SNPS INTO ANCESTRIES
# ------------------------------------------------------------------------------

tp_cols <- c("64538", "64539", "64767", "64843", "64860", "71001", "71014", "81032") # 8 TP (Nile)
as_cols <- c("81309", "81342", "81343", "81344", "81345")                             # 5 AS (Congo)

check_allele_bg <- function(df, tp_cols, as_cols) {
  tp_mat <- as.matrix(df[, tp_cols, drop = FALSE])
  as_mat <- as.matrix(df[, as_cols, drop = FALSE])
  storage.mode(tp_mat) <- "numeric"
  storage.mode(as_mat) <- "numeric"
  
  has_tp <- apply(tp_mat, 1, function(x) any(x > 0, na.rm = TRUE))
  has_as <- apply(as_mat, 1, function(x) any(x > 0, na.rm = TRUE))
  
  category <- ifelse(has_tp & has_as, "both",
                     ifelse(has_tp & !has_as, "ancestor_Tpharyngalis",
                            ifelse(!has_tp & has_as, "ancestor_Astappersi", "neither")))
  return(category)
}

data1$category <- check_allele_bg(data1, tp_cols, as_cols)
data1$MAF <- bg_maf_df$MAF[match(
  paste(data1$CHROM, data1$POS),
  paste(bg_maf_df$CHROM, bg_maf_df$POS)
)]

bg_pool <- data1[!is.na(data1$category) & !is.na(data1$MAF) & data1$MAF > 0 & data1$MAF <= 0.5, ]

cat("\n--- Unthinned Background Dimensions & Counts ---\n")
print(dim(bg_pool))
print(table(bg_pool$category))
print(round(prop.table(table(bg_pool$category)) * 100, 2))

# ------------------------------------------------------------------------------
# 4. RUN UNBIASED EMPIRICAL ENRICHMENT TEST (B = 10,000)
# ------------------------------------------------------------------------------

traits_ordered <- c(
  "cs", "ir", "kd", "tm", "pc1", "pc2", "vb", 
  "mlb", "dlb", "fb", "yf", "rc", "mb", "bb"
)

categories_all <- c("ancestor_Tpharyngalis", "ancestor_Astappersi", "both", "neither")
categories_report <- c("ancestor_Tpharyngalis", "ancestor_Astappersi", "both")

get_counts <- function(d) {
  setNames(
    as.numeric(table(factor(d$category, levels = categories_all))),
    categories_all
  )
}

set.seed(123)
B <- 10000

results <- list()
k <- 0

traits <- unique(Res1$trait)

for (tr in traits) {
  obs <- Res1[Res1$trait == tr, , drop = FALSE]
  counts <- setNames(obs$value[match(categories_all, obs$source)], categories_all)
  counts[is.na(counts)] <- 0
  n <- sum(counts)
  
  if (n == 0 || n > nrow(bg_pool)) {
    warning("Skipping trait: ", tr)
    next
  }
  
  observed_pct <- counts / n * 100
  
  null <- replicate(B, {
    samp <- bg_pool[sample.int(nrow(bg_pool), size = n, replace = FALSE), , drop = FALSE]
    get_counts(samp)
  })
  
  for (cat in categories_report) {
    k <- k + 1
    observed_count <- unname(counts[cat])
    expected_count <- mean(null[cat, ])
    
    # Calculate 75th and 95th percentiles of the null distribution
    null_75_count  <- unname(quantile(null[cat, ], 0.75))
    null_95_count  <- unname(quantile(null[cat, ], 0.95))
    
    p_value <- (1 + sum(null[cat, ] >= observed_count)) / (B + 1)
    
    results[[k]] <- data.frame(
      trait            = tr,
      category         = cat,
      n_observed       = n,
      observed_count   = observed_count,
      observed_percent = unname(observed_pct[cat]),
      expected_count   = expected_count,
      expected_percent = expected_count / n * 100,
      null_75_count    = null_75_count,
      null_95_count    = null_95_count,
      enrichment_ratio = if (expected_count > 0) observed_count / expected_count else NA_real_,
      p_value          = p_value
    )
  }
}

enrichment_unthinned <- do.call(rbind, results)

enrichment_unthinned <- enrichment_unthinned %>%
  group_by(category) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(category, FDR, p_value)

# ------------------------------------------------------------------------------
# 5. STACKED BAR PLOT WITH DUAL THRESHOLD ASTERISKS (Red = 75th, Blue = 95th)
# ------------------------------------------------------------------------------

stack_order <- c(
  "A. stappersi only",    # Congo
  "T. pharyngalis only",  # Nilotic
  "Both",
  "Neither"
)

plot_data <- Res1 %>%
  filter(source != "total") %>%
  mutate(
    source = recode(
      source,
      "ancestor_Tpharyngalis" = "T. pharyngalis only",
      "ancestor_Astappersi"   = "A. stappersi only",
      "both"                  = "Both",
      "neither"               = "Neither"
    ),
    source = factor(as.character(source), levels = stack_order),
    trait  = factor(tolower(trait), levels = traits_ordered)
  ) %>%
  filter(!is.na(trait)) %>%
  group_by(trait) %>%
  mutate(
    percent_plot = if (sum(value, na.rm = TRUE) > 0) {
      value / sum(value, na.rm = TRUE) * 100
    } else {
      NA_real_
    }
  ) %>%
  ungroup()

# Determine significance level for each bar segment
sig_data_classified <- enrichment_unthinned %>%
  mutate(
    source = recode(
      category,
      "ancestor_Astappersi"   = "A. stappersi only",
      "ancestor_Tpharyngalis" = "T. pharyngalis only",
      "both"                  = "Both",
      "neither"               = "Neither"
    ),
    trait = factor(tolower(trait), levels = traits_ordered),
    sig_level = case_when(
      enrichment_ratio > 1.0 & observed_count >= null_95_count ~ "95th",
      enrichment_ratio > 1.0 & observed_count >= null_75_count ~ "75th",
      TRUE ~ "NS"
    )
  ) %>%
  filter(sig_level != "NS") %>%
  dplyr::select(trait, source, sig_level)

sig_positions <- plot_data %>%
  mutate(
    source = as.character(source),
    stack_order_num = match(source, stack_order)
  ) %>%
  arrange(trait, stack_order_num) %>%
  group_by(trait) %>%
  mutate(
    y_bottom = cumsum(percent_plot) - percent_plot,
    y_mid = y_bottom + percent_plot / 2
  ) %>%
  ungroup() %>%
  inner_join(sig_data_classified, by = c("trait", "source"))

p_final <- ggplot(
  plot_data,
  aes(
    x = trait,
    y = percent_plot,
    fill = source,
    pattern = source
  )
) +
  geom_col_pattern(
    width = 0.88,
    colour = "grey80",
    linewidth = 0.4,
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.25,
    pattern_spacing = 0.035,
    pattern_size = 0.3,
    pattern_shape = 16,
    position = position_stack(reverse = TRUE)
  )

if (nrow(sig_positions) > 0) {
  p_final <- p_final +
    geom_text(
      data = sig_positions,
      aes(
        x = trait, 
        y = y_mid, 
        label = "*", 
        colour = sig_level
      ),
      inherit.aes = FALSE,
      size = 7,
      fontface = "bold"
    ) +
    scale_colour_manual(
      name = NULL,
      values = c("75th" = "blue", "95th" = "red"),
      labels = c("75th Percentile (> Null)", "95th Percentile (> Null)")
    )
}

p_final <- p_final +
  scale_fill_manual(
    name = NULL,
    breaks = c("Both", "Neither", "A. stappersi only", "T. pharyngalis only"),
    values = c(
      "Neither"             = "white",
      "Both"                = "lightgrey",
      "T. pharyngalis only" = "white",
      "A. stappersi only"   = "black"
    ),
    labels = c(
      "Both",
      "Neither",
      "Congolese ancestor: A. stappersi",
      "Nilotic ancestor: T. pharyngalis"
    )
  ) +
  scale_pattern_manual(
    values = c(
      "A. stappersi only"   = "none",
      "T. pharyngalis only" = "pch",
      "Both"                = "none",
      "Neither"             = "none"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    labels = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_discrete(
    limits = traits_ordered,
    labels = function(x) {
      result <- toupper(x)
      result[x == "bb"] <- "BF"
      result
    }
  ) +
  labs(x = NULL, y = "Proportion of top SNPs") +
  guides(
    fill = guide_legend(
      override.aes = list(
        pattern = c("none", "none", "none", "pch"),
        pattern_fill = "black",
        pattern_colour = "black",
        pattern_shape = 16,
        fill = c("lightgrey", "white", "black", "white"),
        colour = "grey60"
      )
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10)
  )

print(p_final)
ggsave("ancestral_enrichment_75th_95th_thresholds_mafcorrected_BG_ALT.svg", plot = p_final, width = 9, height = 4, device = "svg")



