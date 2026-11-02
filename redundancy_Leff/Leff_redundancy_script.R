#################### OCT 2026 ################ 

#### GENETIC REDUNDANCY Analysis By Pooja Singh
### C-score and Leff adapted from Yeaman et al 2018 and Soudi et al 2023


####### read data
setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/phylogeny_v2")
dat <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/gwas_colour_phenotypes_v2.lvrs.NA.txt", sep="\t", row.names=1, header=T)
head(dat)

library(dplyr)

# 1. subset to the first 8 colour columns as you intended
dat <- dat[, 1:8]

# 2. Convert binary 0/1 integers into the "A" / "B" character format
dat[dat == 0] <- "A"
dat[dat == 1] <- "B"

# 3. Strip individual suffixes and collapse multiple rows into species consensus
# This mimics your Step 5 phenotype collapse block perfectly
dat_collapsed <- dat %>%
  # Extract individual row IDs into a workable column
  mutate(individual_id = rownames(dat)) %>%
  # Strip individual suffix tags (e.g., "Pundamilia nyererei_5" becomes "Pundamilia nyererei")
  mutate(species_clean = gsub("_[0-9]+$", "", individual_id)) %>%
  # Remove outgroup 106985 if it appears in any string variants
  filter(species_clean != "106985") %>%
  # Group rows belonging to the same biological species together
  group_by(species_clean) %>%
  # Collapse logic: If ANY individual of a species has the trait ("B"), 
  # the species is marked "B". If all are NA, it stays NA. Otherwise "A".
  summarise(across(-individual_id, ~ {
    if (any(.x == "B", na.rm = TRUE)) "B"
    else if (all(is.na(.x))) NA_character_
    else "A"
  }), .groups = 'drop')

# 4. Finalise dat_final structure with unique species row names
dat_final <- as.data.frame(dat_collapsed %>% select(-species_clean))
rownames(dat_final) <- dat_collapsed$species_clean

# 5. Sanity Check: Ensure columns and rows are clean
print("--- DAT_FINAL SUCCESSFULLY BUILT ---")
print(dim(dat_final))
print(head(dat_final, 3))

# Read in alleles for AGRP
af <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ancestor_matching/all_rwas_tophits_pip01.allthegoodstuff.txt", header=T, sep="\t")
gt <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ancestor_matching/all_rwas_tophits_pip01.GT.allsamples.txt", header=T, sep="\t", check.names = FALSE)

library(dplyr)
library(tidyr)

# 1. Standardize column names and create Unified SNP IDs
colnames(af)[1] <- "Locus_ID"
gt$Unified_ID <- paste0(gsub("chr", "", gt[, 1]), "-", gt[, 2])

# 2. Reshape genotype data from wide columns to long format
gt_prepped <- gt %>% 
  select(Unified_ID, everything()) %>%
  select(-2, -3)  # Drops the raw redundant Chr and Pos columns

gt_long <- gt_prepped %>% 
  pivot_longer(cols = -Unified_ID, names_to = "sample_id", values_to = "dosage")

# 3. Map individual sample IDs to species names
# Uses species_df (loaded earlier) to associate each sample with its tree label
gt_long$tree_label <- species_df$species_name[match(gt_long$sample_id, rownames(species_df))]

# Filter out missing records and strip replicate suffix flags (e.g. "_4" to match species names)
gt_long_clean <- gt_long %>% 
  filter(!is.na(tree_label)) %>%
  mutate(species_clean = gsub("_[0-9]+$", "", tree_label))

# 4. Calculate continuous Observed Heterozygosity (Ho) 
# An individual is a heterozygous carrier if dosage is exactly 1
gt_long_clean <- gt_long_clean %>%
  mutate(
    is_het = ifelse(!is.na(dosage) & dosage == 1, 1, 0),
    is_valid = ifelse(!is.na(dosage), 1, 0)
  )

# Aggregate values by biological species groups
species_het_summary <- gt_long_clean %>%
  group_by(species_clean, Unified_ID) %>%
  summarise(
    Valid_Samples = sum(is_valid),
    Het_Samples = sum(is_het),
    Species_Ho = ifelse(Valid_Samples > 0, Het_Samples / Valid_Samples, 0),
    .groups = 'drop'
  )

# 5. Reshape summary table back into wide matrix format
het_matrix_wide <- species_het_summary %>%
  select(species_clean, Unified_ID, Species_Ho) %>%
  pivot_wider(names_from = Unified_ID, values_from = Species_Ho, values_fill = 0)

# Convert to a true numeric R matrix and attach row names
het_matrix_raw <- as.matrix(het_matrix_wide[, -1])
rownames(het_matrix_raw) <- het_matrix_wide$species_clean

# Force any multiple or double spaces to clean single spaces for safety
rownames(het_matrix_raw) <- gsub("\\s+", " ", trimws(rownames(het_matrix_raw)))

# 6. Sanity Verification
print("--- HET_MATRIX_RAW COMPILED SUCCESSFULY ---")
print(dim(het_matrix_raw))
print(het_matrix_raw[1:3, 1:3])

library(ape)
library(dplyr)

# 1. Calculate the raw cophenetic (patristic) distance matrix across all tree branches
tree_dist_matrix <- cophenetic(pruned.tree)

# 2. Clean species names by removing individual replicate number suffixes 
# (e.g., "Pundamilia nyererei_5" collapses to "Pundamilia nyererei")
colnames(tree_dist_matrix) <- gsub("_[0-9]+$", "", colnames(tree_dist_matrix))
rownames(tree_dist_matrix) <- gsub("_[0-9]+$", "", rownames(tree_dist_matrix))

# 3. Reshape wide distance matrix into a long tidy dataframe format
melted_dist <- as.data.frame(as.table(tree_dist_matrix))
colnames(melted_dist) <- c("SpX", "SpY", "Dist")

# 4. Clean formatting and average distances across replicate branches
melted_dist$SpX <- gsub("\\s+", " ", trimws(as.character(melted_dist$SpX)))
melted_dist$SpY <- gsub("\\s+", " ", trimws(as.character(melted_dist$SpY)))

species_dist_summary <- melted_dist %>% 
  group_by(SpX, SpY) %>% 
  summarise(Mean_Phylo_Dist = mean(Dist), .groups = 'drop')

# 5. Sanity Verification
print("--- SPECIES_DIST_SUMMARY DISTANCE MATRIX COMPILED ---")
print(dim(species_dist_summary))
print(head(species_dist_summary, 3))

# ==============================================================================
# FINAL MOLECULAR CANVAS PLOT: AXIS EXTENSION FOR UN-CLIPPED LABELS
# ==============================================================================

cichlid_palette <- c(
  "Flameback" = "#FF0000", "MelanicBody" = "#333333", "MidLateralBand" = "#0000FF",
  "RedChest" = "#E69F00", "BlueBody" = "#56B4E9", "YellowFlank" = "#F0E442",
  "VerticalBars" = "#CC79A7", "DorsalLateralBand" = "#009E73"
)

# Render the final chart with expanded horizontal padding
plt_binary_final <- ggplot(leff_results, aes(x = Phylo_Distance, y = Leff_Est, color = trait)) +
  geom_point(size = 1.4, alpha = 0.18, stroke = 0) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.3) +
  
  # Standard stable text overlays positioned right next to line trajectories
  geom_text(
    data = label_data, 
    aes(x = x + 0.001, y = y, label = label, color = trait),
    hjust = 0, 
    fontface = "bold", 
    size = 3.2, 
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = cichlid_palette) +
  
  # CRITICAL CORRECTION: Forces a 15% empty margin cushion on the right side 
  # of the x-axis so that text labels have plenty of space to display fully
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  
  theme_classic() +
  labs(
    title = "Macroevolutionary Target Space Expansion Across Cichlid Pigmentation Traits",
    x = "Averaged Patristic Distance (Phylogeny Branch Lengths)",
    y = "Effective Loci Pool Size (Leff = (a1 * a2) / ab)",
    color = "Trait Name"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.title = element_text(size = 10, face = "bold"),
    legend.title = element_text(face = "bold", size = 10),
    plot.margin = margin(10, 20, 10, 10) # Clean standardized window border padding
  )

print(plt_binary_final)
#save plot if you like

