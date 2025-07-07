## pooja.singh09@gmail.com
## circos plot to indicate interacting loci
## in this script v3 i have implemented colouring the snp lines based on trait colour
library(circlize)
library(viridis)
library(scales)
library(RColorBrewer)
library(ComplexHeatmap)
library(gridExtra)
library(grid)
library(cowplot)

# Set working directory
setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ld")

# Increase time limits
setTimeLimit(cpu = Inf, elapsed = Inf)

# Load top SNPs
top <- read.table("../fst/allcomparisons_top_gwas_pip01_SNP_0.95fst_tophits_SNPoverlaponly.txt", header=TRUE)

trait_categories <- c(
  "CS" = "eco", "IR"= "eco", "KD"= "eco", "TM"= "eco", 
  "PC1"= "eco", "PC2"= "eco",
  "VB"="pattern", "MLB"="pattern", "DLB"="pattern",
  "FB"="nuptial", "YF"="nuptial", "RC"="nuptial", "MB"="nuptial", "BB"="nuptial"
)

category_colours <- c("eco" = "#4DBBD5", "pattern"="#009E73", "nuptial"="#E64B35")

top$trait_complex <- trait_categories[top$trait]
top$colour <- category_colours[top$trait_complex]

# Load ideogram
custom_ideogram <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/pnye_v3_assembly/pnyev3.fa.fai.ideogram.txt", header=TRUE)

# Load LD SNPs
interacting_snps1 <- read.table("LD_r2_results_allpairs_combined_annotated.txt", header=TRUE)
interacting_snps2 <- interacting_snps1[interacting_snps1$SNP1chr != interacting_snps1$SNP2chr | abs(interacting_snps1$SNP1pos - interacting_snps1$SNP2pos) > 100000, ]
interacting_snps_unique <- interacting_snps2[!duplicated(interacting_snps2), ]

# Color scale
green_palette <- colorRampPalette(c("white", "lightgreen", "darkgreen"))
col_fun <- colorRamp2(
  seq(min(interacting_snps_unique$Rsq), max(interacting_snps_unique$Rsq), length.out = 256), 
  alpha(green_palette(256), 1)
)

# Comparisons
comparisons <- c(
  "entantleter_MG_vs_newdegeni_MG_20kb.windowed.weir.fst_top",
  "ecoprologous_blue_MG_vs_eparopius_MG_20kb.windowed.weir.fst_top",
  "ecinctus_MG_vs_ecoprologous_blue_MG_20kb.windowed.weir.fst_top",
  "yplumbus_vs_ypyrrocephalus_20kb.windowed.weir.fst_top",
  "lscraper_makobe_vs_lyellowchin_makobe_20kb.windowed.weir.fst_top",
  "nomni_makobe_vs_nuniscaper_makobe_20kb.windowed.weir.fst_top",
  "ngigas_makobe_vs_paracyaneus_makobe_20kb.windowed.weir.fst_top",
  "mmbipi_makobe_vs_ppinkanal_makobe_20kb.windowed.weir.fst_top",
  "pnye_kiss_vs_ppun_kiss_20kb.windowed.weir.fst_top",
  "pnye_python_vs_ppun_python_20kb.windowed.weir.fst_top",
  "pnye_ruti_vs_ppun_ruti_20kb.windowed.weir.fst_top",
  "pnye_makobe_vs_ppun_makobe_20kb.windowed.weir.fst_top"
)

# Circos plot function
create_circos_plot <- function(comparison) {
  interacting_snps <- interacting_snps_unique[interacting_snps_unique$comparison == comparison, ]
  interacting_snps$SNP1end <- interacting_snps$SNP1pos
  interacting_snps$SNP2end <- interacting_snps$SNP2pos
  
  # Clean gene names
  interacting_snps$SNP1_gene <- sub(".*Name=([^;]+)(;.*|$)", "\\1", interacting_snps$SNP1_gene)
  interacting_snps$SNP2_gene <- sub(".*Name=([^;]+)(;.*|$)", "\\1", interacting_snps$SNP2_gene)
  interacting_snps$SNP1_gene[interacting_snps$SNP1_gene == ""] <- NA
  interacting_snps$SNP2_gene[interacting_snps$SNP2_gene == ""] <- NA
  
  # Load top SNPs for this comparison
  comparison1 <- sub("_top$", "", comparison)
  top1 <- top[top$comparison == comparison1, ]
  
  # SNP positions for track
  snp_positions <- unique(data.frame(
    chr = c(interacting_snps$SNP1chr, interacting_snps$SNP2chr),
    start = c(interacting_snps$SNP1pos, interacting_snps$SNP2pos),
    end = c(interacting_snps$SNP1end, interacting_snps$SNP2end),
    value = 1
  ))
  
  snp_positions1 <- merge(snp_positions, top1, by.x = c("chr", "start", "end"), by.y = c("seqnames", "start", "end"))
  
  snp_positions2 <- snp_positions1 %>%
    group_by(chr, start, end) %>%
    mutate(colour = ifelse(n() > 1, "grey", colour)) %>%
    distinct()
  
  snp_positions2 <- data.frame(snp_positions2)
  snp_positions2 <- snp_positions2[, c(1:3, 24)]
  snp_positions2$value <- 1
  snp_positions2 <- unique(snp_positions2)
  
  # Title
  plot_title <- sub("_20kb.*", "", comparison)
  
  # Initialize plot
  circos.clear()
  chromosomes <- unique(custom_ideogram$chr)
  n_chr <- length(chromosomes)
  
  circos.par(
    start.degree = 90,
    gap.after = c(rep(2, n_chr - 1), 10),  # Add a final gap after the last sector
    track.margin = c(0.01, 0.01),
    cell.padding = c(0.01, 0.01, 0.01, 0.01),
    canvas.xlim = c(-1.3, 1.3),
    canvas.ylim = c(-1.3, 1.3)
  )
  circos.initialize(factors = custom_ideogram$chr, xlim = cbind(custom_ideogram$start, custom_ideogram$end))
  
  # SNP track
  circos.genomicTrack(snp_positions2, ylim = c(0, 1), panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, col = value$colour, lwd = 3, type = "h")
  }, track.height = 0.1)
  
  # Mean Rsq threshold
  mean_rsq <- 0.01
  interacting_snps <- interacting_snps[order(interacting_snps$Rsq), ]
  
  # Labels + LD lines
  for (i in 1:nrow(interacting_snps)) {
    label1 <- ifelse(
      is.na(interacting_snps$SNP1_gene[i]),
      paste0("(", interacting_snps$SNP1_trait[i], ")"),
      paste0(interacting_snps$SNP1_gene[i], " (", interacting_snps$SNP1_trait[i], ")")
    )
    
    circos.text(
      x = interacting_snps$SNP1pos[i],
      y = 2.5,
      labels = label1,
      sector.index = interacting_snps$SNP1chr[i],
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.6,
      col = "black"
    )
    
    label2 <- ifelse(
      is.na(interacting_snps$SNP2_gene[i]),
      paste0("(", interacting_snps$SNP2_trait[i], ")"),
      paste0(interacting_snps$SNP2_gene[i], " (", interacting_snps$SNP2_trait[i], ")")
    )
    
    circos.text(
      x = interacting_snps$SNP2pos[i],
      y = 2.5,
      labels = label2,
      sector.index = interacting_snps$SNP2chr[i],
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.6,
      col = "black"
    )
    
    # ✅ ADD THE LD LINK:
    if (interacting_snps$Rsq[i] > mean_rsq) {
      circos.link(
        sector.index1 = interacting_snps$SNP1chr[i],
        point1 = interacting_snps$SNP1pos[i],
        sector.index2 = interacting_snps$SNP2chr[i],
        point2 = interacting_snps$SNP2pos[i],
        col = col_fun(interacting_snps$Rsq[i]),
        lwd = 1.5
      )
    }
  }
  
  # Add legend
  heatmap_legend <- Legend(
    col_fun = col_fun,
    title = "Rsq",
    at = seq(min(interacting_snps$Rsq), max(interacting_snps$Rsq), length.out = 5),
    labels = round(seq(min(interacting_snps$Rsq), max(interacting_snps$Rsq), length.out = 5), 2)
  )
  draw(heatmap_legend, x = unit(1, "npc") - unit(1, "cm"), y = unit(1, "npc") - unit(19, "cm"), just = c("right", "top"))
  
  # Add chromosome labels
  for (i in 1:nrow(custom_ideogram)) {
    circos.text(
      x = (custom_ideogram$start[i] + custom_ideogram$end[i]) / 2,
      y = 1.8,
      labels = custom_ideogram$chr[i],
      sector.index = custom_ideogram$chr[i],
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.8,
      col = "black"
    )
  }
  
  # Plot title
  # Plot title
  title(main = plot_title, cex.main = 1)
  
  # ✅ Return the full base R plot (including circos + legend)
  return(recordPlot())
  
}

# Generate all plots
plot_list <- lapply(comparisons, create_circos_plot)




####
pdf("circos_plots_all_comparisons_v3_rsq0.01.pdf", width = 10, height = 10)
for (p in plot_list) {
  replayPlot(p)  # replay full recorded plot
}
dev.off()

## svg
# Specify only the plot indices you want
selected_plots <- c(1,2,3,4,5,6,7,8,9,10,11,12)

for (i in selected_plots) {
  svg_filename <- paste0("circos_plot_", i, ".svg")
  svg(svg_filename, width = 10, height = 10)
  replayPlot(plot_list[[i]])  # replay the recorded plot
  dev.off()
}