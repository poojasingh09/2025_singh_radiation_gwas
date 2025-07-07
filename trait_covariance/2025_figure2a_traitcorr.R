

### pooja.singh09@gmail.com
## covariance among phenotypes


setwd('/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/trait_covar')


## read in trait info

col <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/gwas_colour_phenotypes_v2.lvrs.NA.txt", header=T, sep="\t", row.names=1)
col <- col[,c(-9)]
head(col)

pcs <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_jawshape/gwas_headshape_PC12_edited.lvrs.txt", header=T, sep="\t", row.names=1)
head(pcs)

dent <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_bcftools_pnye3_dentition_repeat/2025_LVcichlid_traits_GWAS.NA.txt", header=T, sep="\t", row.names=1)
head(dent)

a <- merge(col, dent, by="row.names",all = TRUE)
b <- merge(a, pcs, by.x="Row.names", by.y="row.names",all = TRUE)
out <- b[,c(-1)]
rownames(out) <- b$Row.names

head(out)

# plot

library(pheatmap)
library(viridis)
annotation_row <- data.frame(
  Group = factor(c(rep("Male nuptial colour", 5), rep("Melanic stripe pattern", 3), rep("Ecomorphology", 6)))
)

annotation_colors <- list(
  Group = c("Ecomorphology" = "#4DBBD5", "Melanic stripe pattern" ="#009E73",
            "Male nuptial colour" = "#E64B35")
)

rownames(annotation_row) <- rownames(cor(out))

svg("figure2a_traitcorr.svg", w=8, h=8)
pheatmap(cor(out, method="spearman", use = "complete.obs"), 
         color = plasma(256), annotation_row = annotation_row,
         annotation_colors = annotation_colors,border_color = NA)

dev.off()
write.table(cor(out, method="spearman", use = "complete.obs"), "cor_alltraits_spearmansrho.txt", quote=F)
