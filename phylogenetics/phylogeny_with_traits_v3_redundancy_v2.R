

## pooja.singh09@gmail.com
## may 2025
## plotting traits and alleles of top snps on phylogeny to see impact of pseudorpelication

setwd("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/phylogeny_v2/")
library(diversitree)
library(ape)
library(ggtree)
library(ggnewscale)
library(geiger)
library(treeio)
library(ggplot2)
library(phylobase)
library(adephylo)
library(wesanderson)
library(dplyr)
library(gridExtra)
library(grid)

tree <- ape::read.tree("all.filt.vcf.var.1kbthinned.gtr.dna.treefile")
plot(tree,show.node=TRUE)
head(listTips(tree))

#cols = list(YellowFlank = c("grey", "black"), Flameback= c("grey", "black"), RedChest= c("grey", "black"), BlueBody= c("grey", "black"), MelanicBody= c("grey", "black"), VerticalBars= c("grey", "black"), MidLateralBand= c("grey", "black"), DorsalLateralBand= c("grey", "black"))

#dat <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_colour/gwas_colour_v1.NA_binary.lvrs.matchedphylogeny.txt", row.names=1, header=T)
dat <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/gwas_colour_phenotypes_v2.lvrs.NA.txt", sep="\t", row.names=1, header=T)
head(dat)
dat <- dat[,1:8]


species <- c(rownames(dat), "71001", "131282") ## "71001" is T pharyngalis, is Aburtoni 131282
species1 <- species != "106985" # remove pseudocrenilabrus multicolor victoriae 106985
species <- species1

dat[dat == 0] <- "A"
dat[dat == 1] <- "B"

dat1 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_jawshape/gwas_headshape_PC12_edited.lvrs.txt", row.names=1, header=T)
head(dat1)


pruned.tree <-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
plot(pruned.tree)

pruned.tree <- root(pruned.tree, "131282", resolve.root = TRUE)
is.rooted(pruned.tree) 


normal <- ggtree(pruned.tree)
normal_p1 <- gheatmap(normal, dat, offset=.8, width=.2,colnames_angle=95, colnames_offset_y = .25)  + scale_fill_viridis_d(option="D", name="discrete\nvalue")
normal_p1 + geom_tiplab(size=3)

#pdf("normal_phylogeny_withtraits_v3a.pdf", h=50, w=10)
#normal_p1 + geom_tiplab(size=0) +geom_nodelab()
normal_p1 + geom_tiplab(size=3)
#dev.off()

circ <- ggtree(pruned.tree, layout="circular", branch.length='none')
circ_p1 <- gheatmap(circ, dat1, offset=.8, width=.07,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_viridis_c(option="A", name="continuous\nvalue",na.value="white") 

circ_p2 <- circ_p1 + new_scale_fill()
p2 <- gheatmap(circ_p2, dat, offset=3, width=.2,
         colnames_angle=90, colnames_offset_y = 0.25, colnames = F) + scale_fill_manual(values = c("skyblue","orange"),  na.value = "white")
circ_p1 + geom_tiplab(size=2)

p2 

#ggsave(file="circular_phylogeny_withtraits_v2.png", plot=p2, width=10, height=8)
#ggsave("circular_phylogeny_withtraits_v2.svg", device = svglite::svglite,  width=10, height=8)


############### lets add species names to the tips first

# Create a new data frame with species names
tree <- ape::read.tree("all.filt.vcf.var.1kbthinned.gtr.dna.treefile")
plot(tree,show.node=TRUE)
head(listTips(tree))
length(tree$tip.label)
#write.table(tree$tip.label, "all.chr.filt.2kbThinned.raxmlinput.min4.phy.fa.treefile.tre.tips", quote=F, row.names=F)


## read in species names
species_df <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/2023_sampleselection_gwas_singh_manuscriptsubset_OS_PS_removedsamples.txt", sep="\t", header=T, row.names=1)
colnames(species_df) <- c("species_name", "short") 
head(species_df)
dim(species_df)

# Function to add suffixes to duplicates so they are unique
add_suffix_to_duplicates <- function(names) {
  duplicated_names <- names[duplicated(names)]
  for (name in duplicated_names) {
    indices <- which(names == name)
    names[indices] <- paste0(name, "_", seq_along(indices))
  }
  return(names)
}

# Apply the function to the data frame
species_df <- species_df %>%
  mutate(species_name = add_suffix_to_duplicates(species_name))

head(species_df)

species_df1 <- species_df[!(species_df$species_name %in% c("Lipochromis melanopterus_1", 
                                                          "Astatotilapia nubila swamp_6", 
                                                          "Pundamilia yellow azurea_3", 
                                                          "Ptyochromis (Macropleurodus super-lineage) striped rock sheller_3",
                                                          "Pundamilia azurea_2")), , drop = FALSE]


species_df <- species_df1 
head(species_df1)

## now to short names
species_df2 <- species_df1 %>%
  mutate(short= add_suffix_to_duplicates(short))
head(species_df2)
species_df <- species_df2

## subset phylogeny for only tips in my data
# Get the tips to keep

#pruned.tree <- drop.tip(tree,tree$tip.label[-match(rownames(species_df), tree$tip.label)])
#plot(pruned.tree)


#pruned.tree <- root(pruned.tree, "106985", resolve.root = TRUE)
#is.rooted(pruned.tree) 


# Get tips to keep: those that exist in both tree and species_df
tips_to_keep <- intersect(tree$tip.label, rownames(species_df))
# Drop tips not in the intersection
pruned.tree <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))
pruned.tree <- root(pruned.tree, "106985", resolve.root = TRUE) #131282 burtoni
is.rooted(pruned.tree) 


# Match species names to pruned tree tips and reset tiplabels
species_names <- species_df$species_name[match(pruned.tree$tip.label, rownames(species_df))]
pruned.tree$tip.label <- species_names


######## read in traits
dat <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/gwas_colour_phenotypes_v2.lvrs.NA.txt", sep="\t", row.names=1, header=T)
head(dat)
dat_p <- dat[,1:8]
rownames(dat_p) <- rownames(dat)
dat <- dat_p
head(dat)

#species <- c(rownames(dat), "71001") ## "71001" is T pharyngalis
dat[dat == 0] <- "A"
dat[dat == 1] <- "B"

dat0 <- merge(dat, species_df, by="row.names", all.x=TRUE)
dat01 <- dat0[!is.na(dat0$species_name), ]
dat <- dat01[,c(2:9)]
rownames(dat) <- dat01$species_name
head(dat)

dat1 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_jawshape/gwas_headshape_PC12_edited.lvrs.txt", row.names=1, header=T)
dat0 <- merge(dat1, species_df, by="row.names", all.x=TRUE) ## there are some species missing here, come back to this ok?
dat0 <- dat0[!is.na(dat0$species_name), ]
dat1 <- dat0[,c(2:3)]
rownames(dat1) <- dat0$species_name

head(dat1) # this has 2 columns of continuous data PC1 and PC2
head(dat) # this is binary data for 8 traits


# Plot the tree with heatmaps without tip labels

# Create a data frame for tip labels
# Create a data frame for tip labels with the node column
#tip_labels_df <- data.frame(node = 1:length(pruned.tree$tip.label), tip.label = pruned.tree$tip.label, species_name = species_names)


circ <- ggtree(pruned.tree, layout="circular", branch.length='none')
circ_p1 <- gheatmap(circ, dat1, offset=.8, width=.07, colnames_angle=95, colnames_offset_y=.1, colnames=FALSE) +
  scale_fill_viridis_c(option="A", name="continuous\nvalue", na.value="white")

circ_p2 <- circ_p1 + new_scale_fill()
p2 <- gheatmap(circ_p2, dat, offset=3, width=.2, colnames_angle=90, colnames_offset_y=0.25, colnames=FALSE) +
  scale_fill_manual(values=c("skyblue", "orange"), na.value="white")

# Add tip labels using the species names data frame
p3 <- p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 9)

# Display the plot with tip labels
print(p3)


#geom_tiplab(aes(label=species_names), fontface="italic", size=2)
#ggtree(pruned.tree, layout="circular", branch.length='none') + geom_tiplab(fontface="italic", size=2)



####################### now let's write a loop that loops over all traits and plots their top 20 snps

library(ggtree)
library(ggplot2)
library(gridExtra)
library(grid) # Load the grid package

traits <- colnames(dat)
traits_short <- c("yf", "fb", "rc", "bb", "mb", "vb", "mlb", "dlb")


# Open a PDF device to save the plots
pdf("all_plots_colour_20.pdf", width = 15, height = 15)

counter <- 0

for (i in traits) {
  counter <- counter + 1
  dat_mlb <- dat[, i, drop = FALSE]
  circ <- ggtree(pruned.tree, layout="circular", branch.length='none')
  p2 <- gheatmap(circ, dat_mlb, offset=0, width=.1, colnames_angle=90, colnames_offset_y=0.25, colnames=FALSE) +
    scale_fill_manual(values=c("grey90", "red"), na.value="white", labels=c("absent", "present")) +
    labs(fill = traits[counter]) + # Use the actual trait being plotted
    theme(legend.position = "top", legend.justification = "right") # Move legend to top right
  print(p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 4))
  
  # Read in alleles for AGRP
  af <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ancestor_matching/all_rwas_tophits_pip01.allthegoodstuff.txt", header=T, sep="\t")
  gt <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ancestor_matching/all_rwas_tophits_pip01.GT.allsamples.txt", header=T, sep="\t", check.names = FALSE)
  
  focus <- traits_short[counter]
  print(focus)
  
  topsnps <- head(af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ], n=30)
  #topsnps <- af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ]
  print(topsnps)
  
  for (j in 1:nrow(topsnps)) {
    snp1 <- topsnps[j, ]
    agrp <- gt[gt$CHROM %in% snp1$chr & gt$POS == snp1$pos, ]
    agrp1 <- data.frame(t(agrp))
    agrp2 <- agrp1[3:nrow(agrp1), , drop = FALSE]
    colnames(agrp2) <- c("agrp")
    agrp3 <- merge(agrp2, species_df, by="row.names", all.x=TRUE)
    agrp4 <- agrp3[!is.na(agrp3$species_name), ]
    agrp5 <- data.frame(agrp4$agrp)
    rownames(agrp5) <- agrp4$species_name
    colnames(agrp5) <- "agrp"
    
    agrp5[agrp5 == 0] <- "A"
    agrp5[agrp5 == 1] <- "B"
    agrp5[agrp5 == 2] <- "C"
    
    # Create plot p3 for each SNP
    circ_p1 <- gheatmap(circ, dat_mlb, offset = 0.8, width = .1, colnames_angle = 90, colnames_offset_y = 0.25, colnames = FALSE) +
      scale_fill_manual(values = c("grey90", "red"), na.value = "white", labels = c("absent", "present")) + 
      labs(fill = traits[counter]) + # Use the actual trait being plotted
      theme(legend.position = "top", legend.justification = "right") # Move legend to top right
    
    circ_p2 <- circ_p1 + new_scale_fill()
    
    agrp5_matrix <- as.matrix(agrp5)
    
    # Extract the string between "Name=" and ";" from snp1$gene
    gene_name <- ifelse(is.na(snp1$gene), NA, sub(".*Name=([^;]*);.*", "\\1", snp1$gene))
    
    # Determine the legend title based on whether gene_name is NA
    legend_title <- ifelse(is.na(gene_name), 
                           paste("SNP:", snp1$chr, snp1$pos), 
                           paste("SNP:", snp1$chr, snp1$pos, gene_name))
    
    p2 <- gheatmap(circ_p2, agrp5_matrix, offset = 3, width = 0.1, colnames_angle = 90, colnames_offset_y = 0.25, colnames = FALSE) +
      scale_fill_manual(values = c("grey90", "green", "darkgreen"), na.value = "white", labels = c("00", "01", "11")) + 
      labs(fill = legend_title) + # Use the determined legend title
      theme(legend.position = "top", legend.justification = "right") # Move legend to top right
    
    p3 <- p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 9)
    
    # Print the plot for each SNP
    print(p3)
  }
}

# Close the PDF device
dev.off()



######## for traits for pc1,pc2

library(ggtree)
library(ggplot2)
library(gridExtra)
library(grid) # Load the grid package

traits <- colnames(dat1)
traits_short <- c("pc1", "pc2")

# Open a PDF device to save the plots
pdf("all_plots_shape_20.pdf", width = 15, height = 15)

counter <- 0

for (i in traits) {
  counter <- counter + 1
  dat_mlb <- dat1[, i, drop = FALSE]
  circ <- ggtree(pruned.tree, layout="circular", branch.length='none')
  #circ1 <- circ1 + new_scale_fill()
  p2 <- gheatmap(circ, dat_mlb, offset=0, width=0.1, 
                 colnames_angle=90, colnames_offset_y=0.25, colnames=FALSE) +
    scale_fill_viridis_c(option="A", name=traits[counter], na.value="white",direction = -1) +  # Use the trait name directly
    theme(legend.position = "top", 
          legend.justification = "right",
          legend.title.align = 0.5,  # Center the legend title
          legend.spacing.x = unit(0.5, "cm"),  # Add space between legend and plot
          legend.margin = margin(10, 10, 10, 10))  # Add space around the legend
  # Read in alleles for AGRP
  
  af <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ancestor_matching/all_rwas_tophits_pip01.allthegoodstuff.txt", header=T, sep="\t")
  gt <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ancestor_matching/all_rwas_tophits_pip01.GT.allsamples.txt", header=T, sep="\t", check.names = FALSE)
  
  focus <- traits_short[counter]
  print(focus)
  
  topsnps <- head(af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ], n=30)
  #topsnps <- af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ]
  print(topsnps)
  
  for (j in 1:nrow(topsnps)) {
    snp1 <- topsnps[j, ]
    agrp <- gt[gt$CHROM %in% snp1$chr & gt$POS == snp1$pos, ]
    agrp1 <- data.frame(t(agrp))
    agrp2 <- agrp1[3:nrow(agrp1), , drop = FALSE]
    colnames(agrp2) <- c("agrp")
    agrp3 <- merge(agrp2, species_df, by="row.names", all.x=TRUE)
    agrp4 <- agrp3[!is.na(agrp3$species_name), ]
    agrp5 <- data.frame(agrp4$agrp)
    rownames(agrp5) <- agrp4$species_name
    colnames(agrp5) <- "agrp"
    
    #agrp5[agrp5 == 0] <- "A"
    #agrp5[agrp5 == 1] <- "B"
    #agrp5[agrp5 == 2] <- "C"
    
    # Create plot p3 for each SNP
    circ_p1 <- gheatmap(circ, dat_mlb, offset=0, width=0.1, 
                        colnames_angle=90, colnames_offset_y=0.25, colnames=FALSE) +
      scale_fill_viridis_c(option="A", name=traits[counter], na.value="white", direction=-1) +  # Use the trait name directly
      theme(legend.position = "top", 
            legend.justification = "right",
            legend.title.align = 0.5,  # Center the legend title
            legend.spacing.x = unit(0.5, "cm"),  # Add space between legend and plot
            legend.margin = margin(10, 10, 10, 10))  # Add space around the legend
    
    
    circ_p2 <- circ_p1 + new_scale_fill()
    
    agrp5_matrix <- as.matrix(agrp5)
    
    # Extract the string between "Name=" and ";" from snp1$gene
    gene_name <- ifelse(is.na(snp1$gene), NA, sub(".*Name=([^;]*);.*", "\\1", snp1$gene))
    
    # Determine the legend title based on whether gene_name is NA
    legend_title <- ifelse(is.na(gene_name), 
                           paste("SNP:", snp1$chr, snp1$pos), 
                           paste("SNP:", snp1$chr, snp1$pos, gene_name))
    
    p2 <- gheatmap(circ_p2, agrp5_matrix, offset = 3, width = 0.1, colnames_angle = 90, colnames_offset_y = 0.25, colnames = FALSE) +
      scale_fill_manual(values = c("grey90", "green", "darkgreen"), na.value = "white", labels = c("00", "01", "11", "NA")) + 
      labs(fill = legend_title) + # Use the determined legend title
      theme(legend.position = "top", legend.justification = "right") # Move legend to top right
    
    p3 <- p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 9)
    
    # Print the plot for each SNP
    print(p3)
  }
}

# Close the PDF device
dev.off()


####################### now let's write a loop that loops over dentition traits

library(ggtree)
library(ggplot2)
library(gridExtra)
library(grid) # Load the grid package


dat2 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_bcftools_pnye3_dentition_repeat/2025_LVcichlid_traits_GWAS.NA.txt", sep="\t", row.names=1, header=T)
head(dat2)

dat0 <- merge(dat2, species_df, by="row.names", all.x=TRUE) ## there are some species missing here, come back to this ok?
dat0 <- dat0[!is.na(dat0$species_name), ]
head(dat0)
dat2 <- dat0[,c(2:5)]
rownames(dat2) <- dat0$species_name

species <- rownames(dat2) 

traits <- colnames(dat2)
traits_short <- c("cs", "ir", "tm", "kd")

# Open a PDF device to save the plots
pdf("all_plots_dentition_repeat_20.pdf", width = 15, height = 15)

counter <- 0

for (i in traits) {
  counter <- counter + 1
  dat_mlb <- dat2[, i, drop = FALSE]
  circ <- ggtree(pruned.tree, layout="circular", branch.length='none')
  #circ1 <- circ1 + new_scale_fill()
  
  # Plotting with adjusted offset and width for larger gaps
  p2 <- gheatmap(circ, dat_mlb, offset=0, width=0.1, 
                       colnames_angle=90, colnames_offset_y=0.25, colnames=FALSE) +
    scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white", name=traits[counter]) +  # Use the trait name directly
    theme(legend.position = "top", 
          legend.justification = "right",
          legend.title.align = 0.5,  # Center the legend title
          legend.spacing.x = unit(0.5, "cm"),  # Add space between legend and plot
          legend.margin = margin(10, 10, 10, 10))  # Add space around the legend
  
  #print(p2)
  
  
  print(p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 4))
  # Read in alleles for AGRP
  
  af <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ancestor_matching/top_rwas_all_99.9_and_pip01.allthegoodstuff.txt", header=T, sep="\t", check.names = FALSE)
  gt <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ancestor_matching/top_rwas_all_99.9_and_pip01.GT.txt", header=T, sep="\t", check.names = FALSE)
  
  
  focus <- traits_short[counter]
  print(focus)
  
  topsnps <- head(af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ], n=20)
  #topsnps <- af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ]
  print(topsnps)
  
  for (j in 1:nrow(topsnps)) {
    snp1 <- topsnps[j, ]
    agrp <- gt[gt$CHROM %in% snp1$chr & gt$POS == snp1$pos, ]
    agrp1 <- data.frame(t(agrp))
    agrp2 <- agrp1[-c(1, 2), , drop = FALSE]
    colnames(agrp2) <- c("agrp")
    agrp3 <- merge(agrp2, species_df, by="row.names", all.x=TRUE)
    agrp4 <- agrp3[!is.na(agrp3$species_name), ]
    agrp5 <- data.frame(agrp4$agrp)
    rownames(agrp5) <- agrp4$species_name
    colnames(agrp5) <- "agrp"
    
    #agrp5[agrp5 == 0] <- "A"
    #agrp5[agrp5 == 1] <- "B"
    #agrp5[agrp5 == 2] <- "C"
    
    # Create plot p3 for each SNP
    circ_p1 <- gheatmap(circ, dat_mlb, offset=0, width=0.1, 
                        colnames_angle=90, colnames_offset_y=0.25, colnames=FALSE) +
      scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white", name=traits[counter]) + # Use the trait name directly
      theme(legend.position = "top", 
            legend.justification = "right",
            legend.title.align = 0.5,  # Center the legend title
            legend.spacing.x = unit(0.5, "cm"),  # Add space between legend and plot
            legend.margin = margin(10, 10, 10, 10))  # Add space around the legend
    
    
    circ_p2 <- circ_p1 + new_scale_fill()
    
    agrp5_matrix <- as.matrix(agrp5)
    
    # Extract the string between "Name=" and ";" from snp1$gene
    gene_name <- ifelse(is.na(snp1$gene), NA, sub(".*Name=([^;]*);.*", "\\1", snp1$gene))
    
    # Determine the legend title based on whether gene_name is NA
    legend_title <- ifelse(is.na(gene_name), 
                           paste("SNP:", snp1$chr, snp1$pos), 
                           paste("SNP:", snp1$chr, snp1$pos, gene_name))
    
    p2 <- gheatmap(circ_p2, agrp5_matrix, offset = 3, width = 0.1, colnames_angle = 90, colnames_offset_y = 0.25, colnames = FALSE) +
      scale_fill_manual(values = c("grey90", "green", "darkgreen"), na.value = "white", labels = c("00", "01", "11")) + 
      labs(fill = legend_title) + # Use the determined legend title
      theme(legend.position = "top", legend.justification = "right") # Move legend to top right
    
    p3 <- p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 9)
    
    # Print the plot for each SNP
    print(p3)
  }
}

# Close the PDF device
dev.off()





################################## figure 3 ####################### MLB


library(ggtree)
library(ggplot2)
library(gridExtra)
library(grid) # Load the grid package
library(dplyr)

#traits <- colnames(dat)
#traits_short <- c("yf", "fb", "rc", "bb", "mb", "vb", "mlb", "dlb")

traits <- c("MidLateralBand")
traits_short <- c("mlb")

# Open a PDF device to save the plots


counter <- 0

for (i in traits) {
  counter <- counter + 1
  #dat_mlb <- dat[, i, drop = FALSE]
  dat_mlb <- dat[, as.character(i), drop = FALSE]
  circ <- ggtree(pruned.tree, layout="circular", branch.length='none', root.position = 10)

  # Identify descendants of node 491
  # view node labels library(ggtree)
  #circ_p2 + geom_tiplab(fontface = "italic", size = 2, inherit.aes = FALSE, offset = 7.1) + geom_text2(aes(label = node), hjust = -0.3, size = 2, color = "blue")
  
  desc_nodes <- offspring(pruned.tree, 378)
  
  # Color the branches
  circ <- circ + geom_tree(aes(color = ifelse(node %in% desc_nodes, "red", "black"))) +
    scale_color_identity()  # Ensures the colors are used as specified
  
  p2 <- gheatmap(circ, dat_mlb, offset=0, width=.2, colnames_angle=90, colnames_offset_y=0.25, colnames=F) +
    scale_fill_manual(values=c("grey90", "#009E73"), na.value="white", labels=c("absent", "present")) +
    labs(fill = traits[counter]) + # Use the actual trait being plotted
    theme(legend.position = "top", legend.justification = "right") # Move legend to top right
  #print(p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 4))
  p2
  # Read in alleles for AGRP
  af <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ancestor_matching/all_rwas_tophits_pip01.allthegoodstuff.txt", header=T, sep="\t")
  gt <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ancestor_matching/all_rwas_tophits_pip01.GT.allsamples.txt", header=T, sep="\t", check.names = FALSE)
  
  focus <- traits_short[counter]
  print(focus)
  
  topsnps <- head(af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ], n=20)
  #topsnps <- af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ]
  topsnps1 <- topsnps[c(1,2,3,4,9,17,18),]
  dim(topsnps1)
  
}
    head(topsnps1)  
    agrp <- gt %>%
      semi_join(topsnps1, by = c("CHROM" = "chr", "POS" = "pos")) %>%
      arrange(match(paste(CHROM, POS), paste(topsnps1$chr, topsnps1$pos)))
    agrp1 <- data.frame(t(agrp), stringsAsFactors = FALSE)
    colnames(agrp1) <- paste(agrp1[1, ], agrp1[2, ], sep = "_")
    agrp1 <- agrp1[-c(1,2), , drop = FALSE]  # Drop first two rows but keep structure
    agrp3 <- merge(agrp1, species_df, by="row.names", all.x=TRUE)
    agrp4 <- agrp3[!is.na(agrp3$species_name), ]
    rownames(agrp4) <- agrp4$species_name
    agrp5 <- agrp4[,c(2:9)]
    agrp5[] <- lapply(agrp5, function(x) as.numeric(as.character(x)))
    agrp5[agrp5 == 0] <- "A"
    agrp5[agrp5 == 1] <- "B"
    agrp5[agrp5 == 2] <- "C"
    head(agrp5)
    # Create plot p3 for each SNP
    
    agrp5_matrix <- as.matrix(agrp5)
    circ_p1 <- gheatmap(circ, dat_mlb, offset = 0.5, width = .06,
                        colnames_angle = 95, colnames_offset_y = 0.1, colnames = F) +
      scale_fill_manual(values = c("grey90", "#009E73"), na.value = "white", labels = c("absent", "present"))
    circ_p1 
    circ_p1 <- circ_p1 + new_scale_fill()
    circ_p2 <- gheatmap(circ_p1, agrp5_matrix, offset=1.5,  width=.24,
                        colnames_angle=95, colnames_offset_y = .1, colnames = F) + 
                        scale_fill_manual(values = c("grey90", "lightblue", "blue"), na.value = "white", labels = c("00", "01", "11"))
    circ_p2
    

    p3 <- circ_p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 7.1)
    print(p3)


# save figs with tips and legend
svg("figure3_mlb_redundancy_tips_legend_v3.svg",h=8,w=8)
p3
dev.off()

# save figs without tips and legend
svg("figure3_mlb_redundancy_v3.svg",h=8,w=8)
circ_p2 + theme(legend.position = "none")
dev.off()




################################## figure 3 ####################### KD


library(ggtree)
library(ggplot2)
library(gridExtra)
library(grid) # Load the grid package
library(dplyr)

traits <- c("kd")
traits_short <- c("kd")


# Open a PDF device to save the plots


counter <- 0

for (i in traits) {
  counter <- counter + 1
  #dat_mlb <- dat2[, i, drop = FALSE]
  dat_mlb <- dat2[, as.character(i), drop = FALSE]
  circ <- ggtree(pruned.tree, layout="circular", branch.length='none', root.position = 10)
  # Plotting with adjusted offset and width for larger gaps
  p2 <- gheatmap(circ, dat_mlb, offset=0, width=0.1, 
                 colnames_angle=90, colnames_offset_y=0.25, colnames=FALSE) +
    scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white", name=traits[counter]) +  # Use the trait name directly
    theme(legend.position = "top", 
          legend.justification = "right",
          legend.title.align = 0.5,  # Center the legend title
          legend.spacing.x = unit(0.5, "cm"),  # Add space between legend and plot
          legend.margin = margin(10, 10, 10, 10))  # Add space around the legend
  
  #print(p2)
  
  
  print(p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 4))
  # Read in alleles for AGRP
  af <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ancestor_matching/top_rwas_all_99.9_and_pip01.allthegoodstuff.txt", header=T, sep="\t")
  gt <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ancestor_matching/top_rwas_all_99.9_and_pip01.GT.txt", header=T, sep="\t", check.names = FALSE)
  
  focus <- traits_short[counter]
  print(focus)
  
  topsnps <- head(af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ], n=5)
  #topsnps <- af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ]
  topsnps1 <- topsnps
  dim(topsnps1)
  
}
head(topsnps1)  
agrp <- gt %>%
  semi_join(topsnps1, by = c("CHROM" = "chr", "POS" = "pos")) %>%
  arrange(match(paste(CHROM, POS), paste(topsnps1$chr, topsnps1$pos)))
agrp1 <- data.frame(t(agrp), stringsAsFactors = FALSE)
colnames(agrp1) <- paste(agrp1[1, ], agrp1[2, ], sep = "_")
agrp1 <- agrp1[-c(1,2), , drop = FALSE]  # Drop first two rows but keep structure
agrp3 <- merge(agrp1, species_df, by="row.names", all.x=TRUE)
agrp4 <- agrp3[!is.na(agrp3$species_name), ]
rownames(agrp4) <- agrp4$species_name
agrp5 <- agrp4[,c(2:7)]
agrp5[] <- lapply(agrp5, function(x) as.numeric(as.character(x)))
agrp5[agrp5 == 0] <- "A"
agrp5[agrp5 == 1] <- "B"
agrp5[agrp5 == 2] <- "C"
head(agrp5)
# Create plot p3 for each SNP

agrp5_matrix <- as.matrix(agrp5)
circ_p1 <- gheatmap(circ, dat_mlb, offset = 0.5, width = .06,
                    colnames_angle = 95, colnames_offset_y = 0.1, colnames = F) +
  scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p1 
circ_p1 <- circ_p1 + new_scale_fill()
circ_p2 <- gheatmap(circ_p1, agrp5_matrix, offset=1.5,  width=.24,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + 
  scale_fill_manual(values = c("grey90", "lightblue", "blue"), na.value = "white", labels = c("00", "01", "11"))
circ_p2


p3 <- circ_p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 7.1)
print(p3)


# save figs with tips and legend
svg("figure3_kd_redundancy_tips_legend_v2.svg",h=8,w=8)
p3
dev.off()

# save figs without tips and legend
svg("figure3_kd_redundancy_v2.svg",h=8,w=8)
circ_p2 + theme(legend.position = "none")
dev.off()


###################### if u want to use short labels####################

## but if you want to use short labels
short <- species_df$short[match(pruned.tree$tip.label, rownames(species_df))]
pruned.tree$tip.label <- short


#### testing testing
######## read in traits
dat <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/gwas_colour_phenotypes_v2.lvrs.NA.txt", sep="\t", row.names=1, header=T)
head(dat)
dat_p <- dat[,1:8]
rownames(dat_p) <- rownames(dat)
dat <- dat_p
head(dat)

#species <- c(rownames(dat), "71001") ## "71001" is T pharyngalis
dat[dat == 0] <- "A"
dat[dat == 1] <- "B"

dat0 <- merge(dat, species_df, by="row.names", all.x=TRUE)
dat01 <- dat0[!is.na(dat0$short), ]
dat <- dat01[,c(2:9)]
rownames(dat) <- dat01$short
head(dat)

dat1 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_jawshape/gwas_headshape_PC12_edited.lvrs.txt", row.names=1, header=T)
dat0 <- merge(dat1, species_df, by="row.names", all.x=TRUE) ## there are some species missing here, come back to this ok?
dat0 <- dat0[!is.na(dat0$short), ]
dat1 <- dat0[,c(2:3)]
rownames(dat1) <- dat0$short

head(dat1) # this has 2 columns of continuous data PC1 and PC2
head(dat) # this is binary data for 8 traits



dat2 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_bcftools_pnye3_dentition_repeat/2025_LVcichlid_traits_GWAS.NA.txt", sep="\t", row.names=1, header=T)
head(dat2)

dat0 <- merge(dat2, species_df, by="row.names", all.x=TRUE) ## there are some species missing here, come back to this ok?
dat0 <- dat0[!is.na(dat0$short), ]
head(dat0)
dat2 <- dat0[,c(2:5)]
rownames(dat2) <- dat0$short
species <- rownames(dat2) 
head(dat2)

# Plot the tree with heatmaps without tip labels

# Create a data frame for tip labels
# Create a data frame for tip labels with the node column
#tip_labels_df <- data.frame(node = 1:length(pruned.tree$tip.label), tip.label = pruned.tree$tip.label, species_name = species_names)


circ <- ggtree(pruned.tree, layout="circular", branch.length='none')
circ_p1 <- gheatmap(circ, dat1, offset=.8, width=.07, colnames_angle=95, colnames_offset_y=.1, colnames=FALSE) +
  scale_fill_viridis_c(option="A", name="continuous\nvalue", na.value="white")

circ_p2 <- circ_p1 + new_scale_fill()
p2 <- gheatmap(circ_p2, dat, offset=3, width=.2, colnames_angle=90, colnames_offset_y=0.25, colnames=FALSE) +
  scale_fill_manual(values=c("skyblue", "orange"), na.value="white")

# Add tip labels using the species names data frame
p3 <- p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 9)

# Display the plot with tip labels
print(p3)


################################## figure 3 ####################### MLB with short labs ## may 2025




library(ggtree)
library(ggplot2)
library(gridExtra)
library(grid) # Load the grid package
library(dplyr)

#traits <- colnames(dat)
#traits_short <- c("yf", "fb", "rc", "bb", "mb", "vb", "mlb", "dlb")

traits <- c("MidLateralBand")
traits_short <- c("mlb")

# Open a PDF device to save the plots


counter <- 0

for (i in traits) {
  counter <- counter + 1
  #dat_mlb <- dat[, i, drop = FALSE]
  dat_mlb <- dat[, as.character(i), drop = FALSE]
  circ <- ggtree(pruned.tree, layout="circular", branch.length='none', root.position = 10)
  
  # Identify descendants of node 491
  desc_nodes <- offspring(pruned.tree, 491)
  
  # Color the branches
  circ <- circ + geom_tree(aes(color = ifelse(node %in% desc_nodes, "red", "black"))) +
    scale_color_identity()  # Ensures the colors are used as specified
  
  p2 <- gheatmap(circ, dat_mlb, offset=0, width=.2, colnames_angle=90, colnames_offset_y=0.25, colnames=F) +
    scale_fill_manual(values=c("grey90", "#009E73"), na.value="white", labels=c("absent", "present")) +
    labs(fill = traits[counter]) + # Use the actual trait being plotted
    theme(legend.position = "top", legend.justification = "right") # Move legend to top right
  #print(p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 4))
  p2
  # Read in alleles for AGRP
  af <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ancestor_matching/all_rwas_tophits_pip01.allthegoodstuff.txt", header=T, sep="\t")
  gt <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_dentition/ancestor_matching/all_rwas_tophits_pip01.GT.allsamples.txt", header=T, sep="\t", check.names = FALSE)
  
  focus <- traits_short[counter]
  print(focus)
  
  topsnps <- head(af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ], n=20)
  #topsnps <- af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ]
  topsnps1 <- topsnps[c(1,2,3,4,9,17,18),]
  dim(topsnps1)
  
}
head(topsnps1)  
agrp <- gt %>%
  semi_join(topsnps1, by = c("CHROM" = "chr", "POS" = "pos")) %>%
  arrange(match(paste(CHROM, POS), paste(topsnps1$chr, topsnps1$pos)))
agrp1 <- data.frame(t(agrp), stringsAsFactors = FALSE)
colnames(agrp1) <- paste(agrp1[1, ], agrp1[2, ], sep = "_")
agrp1 <- agrp1[-c(1,2), , drop = FALSE]  # Drop first two rows but keep structure
agrp3 <- merge(agrp1, species_df, by="row.names", all.x=TRUE)
agrp4 <- agrp3[!is.na(agrp3$short), ]
rownames(agrp4) <- agrp4$short
agrp5 <- agrp4[,c(2:9)]
agrp5[] <- lapply(agrp5, function(x) as.numeric(as.character(x)))
agrp5[agrp5 == 0] <- "A"
agrp5[agrp5 == 1] <- "B"
agrp5[agrp5 == 2] <- "C"
head(agrp5)
# Create plot p3 for each SNP

agrp5_matrix <- as.matrix(agrp5)
circ_p1 <- gheatmap(circ, dat_mlb, offset = 0.5, width = .06,
                    colnames_angle = 95, colnames_offset_y = 0.1, colnames = F) +
  scale_fill_manual(values = c("grey90", "#009E73"), na.value = "white", labels = c("absent", "present"))
circ_p1 
circ_p1 <- circ_p1 + new_scale_fill()
circ_p2 <- gheatmap(circ_p1, agrp5_matrix, offset=1.5,  width=.24,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + 
  scale_fill_manual(values = c("grey90", "lightblue", "blue"), na.value = "white", labels = c("00", "01", "11"))
circ_p2


p3 <- circ_p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 7.1)
print(p3)


# save figs with tips and legend
svg("figure3_mlb_redundancy_tips_legend_shortnames.svg",h=8,w=8)
p3
dev.off()


################################## figure 3 ####################### KD


library(ggtree)
library(ggplot2)
library(gridExtra)
library(grid) # Load the grid package
library(dplyr)

traits <- c("kd")
traits_short <- c("kd")


# Open a PDF device to save the plots


counter <- 0

for (i in traits) {
  counter <- counter + 1
  #dat_mlb <- dat2[, i, drop = FALSE]
  dat_mlb <- dat2[, as.character(i), drop = FALSE]
  circ <- ggtree(pruned.tree, layout="circular", branch.length='none', root.position = 10)
  # Plotting with adjusted offset and width for larger gaps
  p2 <- gheatmap(circ, dat_mlb, offset=0, width=0.1, 
                 colnames_angle=90, colnames_offset_y=0.25, colnames=FALSE) +
    scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white", name=traits[counter]) +  # Use the trait name directly
    theme(legend.position = "top", 
          legend.justification = "right",
          legend.title.align = 0.5,  # Center the legend title
          legend.spacing.x = unit(0.5, "cm"),  # Add space between legend and plot
          legend.margin = margin(10, 10, 10, 10))  # Add space around the legend
  
  #print(p2)
  
  
  print(p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 4))
  # Read in alleles for AGRP
  af <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ancestor_matching/top_rwas_all_99.9_and_pip01.allthegoodstuff.txt", header=T, sep="\t", check.names = FALSE)
  gt <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/ancestor_matching/top_rwas_all_99.9_and_pip01.GT.txt", header=T, sep="\t", check.names = FALSE)
  
  focus <- traits_short[counter]
  print(focus)
  
  topsnps <- head(af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ], n=5)
  #topsnps <- af[af$trait %in% focus, ][order(-af[af$trait %in% focus, ]$gamma), ]
  topsnps1 <- topsnps
  dim(topsnps1)
  
}
head(topsnps1)  
agrp <- gt %>%
  semi_join(topsnps1, by = c("CHROM" = "chr", "POS" = "pos")) %>%
  arrange(match(paste(CHROM, POS), paste(topsnps1$chr, topsnps1$pos)))
agrp1 <- data.frame(t(agrp), stringsAsFactors = FALSE)
colnames(agrp1) <- paste(agrp1[1, ], agrp1[2, ], sep = "_")
olnames(agrp1[1]) <- "chr11_632363" # because this column header is weeeird
agrp1 <- agrp1[-c(1,2), , drop = FALSE]  # Drop first two rows but keep structure
agrp3 <- merge(agrp1, species_df, by="row.names", all.x=TRUE)
agrp4 <- agrp3[!is.na(agrp3$short), ]
rownames(agrp4) <- agrp4$short
agrp5 <- agrp4[,c(2:7)]
agrp5[] <- lapply(agrp5, function(x) as.numeric(as.character(x)))
agrp5[agrp5 == 0] <- "A"
agrp5[agrp5 == 1] <- "B"
agrp5[agrp5 == 2] <- "C"
head(agrp5)
# Create plot p3 for each SNP

agrp5_matrix <- as.matrix(agrp5)
circ_p1 <- gheatmap(circ, dat_mlb, offset = 0.5, width = .06,
                    colnames_angle = 95, colnames_offset_y = 0.1, colnames = F) +
  scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p1 
circ_p1 <- circ_p1 + new_scale_fill()
circ_p2 <- gheatmap(circ_p1, agrp5_matrix, offset=1.5,  width=.24,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + 
  scale_fill_manual(values = c("grey90", "lightblue", "blue"), na.value = "white", labels = c("00", "01", "11"))
circ_p2


p3 <- circ_p2 + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 7.1)
print(p3)


# save figs with tips and legend
svg("figure3_kd_redundancy_tips_legend_short.svg",h=8,w=8)
p3
dev.off()






                    