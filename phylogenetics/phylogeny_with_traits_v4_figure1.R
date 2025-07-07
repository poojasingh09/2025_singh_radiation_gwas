

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
dat <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/gwas_colour_phenotypes_v2.lvrs.NA_Tp.txt", sep="\t", row.names=1, header=T)
head(dat)
dat <- dat[,1:8]


#species <- c(rownames(dat), "71001") ## "71001" is T pharyngalis

            
dat[dat == 0] <- "A"
dat[dat == 1] <- "B"

dat1 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_jawshape/gwas_headshape_PC12_edited.lvrs.txt", row.names=1, header=T)
head(dat1)


dat2 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_bcftools_pnye3_dentition_repeat/2025_LVcichlid_traits_GWAS.NA.txt", row.names=1, header=T)
head(dat2)
#dat2["71001", ] <- c("NA", "NA", "NA", "NA")


species <- unique(c(rownames(dat), rownames(dat1), rownames(dat2)))
             

missing_species <- species[!species %in% tree$tip.label]
if (length(missing_species) > 0) {
  warning("The following species are not found in the tree: ", paste(missing_species, collapse = ", "))
}

pruned.tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% species)])

#pruned.tree <-drop.tip(tree,tree$tip.label[-match(species, tree$tip.label)])
plot(pruned.tree)

pruned.tree <- root(pruned.tree, "106985", resolve.root = TRUE)
is.rooted(pruned.tree) 


normal <- ggtree(pruned.tree)
normal_p1 <- gheatmap(normal, dat, offset=.8, width=.2,colnames_angle=95, colnames_offset_y = .25)  + scale_fill_viridis_d(option="D", name="discrete\nvalue")
normal + geom_tiplab(size=3)

pdf("normal_newphylogeny_withtraits_v1.pdf", h=50, w=10)
#normal_p1 + geom_tiplab(size=0) +geom_nodelab()
normal_p1 + geom_tiplab(size=3)
dev.off()

circ <- ggtree(pruned.tree, layout="circular", branch.length='none')
circ_p1 <- gheatmap(circ, dat1, offset=.2, width=.3,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_viridis_c(option="A", name="continuous\nvalue",na.value="white") 

circ_p2 <- circ_p1 + new_scale_fill()
p2 <- gheatmap(circ_p2, dat, offset=3, width=.2,
         colnames_angle=90, colnames_offset_y = 0.25, colnames = F) + scale_fill_manual(values = c("skyblue","orange"),  na.value = "white")
circ_p1 + geom_tiplab(size=2)

p2 

circ_p3 <- p2 + new_scale_fill()
p3 <- gheatmap(circ_p3, dat2, offset=7, width=.1,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")

p3 + geom_tiplab(size=2, offset=12)



### this is the final one for the paper!!!!

circ <- ggtree(pruned.tree, layout="circular", branch.length='none')
circ_p1 <- gheatmap(circ, dat2[1], offset=.1,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p1 <- circ_p1 + new_scale_fill()
circ_p1

circ_p2 <- gheatmap(circ_p1, dat2[2], offset=.6,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p2 <- circ_p2 + new_scale_fill()
circ_p2

circ_p3 <- gheatmap(circ_p2, dat2[3], offset=1.1,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p3 <- circ_p3 + new_scale_fill()
circ_p3

circ_p4 <- gheatmap(circ_p3, dat2[4], offset=1.6,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p4 <- circ_p4 + new_scale_fill()
circ_p4

circ_p5 <- gheatmap(circ_p4, dat1, offset=2.1, width=.04,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_viridis_c(option="viridis", name="continuous\nvalue",na.value="white") 
circ_p5 <- circ_p5 + new_scale_fill()
circ_p5

circ_p6 <- gheatmap(circ_p5, dat[c(6,7,8)], offset=3.1, width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("black","grey90"),  na.value = "white")
circ_p6 <- circ_p6 + new_scale_fill()
circ_p6


circ_p7 <- gheatmap(circ_p6, dat[2], offset=4.5,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("red","grey90"),  na.value = "white")
circ_p7 <- circ_p7 + new_scale_fill()
circ_p7 + theme(legend.position = "none")

circ_p8 <- gheatmap(circ_p7, dat[1], offset=5,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("yellow","grey90"),  na.value = "white")
circ_p8 <- circ_p8 + new_scale_fill()
circ_p8 + theme(legend.position = "none")

circ_p9 <- gheatmap(circ_p8, dat[3], offset=5.5,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("red","grey90"),  na.value = "white")
circ_p9 <- circ_p9 + new_scale_fill()
circ_p9 + theme(legend.position = "none")

circ_p10 <- gheatmap(circ_p9, dat[5], offset=6,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("black","grey90"),  na.value = "white")
circ_p10 <- circ_p10 + new_scale_fill()
circ_p10 + theme(legend.position = "none")

circ_p11 <- gheatmap(circ_p10, dat[4], offset=6.5,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("cadetblue1","grey90"),  na.value = "white")
circ_p11 <- circ_p11 + new_scale_fill()
circ_p11 + theme(legend.position = "none")

#ggsave(file="circular_phylogeny_withtraits_v2.png", plot=p2, width=10, height=8)
ggsave("circular_newphylogeny_withtraits_nolegend.svg", device = svglite::svglite,  width=10, height=8)

circ_p11 
ggsave("circular_newphylogeny_withtraits_legend.svg", device = svglite::svglite,  width=10, height=8)

############### lets add species names to the tips first

# Create a new data frame with species names
#tree <- ape::read.tree("all.chr.filt.2kbThinned.raxmlinput.min4.phy.fa.treefile.tre")
plot(tree,show.node=TRUE)
head(listTips(tree))
length(tree$tip.label)
#write.table(tree$tip.label, "all.chr.filt.2kbThinned.raxmlinput.min4.phy.fa.treefile.tre.tips", quote=F, row.names=F)


## read in species names, this file contains fishec id to species name mappingf for samples whose have WGS and/or microCT data, ok?
species_df0 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/phylogeny/species_name.txt", sep="\t", header=T, row.names=1)
colnames(species_df0) <- c("species_name") # make sure "71001" T pharyngalis is in this file as it is root
head(species_df0)
dim(species_df0)
rownames(species_df0) == "71001"
rownames(species_df0) == "106985"
species_df0[rownames(species_df0) == "106641", "species_name"] <- "Labrochromis sp. demersal" #otherwise it does not print!

## 1 ind per species (i only have 124, ole expects 136?!) 

unique <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/phylogeny/1_ind_per_species.txt", sep="\t", header=T, check.names = F)
head(unique)
dim(unique)
unique$sample == "71001"
unique$sample == "106985"

species_df_1 <- species_df0 %>% filter(rownames(species_df0) %in% unique$sample)
species_df_11 <- species_df0 %>% filter(species_df0$species_name %in% species_df_1$species_name)
dim(unique)
dim(species_df_1)
dim(species_df_11)
#dim(species_df)
rownames(species_df_11) == "71001"
rownames(species_df_11) == "106985"


# Function to add suffixes to duplicates so they are unique, only second onwards are labelled
add_suffix_to_duplicates <- function(names) {
  seen <- list()
  
  for (i in seq_along(names)) {
    name <- names[i]
    
    if (name %in% names[1:(i-1)]) {
      seen[[name]] <- seen[[name]] + 1
      names[i] <- paste0(name, "_", seen[[name]])
    } else {
      seen[[name]] <- 1  # Initialize counter for new names
    }
  }
  
  return(names)
}

# Apply the function to the data frame
species_df_2 <- species_df_1 %>%
  mutate(species_name = add_suffix_to_duplicates(species_name))


species_df_2 <- species_df0 %>%
  mutate(species_name = add_suffix_to_duplicates(species_name))


head(species_df_2)
rownames(species_df_2) == "71001"
rownames(species_df_2) == "106985"

species_df <- species_df_2
species_df[1,1] <- "Double-stripe group tanaos_1" #for some reason this has an extra underscore so i renamed it

write.table(species_df, "species_names_newphylogeny_numbered.txt", quote=F)

## subset phylogeny for only tips in my data
# Get the tips to keep

pruned.tree <- drop.tip(tree,tree$tip.label[-match(rownames(species_df), tree$tip.label)])
plot(pruned.tree)

############# new additions may 2025
# Find matching indices (can include NAs if not found)
matched_indices <- match(rownames(species_df), tree$tip.label)

# Remove NAs (i.e., rownames that didn't match any tip label)
matched_indices <- matched_indices[!is.na(matched_indices)]

# Drop tips that are NOT in the matched indices
pruned.tree <- drop.tip(tree, tree$tip.label[-matched_indices])
##############


pruned.tree <- root(pruned.tree, "84961", resolve.root = TRUE) #84961 is A. bloyeti and 106985 is T pharyngalis
is.rooted(pruned.tree) 
plot(pruned.tree)

# Match species names to pruned tree tips and reset tiplabels
species_names <- species_df$species_name[match(pruned.tree$tip.label, rownames(species_df))]
pruned.tree$tip.label <- species_names

######## adapt traits

#dat <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_colour/gwas_colour_v1.NA_binary.lvrs.matchedphylogeny.txt", row.names=1, header=T)
dat <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/gwas_colour_phenotypes_v2.lvrs.NA_Tp.txt", sep="\t", row.names=1, header=T)
head(dat)
dat <- dat[,1:8]

dat[dat == 0] <- "A"
dat[dat == 1] <- "B"


dat0 <- merge(dat, species_df, by="row.names", all.x=FALSE)
dat <- dat0[,c(2:9)]
rownames(dat) <- dat0$species_name

dat1 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_jawshape/gwas_headshape_PC12_plottingedited_dupsremoved.txt", row.names=1, header=T, sep="\t")
#dat1 <- rbind(dat1, data.frame(PC1 = -0.040229339, PC2 = 0.054773309, row.names = "71001")) ## data for outgroup pharyngalis
head(dat1)
dat1 <- dat1[,c(7,8)]
a <- merge(dat1, species_df0, by="row.names", all.x=FALSE)
a_unique <- a[!duplicated(a$species_name), ]
dat0 <- merge(a_unique, species_df, by = "species_name", all.x = FALSE)
dat1 <- dat0[,c(3:4)]
rownames(dat1) <- dat0$species_name

dat2 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_bcftools_pnye3_dentition_repeat/2025_LVcichlid_traits_GWAS.NA.txt", row.names=1, header=T)
dat0 <- merge(dat2, species_df, by="row.names", all.x=FALSE) ## there are some species missing here, come back to this ok?
dat2 <- dat0[,c(2:5)]
rownames(dat2) <- dat0$species_name


head(dat2) # this has dentision data
head(dat1) # this has 2 columns of continuous data PC1 and PC2
head(dat) # this is binary data for 8 traits

dim(dat2)
dim(dat1)
dim(dat)



# circular dendo

circ <- ggtree(pruned.tree, layout="circular", branch.length='none')
circ
circ_p1 <- gheatmap(circ, dat2[1], offset=.1,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p1 <- circ_p1 + new_scale_fill()
circ_p1

# Get coordinates for tips (which we will use for drawing lines)
tip_data <- as_tibble(circ_p1$data) %>% filter(isTip == TRUE)

# Add dotted lines connecting tips to annotation rows
circ_p1 + 
  geom_segment(data = tip_data, 
               aes(x = x, y = y, xend = x + 0.1, yend = y),  # Adjust xend and yend for annotation row positions
               linetype = "dotted", color = "black", size = 0.5) + 
  theme(legend.position = "none")

circ_p2 <- gheatmap(circ_p1, dat2[2], offset=.4,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p2 <- circ_p2 + new_scale_fill()
circ_p2

circ_p3 <- gheatmap(circ_p2, dat2[3], offset=.7,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p3 <- circ_p3 + new_scale_fill()
circ_p3

circ_p4 <- gheatmap(circ_p3, dat2[4], offset=1,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p4 <- circ_p4 + new_scale_fill()
circ_p4

circ_p5 <- gheatmap(circ_p4, dat1, offset=1.3, width=.04,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_viridis_c(option="viridis", name="continuous\nvalue",na.value="white", direction = -1) 
circ_p5 <- circ_p5 + new_scale_fill()
circ_p5


circ_p6 <- gheatmap(circ_p5, dat[c(6,7,8)], offset=1.9, width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90","black"),  na.value = "white")
circ_p6 <- circ_p6 + new_scale_fill()
circ_p6


circ_p7 <- gheatmap(circ_p6, dat[2], offset=2.8,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "red"),  na.value = "white")
circ_p7 <- circ_p7 + new_scale_fill()
circ_p7 + theme(legend.position = "none")

circ_p8 <- gheatmap(circ_p7, dat[1], offset=3.1,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "yellow"),  na.value = "white")
circ_p8 <- circ_p8 + new_scale_fill()
circ_p8 + theme(legend.position = "none")

circ_p9 <- gheatmap(circ_p8, dat[3], offset=3.4,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "red"),  na.value = "white")
circ_p9 <- circ_p9 + new_scale_fill()
circ_p9 + theme(legend.position = "none")

circ_p10 <- gheatmap(circ_p9, dat[5], offset=3.7,  width=.02,
                     colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "black"),  na.value = "white")
circ_p10 <- circ_p10 + new_scale_fill()
circ_p10 + theme(legend.position = "none")

circ_p11 <- gheatmap(circ_p10, dat[4], offset=4,  width=.02,
                     colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90","cadetblue1"),  na.value = "white")
circ_p11 <- circ_p11 + new_scale_fill()
circ_p11 + theme(legend.position = "none") + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE, offset = 5)



#ggsave(file="circular_phylogeny_withtraits_v2.png", plot=p2, width=10, height=8)
ggsave("circular_newphylogeny_withtraits_nolegend_tips_space.svg", device = svglite::svglite,  width=14, height=14)


circ_p11 + 
  theme(legend.position = "none") + 
  geom_tiplab(fontface = "italic", size = 3, inherit.aes = FALSE, offset = 4.5)

#ggsave(file="circular_phylogeny_withtraits_v2.png", plot=p2, width=10, height=8)
ggsave("circular_newphylogeny_withtraits_nolegend_tips.svg", device = svglite::svglite,  width=20, height=20)


circ_p11 + geom_tiplab(fontface="italic", size=1.5, inherit.aes=FALSE, offset = 7)
ggsave("circular_newphylogeny_withtraits_legend_tips.svg", device = svglite::svglite,  width=10, height=8)

## circular with branch length

circ <- ggtree(pruned.tree, layout="circular")
circ_p1 <- gheatmap(circ, dat2[1], offset=.01,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p1 <- circ_p1 + new_scale_fill()
circ_p1


# Create the circular tree plot
circ <- ggtree(pruned.tree, layout = "circular")

# Add the heatmap annotation
circ_p1 <- gheatmap(circ, dat2[1], offset = 0.01, width = .06,
                    colnames_angle = 95, colnames_offset_y = 0.1, colnames = FALSE) +
  scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p1 

circ_p2 <- gheatmap(circ_p1, dat2[2], offset=.04,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p2 <- circ_p2 + new_scale_fill()
circ_p2

circ_p3 <- gheatmap(circ_p2, dat2[3], offset=.07,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p3 <- circ_p3 + new_scale_fill()
circ_p3

circ_p4 <- gheatmap(circ_p3, dat2[4], offset=0.1,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p4 <- circ_p4 + new_scale_fill()
circ_p4

circ_p5 <- gheatmap(circ_p4, dat1, offset=0.13, width=.12,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_viridis_c(option="A", name="continuous\nvalue",na.value="white", direction = -1) 
circ_p5 <- circ_p5 + new_scale_fill()
circ_p5


circ_p6 <- gheatmap(circ_p5, dat[c(6,7,8)], offset=0.19, width=.18,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90","#009E73"),  na.value = "white")
circ_p6 <- circ_p6 + new_scale_fill()
circ_p6

circ_p7 <- gheatmap(circ_p6, dat[2], offset=0.28,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "red"),  na.value = "white")
circ_p7 <- circ_p7 + new_scale_fill()
circ_p7 + theme(legend.position = "none")

circ_p8 <- gheatmap(circ_p7, dat[1], offset=0.31,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "yellow"),  na.value = "white")
circ_p8 <- circ_p8 + new_scale_fill()
circ_p8 + theme(legend.position = "none")

circ_p9 <- gheatmap(circ_p8, dat[3], offset=0.34,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "red"),  na.value = "white")
circ_p9 <- circ_p9 + new_scale_fill()
circ_p9 + theme(legend.position = "none")

circ_p10 <- gheatmap(circ_p9, dat[5], offset=0.37,  width=.06,
                     colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "black"),  na.value = "white")
circ_p10 <- circ_p10 + new_scale_fill()
circ_p10 + theme(legend.position = "none")

circ_p11 <- gheatmap(circ_p10, dat[4], offset=0.40,  width=.06,
                     colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90","cadetblue1"),  na.value = "white")
circ_p11 <- circ_p11 + new_scale_fill()

circ_p11 + 
  theme(legend.position = "none") + 
  geom_tiplab(aes(x = 1, y = y), fontface = "italic", size = 5, inherit.aes = FALSE)

#ggsave(file="circular_phylogeny_withtraits_v2.png", plot=p2, width=10, height=8)
ggsave("circular_newphylogeny_withtraits_nolegend_tips.svg", device = svglite::svglite, width=24, height=23)

circ_p11 + 
  theme(legend.position = "none")
ggsave("circular_newphylogeny_withtraits_nolegend_notips.svg", device = svglite::svglite,  width=25, height=25)

circ_p11
ggsave("circular_phylogeny_withtraits_legend_notips.svg", device = svglite::svglite,  width=25, height=25)




## rectangular with branch length is dendogram
circ <- ggtree(pruned.tree, layout="rectangular")
circ_p1 <- gheatmap(circ, dat2[1], offset = 0.01, width = .02,
                    colnames_angle = 95, colnames_offset_y = 0.1, colnames = FALSE) +
  scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p1 

circ_p2 <- gheatmap(circ_p1, dat2[2], offset=.02,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p2 <- circ_p2 + new_scale_fill()
circ_p2

circ_p3 <- gheatmap(circ_p2, dat2[3], offset=.03,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p3 <- circ_p3 + new_scale_fill()
circ_p3

circ_p4 <- gheatmap(circ_p3, dat2[4], offset=0.04,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p4 <- circ_p4 + new_scale_fill()
circ_p4

circ_p5 <- gheatmap(circ_p4, dat1, offset=0.05, width=.04,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_viridis_c(option="A", name="continuous\nvalue",na.value="white", direction = -1) 
circ_p5 <- circ_p5 + new_scale_fill()
circ_p5


circ_p6 <- gheatmap(circ_p5, dat[c(6,7,8)], offset=0.07, width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90","#009E73"),  na.value = "white")
circ_p6 <- circ_p6 + new_scale_fill()
circ_p6

circ_p7 <- gheatmap(circ_p6, dat[2], offset=0.102,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "red"),  na.value = "white")
circ_p7 <- circ_p7 + new_scale_fill()
circ_p7 + theme(legend.position = "none")

circ_p8 <- gheatmap(circ_p7, dat[1], offset=0.112,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "yellow"),  na.value = "white")
circ_p8 <- circ_p8 + new_scale_fill()
circ_p8 + theme(legend.position = "none")

circ_p9 <- gheatmap(circ_p8, dat[3], offset=0.122,  width=.02,
                    colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "red"),  na.value = "white")
circ_p9 <- circ_p9 + new_scale_fill()
circ_p9 + theme(legend.position = "none")

circ_p10 <- gheatmap(circ_p9, dat[5], offset=0.132,  width=.02,
                     colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90", "black"),  na.value = "white")
circ_p10 <- circ_p10 + new_scale_fill()
circ_p10 + theme(legend.position = "none")

circ_p11 <- gheatmap(circ_p10, dat[4], offset=0.142,  width=.02,
                     colnames_angle=95, colnames_offset_y = .1, colnames = F) + scale_fill_manual(values = c("grey90","cadetblue1"),  na.value = "white")
circ_p11 <- circ_p11 + new_scale_fill()

circ_p11 + 
  theme(legend.position = "none") + 
  geom_tiplab(aes(x = 0.705, y = y), fontface = "italic", size = 2, inherit.aes = FALSE) + coord_cartesian(clip="off")

#ggsave(file="circular_phylogeny_withtraits_v2.png", plot=p2, width=14, height=8)
ggsave("circular_newphylogeny_withtraits_nolegend_tips_rect.svg", device = svglite::svglite,  width=36, height=14)

circ_p11 + geom_tiplab(aes(x = 0.7, y = y), fontface = "italic", size = 2, inherit.aes = FALSE) + coord_cartesian(clip="off")
ggsave("circular_newphylogeny_withtraits_legend_tips_rect.svg", device = svglite::svglite,  width=36, height=14)


## circular with branch length and bootstrap support

# Create the circular tree plot
#circ <- ggtree(pruned.tree, layout = "circular")

circ <- ggtree(pruned.tree, layout = "circular") +
  geom_point(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 50, 
                 x = x, y = y, size = as.numeric(label)), 
             shape = 21, fill = "black", color = "black", alpha = 0.4) +
  scale_size_continuous(range = c(1, 3)) +  # Reduce circle sizes
  theme(legend.position = "right")

print(circ)


#reset colnames for appropriate plotting
# Add the heatmap annotation
# First heatmap
circ_p1 <- gheatmap(circ, dat2[1], offset = 0.01, width = .06, 
                    colnames_angle = 95, colnames_offset_y = 0.2, colnames = T) + 
  scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")


circ_p1
# Reset scale before adding second heatmap
circ_p2 <- circ_p1 + new_scale_fill()  

# Second heatmap with independent scale
circ_p2 <- gheatmap(circ_p2, dat2[2], offset = .04, width = .06, 
                    colnames_angle = 95, colnames_offset_y = .1, colnames = TRUE) + 
  scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")

circ_p2

circ_p3 <- circ_p2 + new_scale_fill() 
circ_p3 <- gheatmap(circ_p3, dat2[3], offset=.07,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p3

circ_p4 <- circ_p3 + new_scale_fill()
circ_p4 <- gheatmap(circ_p4, dat2[4], offset=0.1,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p4

circ_p5 <- circ_p4 + new_scale_fill()
circ_p5 <- gheatmap(circ_p5, dat1, offset=0.13, width=.12,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_viridis_c(option="A", name="continuous\nvalue",na.value="white", direction = -1) 
circ_p5

circ_p6 <- circ_p5 + new_scale_fill()
circ_p6 <- gheatmap(circ_p6, dat[c(6,7,8)], offset=0.19, width=.18,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("grey90","#009E73"),  na.value = "white")
circ_p6

circ_p7 <- circ_p6 + new_scale_fill()
circ_p7 <- gheatmap(circ_p7, dat[2], offset=0.28,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("grey90", "red"),  na.value = "white")
circ_p7 + theme(legend.position = "none")

circ_p8 <- circ_p7 + new_scale_fill()
circ_p8 <- gheatmap(circ_p8, dat[1], offset=0.31,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("grey90", "yellow"),  na.value = "white")
circ_p8 + theme(legend.position = "none")

circ_p9 <- circ_p8 + new_scale_fill()
circ_p9 <- gheatmap(circ_p9, dat[3], offset=0.34,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("grey90", "red"),  na.value = "white")
circ_p9 + theme(legend.position = "none")

circ_p10 <- circ_p9 + new_scale_fill()
circ_p10 <- gheatmap(circ_p10, dat[5], offset=0.37,  width=.06,
                     colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("grey90", "black"),  na.value = "white")
circ_p10 + theme(legend.position = "none")

circ_p11 <- circ_p10 + new_scale_fill()
circ_p11 <- gheatmap(circ_p11, dat[4], offset=0.40,  width=.06,
                     colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("grey90","cadetblue1"),  na.value = "white")
circ_p11 + 
  theme(legend.position = "none") + 
  geom_tiplab(aes(x = 1, y = y), fontface = "italic", size = 5, inherit.aes = FALSE)

#ggsave(file="circular_phylogeny_withtraits_v2.png", plot=p2, width=10, height=8)
ggsave("circular_newphylogeny_withtraits_nolegend_tips.svg", device = svglite::svglite, width=20, height=20)

circ_p11 + 
  theme(legend.position = "none")
ggsave("circular_newphylogeny_withtraits_nolegend_notips.svg", device = svglite::svglite,  width=20, height=20)

circ_p11
ggsave("circular_newphylogeny_withtraits_legend_notips.svg", device = svglite::svglite,  width=20, height=20)





####
####
####
###
###
##
##
#
#




## circular with branch length and bootstrap support########################## now with SHORT names for your convenience
############### lets add species names to the tips first

# Create a new data frame with species names
tree <- ape::read.tree("all.filt.vcf.var.1kbthinned.gtr.dna.treefile")
plot(tree,show.node=TRUE)
head(listTips(tree))
length(tree$tip.label)
#write.table(tree$tip.label, "all.chr.filt.2kbThinned.raxmlinput.min4.phy.fa.treefile.tre.tips", quote=F, row.names=F)


## read in species names
species_df_o <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/2023_sampleselection_gwas_singh_manuscriptsubset_OS_PS_removedsamples.txt", sep="\t", header=T, row.names=1)
colnames(species_df_o) <- c("species_name", "short") 
head(species_df_o)
dim(species_df_o)

## read in species names, this file contains fishec id to species name mappingf for samples whose have WGS and/or microCT data, ok?
species_df0 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/phylogeny/species_name.txt", sep="\t", header=T, row.names=1)
colnames(species_df0) <- c("species_name") # make sure "71001" T pharyngalis is in this file as it is root
head(species_df0)
dim(species_df0)
rownames(species_df0) == "71001"
rownames(species_df0) == "131282"
species_df0[rownames(species_df0) == "106641", "species_name"] <- "Labrochromis sp. demersal" #otherwise it does not print!

## 1 ind per species (i only have 124, ole expects 136?!) 

unique <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_overlap_pip01_q99.9/phylogeny/1_ind_per_species.txt", sep="\t", header=T, check.names = F)
head(unique)
dim(unique)
unique$sample == "71001"
unique$sample == "131282"

unique <- c(unique$sample, "71001", "131282") ## "71001" is T pharyngalis, is Aburtoni 131282
unique1 <- unique[unique != "106985"] # remove pseudocrenilabrus multicolor victoriae 106985
unique <- unique1
unique

species_df_1 <- species_df0 %>% filter(rownames(species_df0) %in% unique)
species_df_11 <- merge(species_df_1, species_df_o, by="row.names", all.x=T, all.y=F)
species_df_11[34,] <- c("11965", "Hoplotilapia retrodens","Hoplotilapia retrodens","HopRet") # these are NAs
species_df_11[87,] <- c("79628","Macropleurodus bicolor","Macropleurodus bicolor","MacBic") # be very careful
species_df_11[60,] <- c("131282","Astatotilapia burtoni","Astatotilapia burtoni","AstBur") # renaming things double check like crazy
species_df_11[90,] <- c("81343","Astatotilapia stappersi","Astatotilapia stappersi","AstSta")
dim(unique)
dim(species_df_1)
dim(species_df_11)
dim(species_df_o)
species_df_11$Row.names == "71001"
species_df_11$Row.names == "131282"

## check species_df11 here 

# Function to add suffixes to duplicates so they are unique, only second onwards are labelled
add_suffix_to_duplicates <- function(names) {
  seen <- list()
  
  for (i in seq_along(names)) {
    name <- names[i]
    
    if (name %in% names[1:(i-1)]) {
      seen[[name]] <- seen[[name]] + 1
      names[i] <- paste0(name, "_", seen[[name]])
    } else {
      seen[[name]] <- 1  # Initialize counter for new names
    }
  }
  
  return(names)
}

# Apply the function to the data frame
species_df_2 <- species_df_11 %>%
  mutate(short = add_suffix_to_duplicates(short))

head(species_df_2)
rownames(species_df_2) == "71001"
rownames(species_df_2) == "131282"

species_df <- species_df_2
head(species_df)
species_df[1,4] <- "PunNye"
rownames(species_df) <- species_df$Row.names

write.table(species_df[,c(1,2,4)], "species_names_newphylogeny.txt", quote=F, sep="\t")

## subset phylogeny for only tips in my data
# Get the tips to keep

#pruned.tree <- drop.tip(tree,tree$tip.label[-match(rownames(species_df), tree$tip.label)]) # this code gives issues
#plot(pruned.tree)

# Get tips to keep: those that exist in both tree and species_df
tips_to_keep <- intersect(tree$tip.label, rownames(species_df))
# Drop tips not in the intersection
pruned.tree <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))
pruned.tree <- root(pruned.tree, "131282", resolve.root = TRUE) #131282 burtoni
is.rooted(pruned.tree) 



# Match species names to pruned tree tips and reset tiplabels
short<- species_df$short[match(pruned.tree$tip.label, rownames(species_df))]
pruned.tree$tip.label <- short

######## adapt traits

#dat <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_colour/gwas_colour_v1.NA_binary.lvrs.matchedphylogeny.txt", row.names=1, header=T)
dat <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2024_bcftools_pnye3_repeatcolour/gwas_colour_phenotypes_v2.lvrs.NA_Tp.txt", sep="\t", row.names=1, header=T)
head(dat)
dat <- dat[,1:8]

dat[dat == 0] <- "A"
dat[dat == 1] <- "B"


dat0 <- merge(dat, species_df, by="row.names",all.x=FALSE)
dat <- dat0[,c(2:9)]
rownames(dat) <- dat0$short


### PCs are special because we don't have have the microct data for those that are sequenced, but we want the sepcies value for plotting so lets do this
dat1 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2023_bcftools_pnye3/gwas_jawshape/gwas_headshape_PC12_plottingedited_dupsremoved.txt", row.names=1, header=T, sep="\t")
#dat1 <- rbind(dat1, data.frame(PC1 = -0.040229339, PC2 = 0.054773309, row.names = "71001")) ## data for outgroup pharyngalis
head(dat1)
dat1 <- dat1[,c(7,8)]
a <- merge(dat1, species_df0, by="row.names", all.x=FALSE)
a_unique <- a[!duplicated(a$species_name), ]
a_unique$species_name[15] <- "Labrochromis sp. demersal MG"
dat0 <- merge(a_unique, species_df, by.x = "species_name",by.y = "species_name.x",  all.x = FALSE)
dat1 <- dat0[,c(3:4)]
rownames(dat1) <- dat0$short

dat2 <- read.table("/Users/singhpoo/Desktop/postdoc_bern_2021/postdoc_researchprojects/genomics/2022_gwas_project/2025_bcftools_pnye3_dentition_repeat/2025_LVcichlid_traits_GWAS.NA.txt", row.names=1, header=T)
dat0 <- merge(dat2, species_df, by="row.names", all.x=FALSE) ## there are some species missing here, come back to this ok?
dat2 <- dat0[,c(2:5)]
rownames(dat2) <- dat0$short


head(dat2) # this has dentision data
head(dat1) # this has 2 columns of continuous data PC1 and PC2
head(dat) # this is binary data for 8 traits

dim(dat2)
dim(dat1)
dim(dat)

write.tree(pruned.tree, file = "prunedtree.tree")

# Create the circular tree plot
#circ <- ggtree(pruned.tree, layout = "circular")
library(stringr)

#circ <- ggtree(pruned.tree, layout = "circular") +
  #geom_point(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 50, 
   #              x = x, y = y, size = as.numeric(label)), 
    #         shape = 21, fill = "black", color = "black", alpha = 0.4) +
  #scale_size_continuous(range = c(1, 3)) +  # Reduce circle sizes
  #theme(legend.position = "right")

#print(circ)
# Step 1: Plot the circular tree
circ <- ggtree(pruned.tree, layout = "circular")

# Step 2: Extract the ggtree data
tree_data <- circ$data

# Step 3: Filter internal nodes and extract the first number from label
tree_data <- subset(tree_data, isTip == FALSE)
tree_data$support <- as.numeric(str_extract(tree_data$label, "^[0-9.]+"))

# Step 4: Filter support values >= 50
filtered_nodes <- subset(tree_data, !is.na(support) & support >= 50)

# Step 5: Add circles to the plot using filtered data
circ <- circ +
  geom_point(data = filtered_nodes,
             aes(x = x, y = y, size = support),
             shape = 21, fill = "black", color = "black", alpha = 0.4) +
  scale_size_continuous(range = c(1, 3)) +  # Reduce circle sizes
  theme(legend.position = "right")

# Step 6: Show the plot
print(circ)



#reset colnames for appropriate plotting
# Add the heatmap annotation
# First heatmap
circ_p1 <- gheatmap(circ, dat2[1], offset = 0.11, width = .06, 
                    colnames_angle = 95, colnames_offset_y = 0.2, colnames = T) + 
  scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")


circ_p1


# Reset scale before adding second heatmap
circ_p2 <- circ_p1 + new_scale_fill()  

# Second heatmap with independent scale
circ_p2 <- gheatmap(circ_p2, dat2[2], offset = .14, width = .06, 
                    colnames_angle = 95, colnames_offset_y = .1, colnames = TRUE) + 
  scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")

circ_p2

circ_p2 + 
  theme(legend.position = "none") + 
  geom_tiplab(aes(x = 0.55, y = y), fontface = "italic", size = 2, inherit.aes = FALSE)



circ_p3 <- circ_p2 + new_scale_fill() 
circ_p3 <- gheatmap(circ_p3, dat2[3], offset=.17,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p3

circ_p4 <- circ_p3 + new_scale_fill()
circ_p4 <- gheatmap(circ_p4, dat2[4], offset=0.20,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_gradient(low = "grey90", high = "grey30", na.value = "white")
circ_p4

circ_p5 <- circ_p4 + new_scale_fill()
circ_p5 <- gheatmap(circ_p5, dat1, offset=0.23, width=.12,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_viridis_c(option="A", name="continuous\nvalue",na.value="white", direction = -1) 
circ_p5

circ_p6 <- circ_p5 + new_scale_fill()


circ_p6 <- gheatmap(circ_p6, dat[c(6,7,8)], offset=0.30, width=.18,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("grey90","#009E73"),  na.value = "white")
circ_p6

circ_p7 <- circ_p6 + new_scale_fill()
circ_p7 <- gheatmap(circ_p7, dat[2], offset=0.405,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("grey90", "red"),  na.value = "white")
circ_p7 + theme(legend.position = "none")

circ_p8 <- circ_p7 + new_scale_fill()
circ_p8 <- gheatmap(circ_p8, dat[1], offset=0.441,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("cadetblue1", "yellow"),  na.value = "white")
circ_p8 + theme(legend.position = "none")

circ_p9 <- circ_p8 + new_scale_fill()
circ_p9 <- gheatmap(circ_p9, dat[3], offset=0.475,  width=.06,
                    colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("grey90", "red"),  na.value = "white")
circ_p9 + theme(legend.position = "none")

circ_p10 <- circ_p9 + new_scale_fill()
circ_p10 <- gheatmap(circ_p10, dat[5], offset=0.51,  width=.06,
                     colnames_angle=95, colnames_offset_y = .1, colnames = T) + scale_fill_manual(values = c("grey90", "black"),  na.value = "white")
circ_p10 + theme(legend.position = "none")

circ_p11 <- circ_p10 

circ_p11 + 
  theme(legend.position = "none") + 
  geom_tiplab(aes(x = 0.595, y = y), fontface = "italic", size = 3.8, inherit.aes = FALSE)


#ggsave(file="circular_phylogeny_withtraits_v2.png", plot=p2, width=10, height=8)
ggsave("circular_newphylogeny_withtraits_nolegend_tips_short.svg", device = svglite::svglite, width=20, height=20)

circ_p11 + 
  theme(legend.position = "none")
ggsave("circular_newphylogeny_withtraits_nolegend_notips_short.svg", device = svglite::svglite,  width=20, height=20)

circ_p11
ggsave("circular_newphylogeny_withtraits_legend_notips_short.svg", device = svglite::svglite,  width=20, height=20)


####### testing testing

### this is to root with multiple nodes #####
mrca_node <- getMRCA(tree, c("T_B","78584","78357", "106985")) # thorachochromis buysi, S. checkerboard, O. red-cheek, P. multicolor victoriae
tree1 <- root(tree, node = mrca_node, resolve.root = TRUE)

species_df0$species_name <- paste(rownames(species_df0), species_df0$species_name)
full<- species_df0$short[match(tree$tip.label, rownames(species_df0))]
tree$tip.label <- full


circ <- ggtree(tree, layout = "rectangular", branch.length = TRUE)
circ + theme(legend.position = "none") + geom_tiplab(fontface="italic", size=2, inherit.aes=FALSE)
circ + 
  theme(legend.position = "none") + 
  geom_tiplab(aes(label = label), fontface = "italic", size = 2, inherit.aes = FALSE) +
  geom_text2(aes(label = branch.length), size = 2, vjust = -0.5)
                    