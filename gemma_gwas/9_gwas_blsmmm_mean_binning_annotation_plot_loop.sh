for i in `cat trait`;

do

Rscript --vanilla 8_gwas_blsmmm_mean_binning_annotation_plot.R $i

done
