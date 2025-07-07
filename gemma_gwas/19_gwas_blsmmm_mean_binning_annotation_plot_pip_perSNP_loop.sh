for i in `cat trait`;

do

Rscript --vanilla 18_gwas_blsmmm_mean_binning_annotation_plot_pip_perSNP.R $i

done
