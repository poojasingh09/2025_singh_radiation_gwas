for i in `cat trait`;

do

Rscript --vanilla 7_gwas_blsmmm_mean_binning_annotation.R $i

done
