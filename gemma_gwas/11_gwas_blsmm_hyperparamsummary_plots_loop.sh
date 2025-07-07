for j in `cat trait`;

do

Rscript --vanilla 10_gwas_blsmm_hyperparamsummary_plots.R $j ""

done
