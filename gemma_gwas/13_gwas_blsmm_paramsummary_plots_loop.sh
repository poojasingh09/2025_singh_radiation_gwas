for j in `cat trait`;

do

	Rscript --vanilla 12_gwas_blsmm_paramsummary_plots.R $j

done
