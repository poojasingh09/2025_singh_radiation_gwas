cat *fst_LD_r2_results.txt | grep -v Rsq > LD_r2_results_allpairs_combined.txt 
less LD_r2_results_allpairs_combined.txt
head -n 1 ecinctus_MG_vs_ecoprologous_blue_MG_20kb.windowed.weir.fst_LD_r2_results.txt > h
cat h LD_r2_results_allpairs_combined.txt  > a
mv a LD_r2_results_allpairs_combined.txt
