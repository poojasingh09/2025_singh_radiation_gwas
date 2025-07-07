while read line
do
   var1=`eval echo $line | awk '{print $1}'`
   var2=`eval echo $line | awk '{print $2}'`
   echo $var1
   echo $var2

   sed "s/lilo/$var1/g" 2025_ld_r2_generic.R | sed "s/stitch/$var2/g" > ld_r2_generic_$var1.R
   Rscript --vanilla ld_r2_generic_$var1.R &

done < ld_pairs




