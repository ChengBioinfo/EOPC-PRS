# run for General-, EO- and LO-population
for i in $(seq 1 22);do
plink --bfile ukb_imp_chr$i --clump cox_result.txt --clump-p1 1e-5 --clump-p2 1e-5 --clump-r2 0.8 --out clump_chr$i;
done