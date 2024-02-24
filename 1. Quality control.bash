# run for General-, EO- and LO-population
for i in $(seq 1 22);do
plink --allow-no-sex --bfile ukb_imp_chr$i --geno 0.05 --maf 0.01 --hwe 1e-6 --out chr$i --make-bed --keep 0.General_IDs.txt;
done