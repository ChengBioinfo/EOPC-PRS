plink --bfile /database/database1/nas/public_data/5.PLCO/PCa/merge_CGEMS_PEGASUS/merge_7776 --exclude dupSNP.txt --make-bed --out ../PLCO/merge_7776_dedup
plink --bfile  /database/database1/nas/public_data/4.TCGA/TCGA_PRAD/SNP_impute/bmerge_atcg_info --exclude dupSNP.txt --make-bed --out ../TCGA/bmerge_atcg_info_dedup
plink --bfile /data2/cyf/EOPC_rep/TCGA/bmerge_atcg_info_dedup --keep ../TCGA/white.txt --make-bed --out TCGA_white

#merge
plink --bfile /data2/cyf/EOPC_rep/PLCO/merge_7776_dedup --bmerge TCGA_white.bed TCGA_white.bim TCGA_white.fam --make-bed --out merge
plink --bfile merge_7776_dedup --bmerge TCGA_white.bed TCGA_white.bim TCGA_white.fam --make-bed --out merge

#extract SNP
plink --bfile merge --write-snplist --out snp
plink --bfile /database/database1/nas/public_data/2.1000Genome_exp/1KG_Data_2012_March_Yu_Hongjie/1000genome --extract snp.snplist --out 1KG_pca --make-bed
plink --bfile 1KG_pca --write-snplist --out 1KG_pca
plink --bfile merge --extract 1KG_pca.snplist --out cal_pca --make-bed

#flip and merge
plink --bfile cal_pca --bmerge 1KG_pca.bed 1KG_pca.bim 1KG_pca.fam --out pca_merge --make-bed
plink --bfile 1KG_pca --exclude pca_merge-merge.missnp --out 1KG_pca_delt --make-bed
plink --bfile cal_pca --bmerge 1KG_pca_delt.bed 1KG_pca_delt.bim 1KG_pca_delt.fam --out pca_merge --make-bed

#thinld
plink --bfile pca_merge --maf 0.25 --out pca_merge_thin --make-bed
plink --bfile pca_merge_thin --indep-pairwise 1500 150 0.05 --out pca_merge_thinld
plink --bfile pca_merge_thin --extract pca_merge_thinld.prune.in --out pca_merge_thinld --make-bed

plink --bfile pca_merge_thinld --pca 10 --out PCA_results 

#PCA
plink --threads 20 --bfile merge --pca 10 --out pca --allow-no-sex --pheno T_pheno --pheno-name HDL_High_density_lipoprotein__mmol_l 



###EIGENSOFT extract outliers

#QC
plink --bfile ../merge --maf 0.25 --out merge_thin --make-bed
plink --bfile merge_thin --indep-pairwise 1500 150 0.05 --out merge_thinld
plink --bfile merge_thin --extract merge_thinld.prune.in --out merge_thinld --make-bed
plink --bfile merge_thinld --exclude del_SNP.txt --out merge_thinld_del --recode --make-bed

cp merge_thinld_del.bim merge_thinld_del.pedsnp
cp merge_thinld_del.fam merge_thinld_del.pedind

smartpca.perl -i merge_thinld_del.ped -a merge_thinld_del.pedsnp -b merge_thinld_del.pedind -o merge_thinld_del.pca -p merge_thinld_del.plot -e merge_thinld_del.eval -l merge_thinld_del > merge_thinld_del.out

perl ploteig -i merge_thinld_del_pop.pca.evec -c 1:2 -p PLCO_Control:PLCO_Case:TCGA:Case -o merge_thinld_del_pop.pca.plot

###extract SNPs
cd extract_SNPs
cp merge.fam merge.fam.old
------------------------------------

R
fam <- read.table("../merge.fam.old", sep=" ", header = F)
fam$V1 <- as.character(fam$V1)
fam$V2 <- as.character(fam$V2)
fam$V1[grep("TCGA", fam$V1)] <- paste0(substr(fam$V1[grep("TCGA", fam$V1)], 1, 12),
                                       substr(fam$V1[grep("TCGA", fam$V1)], 26, 28))
fam$V2[grep("TCGA", fam$V2)] <- paste0(substr(fam$V2[grep("TCGA", fam$V2)], 1, 12),
                                       substr(fam$V2[grep("TCGA", fam$V2)], 26, 28))
write.table(fam, "../merge.fam", sep = " ", quote = F, col.names = F, row.names = F)
---------------------------------------------------------------------------

plink --bfile ../merge --extract General_pos.txt --out General --make-bed --recode A --keep ../EIGENSOFT/retain_individual.txt
plink --bfile ../merge --extract EOPC_pos.txt --out EOPC --make-bed --recode A --keep ../EIGENSOFT/retain_individual.txt
plink --bfile ../merge --extract General_pos.txt --out LOPC --make-bed --recode A --keep ../EIGENSOFT/retain_individual.txt