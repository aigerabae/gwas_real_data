# yay! quality control now!
cd ../qc
$ cp ../prepared_ped_map/kaz5.ped ./kaz.ped
$ cp ../prepared_ped_map/kaz3.map ./kaz.map

# 1) map/ped to binary
plink --file kaz --make-bed --out kaz1
# Warning: 532649 het. haploid genotypes present (see kaz1.hh ); many commands
#treat these as missing.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
#treat these as missing.

# 2) checking sex
# unnecesary - already have XY region
plink --bfile kaz1 --split-x b38 --make-bed --out kaz2

# is the problem from X or Y?
plink --bfile kaz1 --chr 24 --freq
plink --bfile kaz1 --chr 23 --freq
# the problem is partly on Y chromsome but mostly on X

# let's impute sex and change the threshold for calling males
plink --bfile kaz1 --impute-sex 0.2 0.4 --make-bed --out kaz2
awk '$5 == "PROBLEM" {print $1, $2}' kaz2.sexcheck > problem_individuals.txt
plink --bfile kaz2 --remove problem_individuals.txt --make-bed --out kaz3

# let's remove all non-autosomal regions
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' kaz3.bim > snp_1_22.txt
plink --bfile kaz3 --extract snp_1_22.txt --make-bed --out kaz4

# 2) remove missing
plink --bfile kaz4 --geno 0.02 --make-bed --out kaz5
plink --bfile kaz5 --mind 0.02 --make-bed --out kaz6

# 2.5) remove underweight
awk '{ if ($6 >= 3) print $1, $2}' kaz6.fam > underweight.txt
plink --bfile kaz6 --remove underweight.txt --make-bed --out kaz7

# 3) MAF
plink --bfile kaz7 --maf 0.05 --make-bed --out kaz8

# 4) HW
plink --bfile kaz8 --hwe 1e-6 --make-bed --out kaz9
plink --bfile kaz9 --hwe 1e-10 --hwe-all --make-bed --out kaz10

# 5) LD pruning
plink --bfile kaz10 --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile kaz10 --extract indepSNP.prune.in --het --out R_check
Rscript --no-save heterozygosity_outliers_list.R
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
plink --bfile kaz10 --remove het_fail_ind.txt --make-bed --out kaz11

# 6) crytic relatedness
plink --bfile kaz11 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
plink --bfile kaz11 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
awk 'NR==FNR{a[$1,$2];next} (($1,$2) in a) || (($3,$4) in a)' related_pairs.txt missing_report.imiss | sort -k3,3n | awk '!seen[$1]++ {print $1, $2 > "0.2_low_call_rate_pihat.txt"}'
plink --bfile kaz11 --remove 0.2_low_call_rate_pihat.txt --make-bed --out kaz12

# 7) PCA
plink2 --bfile kaz12 --pca 10 --out qcvcf 

#8 ) Association testing
plink --bfile kaz12 --covar qcvcf.eigenvec --logistic --hide-covar --thread-num 8 --out simple_logistic
awk '!/'NA'/' logistic_results.assoc.logistic > logistic_results.assoc.logistic1

plink --bfile kaz12 --covar qcvcf.eigenvec --logistic --dominant --hide-covar --thread-num 8 --out dominant_results
awk '!/'NA'/' dominant_results.assoc.logistic > dominant_results.assoc.logistic1

plink --bfile kaz12 --covar qcvcf.eigenvec --logistic --recessive --hide-covar --thread-num 8 --out recessive_results
awk '!/'NA'/' recessive_results.assoc.logistic > recessive_results.assoc.logistic1
