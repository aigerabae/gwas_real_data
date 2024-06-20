This file contains further processing the resulting ped and map files; in this analysis, sex was imputed for some patients; patients and SNPs with high missing rates were excluded; MAFs were calculated

1) ped map to binary
plink --file kaz --make-bed --out kaz1
# Warning: 532649 het. haploid genotypes present (see kaz1.hh ); many commands
# treat these as missing.
# Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
# treat these as missing.
This data file has wrongfully assigned phenotypes

2) let's view X chromosome inbreeding (homozygosity) estimate F, plot it, and then impute sex
plink --bfile kaz1 --check-sex
Rscript --no-save gender_check.R

plink --bfile kaz1 --impute-sex 0.2 0.4 --make-bed --out kaz2
awk '$5 == "PROBLEM" {print $1, $2}' kaz2.sexcheck > problem_individuals.txt
plink --bfile kaz2 --remove problem_individuals.txt --make-bed --out kaz3

4) let's remove all non-autosomal regions
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' kaz3.bim > snp_1_22.txt
plink --bfile kaz3 --extract snp_1_22.txt --make-bed --out kaz4

5) remove missing
plink --bfile kaz4 --geno 0.02 --make-bed --out kaz5
plink --bfile kaz5 --mind 0.02 --make-bed --out kaz6

6) remove low MAFs; create table of MAFs
plink --bfile kaz7 --maf 0.05 --make-bed --out kaz8
plink  --bfile kaz7  --freq --out maf_kaz7

7) crytic relatedness
plink --bfile kaz11 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
plink --bfile kaz11 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
awk 'NR==FNR{a[$1,$2];next} (($1,$2) in a) || (($3,$4) in a)' related_pairs.txt missing_report.imiss | sort -k3,3n | awk '!seen[$1]++ {print $1, $2 > "0.2_low_call_rate_pihat.txt"}'
plink --bfile kaz11 --remove 0.2_low_call_rate_pihat.txt --make-bed --out kaz12

8) 
