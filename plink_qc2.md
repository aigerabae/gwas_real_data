This file contains further processing the resulting ped and map files; in this analysis, sex was imputed for some patients; patients and SNPs with high missing rates were excluded; MAFs were calculated

1) ped map to binary

```bash
plink --file kaz --make-bed --out kaz1
Warning: 532649 het. haploid genotypes present (see kaz1.hh ); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.
This warning means this data file has wrongfully assigned phenotypes
```

2) let's view X chromosome inbreeding (homozygosity) estimate F, plot it, and then impute sex
```bash
plink --bfile kaz1 --check-sex
Rscript --no-save gender_check.R
```

some individuals are clearly just misgendered and 3 are unclear (intermediate values)
Let's impute those who are just assigned wrong sex
```bash
plink --bfile kaz1 --impute-sex --make-bed --out kaz2
```

```bash
plink --bfile kaz2 --check-sex --out kaz2
awk '$5 == "PROBLEM" {print $1, $2}' kaz2.sexcheck > problem_individuals.txt
plink --bfile kaz2 --remove problem_individuals.txt --make-bed --out kaz3
```

3) individuals removed; the rest should have their sex assigned correctly


4) let's remove all non-autosomal regions
```bash
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' kaz3.bim > snp_1_22.txt
plink --bfile kaz3 --extract snp_1_22.txt --make-bed --out kaz4
```

5) remove missing
```bash
plink --bfile kaz4 --geno 0.02 --make-bed --out kaz5
plink --bfile kaz5 --mind 0.02 --make-bed --out kaz6
```

6) remove low MAFs; create table of MAFs
```bash
plink --bfile kaz6 --maf 0.01 --make-bed --out kaz7
plink  --bfile kaz7  --freq --out maf_kaz7
```

7) crytic relatedness
```bash
plink --bfile kaz7 --genome --min 0.2 --out pihat_min0.2
plink --bfile kaz7 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
```

I manually sorted through that file to keep as many individuals with the lowest missingness scores as possible while removing relatives; I put relatives that should be removed in 0.2_low_call_rate_pihat.txt
search for F value: 

```bash
id=WE002 
awk -v id="$id" '$2 == id {print $6}' missing_report.imiss
```

I removed relatives from that list I created manually
```bash
plink --bfile kaz7 --remove 0.2_low_call_rate_pihat.txt --make-bed --out kaz8
```

8) remove non-acgt nucleotides; plink binary to vcf
 ```bash
plink --bfile kaz8 --snps-only 'just-acgt' --make-bed --out kaz9
plink --bfile kaz9 --recode vcf
```
