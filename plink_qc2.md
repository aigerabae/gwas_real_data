This file contains further processing the resulting ped and map files; in this analysis, sex was imputed for some patients; patients and SNPs with high missing rates were excluded; MAFs were calculated

1) ped map to binary

```bash
plink --file kaz --make-bed --out kaz
```
"Warning: 532649 het. haploid genotypes present (see kaz1.hh ); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing.
This warning means this data file has wrongfully assigned phenotypes"

- create histogram of calling rate and remove individuals with low calling rate using phenotypes.tsv
```bash
awk 'FNR==NR {fam[$1]; next} $2 in fam {print $2, $5}' kaz.fam phenotypes.tsv > calling_rate.txt
```

Let's write a script that would plot the histogram for us in a pdf file
```bash
nano create_histogram.py 
```

```python
#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import sys

def main():
    # Read the data from the file
    data = pd.read_csv('calling_rate.txt', delim_whitespace=True, header=None, names=['Column1', 'Column2'])

    # Create a histogram of the values in column 1
    plt.figure(figsize=(10, 6))
    plt.hist(data['Column2'], bins=50, edgecolor='black')
    plt.xlabel('Samples')
    plt.ylabel('Frequency')
    plt.title('Histogram of Calling Rate of Samples')

    # Save the histogram as a PDF
    plt.savefig('calling_rate_hist.pdf')
    plt.close()
```

```bash
chmod +x create_histogram.py
./create_histogram.py
```

Viewing the histogram - there are some low calling individuals; let's remove them from out binary kaz files
awk '$2 < 0.9 {print $1,$1}' calling_rate.txt > low_calling_rate.txt
plink --bfile kaz --remove low_calling_rate.txt --make-bed --out kaz1


2) let's view X chromosome inbreeding (homozygosity) estimate F, plot it, and then impute sex

```bash
plink --bfile kaz1 --check-sex
Rscript --no-save gender_check.R
```

3) some individuals are clearly just misgendered and 3 are unclear (intermediate values); let's get misgendered individuals removed; the rest should have their sex assigned correctly

```bash
plink --bfile kaz1 --impute-sex --make-bed --out kaz2
```

```bash
plink --bfile kaz2 --check-sex --out kaz2
awk '$5 == "PROBLEM" {print $1, $2}' kaz2.sexcheck > problem_individuals.txt
plink --bfile kaz2 --remove problem_individuals.txt --make-bed --out kaz3
```

4) remove missing
```bash
plink --bfile kaz3 --geno 0.02 --make-bed --out kaz4
plink --bfile kaz4 --mind 0.02 --make-bed --out kaz5
```
5) remove low MAFs
plink --bfile kaz5 --maf 0.001 --make-bed --out kaz6

6) crytic relatedness
```bash
plink --bfile kaz6 --genome --min 0.2 --out pihat_min0.2
plink --bfile kaz6 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
```
I manually sorted through that file to keep as many individuals with the lowest missingness scores as possible while removing relatives; I put relatives that should be removed in 0.2_low_call_rate_pihat.txt: 
5	5
6	6
8	8
16	16
20	20
22	22
25	25
73	73
81	81
82	82
100	100
104	104
105	105
107	107
108	108
117	117
122	122
125	125
126	126
142	142
WE019	WE019
WE058	WE058
WE070	WE070
WE091	WE091
WE092	WE092

search for F value: 

```bash
id=WE002
# example ID
awk -v id="$id" '$2 == id {print $6}' missing_report.imiss
# prints its respective missingness rate
```

I removed relatives from that list I created manually
```bash
plink --bfile kaz6 --remove 0.2_low_call_rate_pihat.txt --make-bed --out kaz7
```

7) remove non-acgt nucleotides
```bash
plink --bfile kaz7 --snps-only 'just-acgt' --make-bed --out kaz8
```

8) figuring out names of SNPs:

Let's convert them to standard names using Illumina Infinium Global Screening Array v2.0 Loci Name to rsID Conversion File (https://support.illumina.com/downloads/infinium-global-screening-array-v2-0-support-files.html) and exclude those SNPs that couldn't be converted to standard names

```bash
plink --bfile kaz8 --update-name GSA-24v2-0_A1_b150_rsids.txt --make-bed --out kaz9
awk '$2 !~ /^rs/ {print}' kaz9.bim | sort -k2,2 > non_rs_SNP.txt
plink --bfile kaz9 --exclude non_rs_SNP.txt --make-bed --out kaz10
```

9) let's remove all non-autosomal regions
```bash
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' kaz10.bim > snp_1_22.txt
plink --bfile kaz10 --extract snp_1_22.txt --make-bed --out kaz10_autosomal
```

let's take out Y chromosome and mitochondrial SNPs; keep in mind that X=23,Y=24,XY=25,MT=26:
```bash
awk '{ if ($1 == 26) print $2 }' kaz10.bim > snp_mitoch.txt
plink --bfile kaz10 --extract snp_mitoch.txt --make-bed --out kaz10_mitoch
awk '{ if ($1 == 24) print $2 }' kaz10.bim > snp_y.txt
plink --bfile kaz10 --extract snp_y.txt --make-bed --out kaz10_y_chr
```
Warning: 1 het. haploid genotype present (see kaz9_y_chr.hh ); many commands
treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.

10) create table of MAFs
```bash
plink  --bfile kaz10  --freq --out maf_kaz10
plink  --bfile kaz10_autosomal  --freq --out maf_kaz10_autosomal
```

11) plink binary to vcf
 ```bash
plink --bfile kaz10_autosomal --recode vcf
plink --bfile kaz10_mitoch --recode vcf
plink --bfile kaz10_y_chr --recode vcf
```
