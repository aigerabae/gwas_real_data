This file contains further processing the resulting ped and map files; in this analysis, sex was imputed for some patients; patients and SNPs with high missing rates were excluded; MAFs were calculated

1) ped map to binary

```bash
plink --file kaz --make-bed --out kaz1
```
"Warning: 52344 het. haploid genotypes present (see kaz1.hh ); many commands treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands treat these as missing"

This warning means this data file has wrongfully assigned phenotypes (or low quality reads)

- create histogram of calling rate and remove individuals with low calling rate using phenotypes.tsv
```bash
awk 'FNR==NR {fam[$1]; next} $2 in fam {print $2, $5}' kaz1.fam phenotypes.tsv > calling_rate.txt
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
    if len(sys.argv) != 2:
        print("Usage: ./script_name.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    # Read the data from the file
    data = pd.read_csv(input_file, delim_whitespace=True, header=None, names=['Column1', 'Column2'])

    # Create a histogram of the values in column 2
    plt.figure(figsize=(10, 6))
    plt.hist(data['Column2'], bins=50, edgecolor='black')
    plt.xlabel('Samples')
    plt.ylabel('Frequency')
    plt.title('Histogram of Calling Rate of Samples')

    # Save the histogram as a PDF
    plt.savefig('calling_rate_hist.pdf')
    plt.close()

if __name__ == "__main__":
    main()
```

```bash
chmod +x create_histogram.py
./create_histogram.py calling_rate.txt
```

Viewing the histogram - there are some low calling individuals; let's remove them from out binary kaz files
```bash
awk '$2 < 0.9 {print $1,$1}' calling_rate.txt > low_calling_rate.txt
plink --bfile kaz1 --remove low_calling_rate.txt --make-bed --out kaz2
```

2) let's view X chromosome inbreeding (homozygosity) estimate F, plot it, and then impute sex

```bash
plink --bfile kaz2 --check-sex
Rscript --no-save gender_check.R
```

If we use ycount option we can see that females (as seen by their X chromosome F) have non-zero Y chromosome count (around 1000) while males have it around 6000

Let's impute sex for those 7 who were misgendered
```bash
plink --bfile kaz2 --impute-sex --make-bed --out kaz3
```

3) remove missing
```bash
plink --bfile kaz3 --geno 0.02 --make-bed --out kaz4
plink --bfile kaz4 --mind 0.02 --make-bed --out kaz5
```

4) remove low MAFs
```bash
plink --bfile kaz5 --maf 0.001 --make-bed --out kaz6
```

5) crytic relatedness
```bash
plink --bfile kaz6 --genome --min 0.2 --out pihat_min0.2
plink --bfile kaz6 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
```
I manually sorted through that file to keep as many individuals with the lowest missingness scores as possible while removing relatives; I put relatives that should be removed in 0.2_low_call_rate_pihat.txt: 

```bash
nano 0.2_low_call_rate_pihat.txt
```

```bash
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
```

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

6) remove non-acgt nucleotides
```bash
plink --bfile kaz7 --snps-only 'just-acgt' --make-bed --out kaz8
```

7) figuring out names of SNPs:

Let's convert them to standard names using Illumina Infinium Global Screening Array v2.0 Loci Name to rsID Conversion File (https://support.illumina.com/downloads/infinium-global-screening-array-v2-0-support-files.html) and exclude those SNPs that couldn't be converted to standard names; but first let's remove duplicates and names with more than 1 rsID (separated by comma)

```bash
awk -F'\t' '!seen[$1]++' GSA-24v2-0_A1_b150_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > GSA-dictionary.txt
```

```bash
plink --bfile kaz8 --update-name GSA-dictionary.txt --make-bed --out kaz9
awk '$2 !~ /^rs/' kaz9.bim | sort -k2,2 > non_rs_SNP.txt
plink --bfile kaz9 --exclude non_rs_SNP.txt --make-bed --out kaz10
```

checking if any non standard names remain
```bash
awk '$2 !~ /^rs/ {print}' kaz10.bim | sort -k2,2 
```
Nope! all clean

let's remove any duplicates and return it in plink1.9 version
```bash
plink2 --bfile kaz10 --rm-dup exclude-all --make-bed --out kaz11
plink --bfile kaz11 --make-bed --out kaz12
```

8) let's remove all non-autosomal regions
```bash
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' kaz12.bim > snp_1_22.txt
plink --bfile kaz12 --extract snp_1_22.txt --make-bed --out kaz12_autosomal
```

let's take out Y chromosome and mitochondrial SNPs; keep in mind that X=23,Y=24,XY=25,MT=26:
```bash
awk '{ if ($1 == 26) print $2 }' kaz12.bim > snp_mitoch.txt
plink --bfile kaz12 --extract snp_mitoch.txt --make-bed --out kaz12_mitoch
awk '{ if ($1 == 24) print $2 }' kaz12.bim > snp_y.txt
plink --bfile kaz12 --extract snp_y.txt --make-bed --out kaz12_y_chr
```

9) create table of MAFs
```bash
plink2 --bfile kaz12  --freq --out maf_kaz12
plink2 --bfile kaz12_autosomal  --freq --out maf_kaz12_autosomal
plink2 --bfile kaz12_mitoch  --freq --out maf_kaz12_mitoch
plink2 --bfile kaz12_y_chr --freq --out maf_kaz12_y_chr
```

10) plink binary to vcf
 ```bash
plink --bfile kaz12_autosomal --recode vcf --out kaz12_autosomal
plink --bfile kaz12_mitoch --recode vcf --out kaz12_mitoch
plink --bfile kaz12_y_chr --recode vcf --out kaz12_y_chr
plink --bfile kaz12 --recode vcf --out kaz12_all
```

11) adding additional info to vcf file (MAF and allele count)
```bash
bcftools view -h kaz12_autosomal.vcf > kaz_a1.vcf
bcftools view -H kaz12_autosomal.vcf > kaz_a2.vcf
cat maf_kaz12_autosomal.afreq | tail -n +2 | cut -f 6,7 > added_info.txt
cat maf_kaz12_autosomal.afreq | head -n 1 | cut -f 6,7 > added_header.txt
(cat kaz_a1.vcf | sed '$d'; paste <(tail -n 1 kaz_a1.vcf) added_header.txt) > kaz_a4.vcf
paste kaz_a2.vcf added_info.txt > kaz_a3.vcf
cat kaz_a4.vcf kaz_a3.vcf > kaz_a5.vcf
```

now the header has bcftools in it.. need to remove? 

12) runs of homozygosity (ROH)

plink --bfile kaz12_autosomal --homozyg-density 60 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-window-het 0

make HGDP into plink binary file
Build 36.1 to 38

Ped file
a) getting p1 from phenotypes metadata, encoding male and female as 1 and 2, sorting
```bash
awk '{print $1"\t" $1"\t" "0\t" "0\t" $2"\t" "-9\t"}' metadata.txt | tail -n +2 > HGDP_p1_1.txt
awk -F'\t' '{ if ($5 == "male") $5 = 1; else if ($5 == "female") $5 = 2; print }' OFS='\t' HGDP_p1_1.txt > HGDP_p1_2.txt
sort -k1,1 HGDP_p1_2.txt > HGDP_p1_3.txt
```

b) Transposing for p2 (genotypes)
```bash
split -l 90000 HGDP1.txt file_
```

I renamed them to gtReport1 to gtReport8:
```bash
mv file_aa ./gtReport1
mv file_ab ./gtReport2
mv file_ac ./gtReport3
mv file_ad ./gtReport4
mv file_ae ./gtReport5
mv file_af ./gtReport6
mv file_ag ./gtReport7
mv file_ah ./gtReport8
```

then using the code below removed spaces and replaced them with tabs, then transposed it, sorted the file, removed empty lines from genotypes file, removed all intermediate files
(takes about 1 hour to run)

```bash
process_file() {
  file=$1
  awk '{printf $1"\t"; for(i=2;i<=NF;i++){printf $i"\t"}; print ""}' "$file" > "${file}_converted.tsv"
  awk -F'\t' 'NF' "${file}_converted.tsv" > "${file}_converted_cleaned.tsv"
  datamash transpose < "${file}_converted_cleaned.tsv" > "${file}_transposed.tsv"
}
export -f process_file
parallel -j 8 process_file ::: gtReport{1..8}
for i in {1..8}; do sed -i '' 's/-/0/g' gtReport${i}_transposed.tsv; done
paste gtReport{1..8}_transposed.tsv > HGDP_transposed.txt
sort -k1,1 HGDP_transposed.txt | tail -n +2 > HGDP_p2_1.txt
sed -E 's/([^\t]+)\t/\1\t/;s/\t(.)/\t\1\t/g' HGDP_p2_1.txt > HGDP_p2_2.txt
rm *gtRe*
```

c) joining p1 and p2
```bash
join -t$'\t' HGDP_p1_3.txt HGDP_p2_2.txt | sed 's/\t\t*/\t/g' > HGDP.ped
```

d) Map file
```bash
awk '{print $2"\t" $1 "\t0\t" $3}' HGDP_Map.txt > HGDP.map
```

e) Ped/map to plink binary
```bash
plink --file HGDP --missing-code -9,0,NA,na,- --make-bed --out HGDP
```

f) remap 36.1 to 38 build (for some reason cut wasnt working properly with chromosomes; i dont know why)
```bash
comm -12 <(awk '{print $1}' dictionary_pos | sort) <(awk '{print $2}' HGDP.bim | sort) > common_snps.txt
cat kaz12_autosomal.bim | awk '{print $2"\t" $1}' > dictionary_chr
cat kaz12_autosomal.bim | cut -f 2,4 > dictionary_pos
plink --bfile HGDP --extract common_snps.txt --make-bed --out HGDP1
plink --bfile HGDP1 --update-chr dictionary_chr --make-bed --out HGDP2
plink --bfile HGDP2 --update-map dictionary_pos --make-bed --out HGDP3
```

?
sort -k4,4n -t$'\t' HGDP3.bim > HGDP4.bim
cp HGDP3.bed ./HGDP4.bed
cp HGDP3.fam ./HGDP4.fam
?

Problem! Need to update the binary bed file using info in bed/fam files? Or just keep it unsorted? 

Now need to do proper QC of these data and then do ROH and Fst analysis and also extract their ethicities and red PCA...

use all populations from HDGP and SGDP for ROH and Fst and only selected eurasian populations for PCA

