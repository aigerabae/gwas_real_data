Data accessed at https://gnomad.broadinstitute.org/downloads#v4-core-dataset
with 
gsutil cp gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/f2_fst/hgdp_tgp.bim ./
gsutil cp gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/f2_fst/hgdp_tgp.bed ./
gsutil cp gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/f2_fst/hgdp_tgp.fam ./

1) remove families and get them into separate metadata file
cat hgdp_tgp.fam | cut -f 1,2 > metadata.txt
plink --bfile hgdp_tgp --recode vcf --out gnomad1
plink2 --vcf gnomad1.vcf --const-fid 0 --make-bed --out gnomad2


3) QC of HGDP data
```bash
plink --bfile hgdp_tgp --geno 0.02 --make-bed --out gnomad1
plink --bfile gnomad1 --mind 0.02 --make-bed --out gnomad2
plink --bfile gnomad2 --maf 0.001 --make-bed --out gnomad3
plink --bfile gnomad3 --snps-only 'just-acgt' --make-bed --out gnomad4
```  

2) PCA
Only keeping selected ethicities
```bash
cut -f 1,5 metadata.txt > ethnicities1.txt
awk 'NR==FNR {ids[$1]; next} $1 in ids {print $1"\t" $2}' HGDP8.fam ethnicities.txt > ethnicities2.txt
cat ethnicities2.txt <(awk '{print $1 "\tKazakh"}' kaz12_autosomal.fam) > ethnicities3.txt
cat ethnicities3.txt | grep -e "Kazakh" -e "Uygur" -e "Hazara" -e "Russian" -e "French" -e "Basque" -e "Bergamo" -e "Pathan" -e "Sindhi" -e "Kalash" -e "Adygei" -e "Bedouin" -e "Mozabite" -e "Japanese" -e "Northern" -e "Mongolian" -e "Yakut" -e "Han" | awk '{print $1"\t" $2}' > ethnicities4.txt
plink --bfile HGDP8 --keep selected_ethnicities.txt --biallelic-only strict --make-bed --out HGDP9
```

Merging kazakh and HGDP data (first - deal with multiallelic variants) and doing PCA;
```bash

plink --bfile kaz12_autosomal --bmerge HGDP9.bed HGDP9.bim HGDP9.fam --make-bed --out merged1
plink --bfile HGDP9 --exclude merged_dataset-merge.missnp --biallelic-only strict --make-bed --out HGDP10
plink --bfile kaz12_autosomal --exclude merged_dataset-merge.missnp --biallelic-only strict --make-bed --out kaz13_autosomal
plink --bfile kaz13_autosomal --bmerge HGDP10.bed HGDP10.bim HGDP10.fam --make-bed --out merged2
plink --bfile merged2 --geno 0.02 --make-bed --out merged3
plink --bfile merged3 --mind 0.02 --make-bed --out merged4
echo -e "rs9267522\nrs11229\nrs75412302\nrs12660719" > duplicates.txt
plink --bfile merged4 --exclude duplicates.txt --make-bed --out merged5
```

```bash
nano outliers_kazakh.txt
```

Attention: This includes potential outliers (starting with number 3) that I am thinking about excluding based on PCA but haven't decided yet
```bash
1210510 1210510
1302810  1302810
1413810  1413810  
WE016  WE016
WE128  WE128
211  211
88  88
1407010  1407010
2  2
1217910  1217910
113  113
33  33
240  240
```

```bash
nano outliers_HGDP.txt
```

```bash
HGDP00621       HGDP00621
HGDP01270       HGDP01270
HGDP01271       HGDP01271
HGDP00175       HGDP00175
HGDP00953       HGDP00953
HGDP00949       HGDP00949
HGDP00969       HGDP00969
HGDP00959       HGDP00959
```

```bash
cat outliers_kazakh.txt outliers_HGDP.txt > outliers.txt
```

```bash
plink --bfile merged5 --remove outliers.txt --make-bed --out merged6
plink2 --bfile merged6 --pca 10 --out all_pca 
```

```bash
python plot_eigenvec.py all_pca.eigenvec ethnicities4.txt
```

6) runs of homozygosity (ROH)
```bash
plink --bfile kaz12_autosomal --homozyg-density 60 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-window-het 0
```

7) Fst
```bash
cut -f 1,5 metadata.txt > metadata1.txt
awk '{print $2, "Kazakh"}' kaz12_autosomal.fam | sort >> metadata1.txt
awk 'NR==FNR {ids[$1]; next} $1 in ids {print $1"\t" $2}' merged6.fam metadata1.txt > metadata2.txt
awk '{print $1"\t" $1"\t" $2}' metadata2.txt > metadata3.txt
plink2 --bfile merged6 --fst CATPHENO --within metadata3.txt --double-id --out fst_output
chmod +x plot_fst_heatmap.py
./plot_fst_heatmap.py fst_output.fst.summary
```

8) admixture
First - LD pruning:
```bash
plink --bfile merged6 --indep-pairwise 1000 150 0.4 --out pruned_data
plink --bfile merged6 --extract pruned_data.prune.in --make-bed --out merged7

awk 'BEGIN {OFS="\t"} {if ($2 == "Kazakh" || $2 == "Hazara" || $2 == "Uygur") region = "Central_asia"; else if ($2 == "Bergamo Italia" || $2 == "French" || $2 == "Basque" || $2 == "Russian") region = "Europe"; else if ($2 == "Adygei") region = "Caucasus"; else if ($2 == "Pathan" || $2 == "Sindhi" || $2 == "Kalash") region = "South_asia"; else if ($2 == "Han" || $2 == "Northern" || $2 == "Japanese" || $2 == "Mongolian" || $2 == "Yakut") region = "East_asia"; else if ($2 == "Bedouin" || $2 == "Mozabite") region = "Middle_east"; else region = "Unknown"; print $1, $2, region;}' ethnicities4.txt > ethnicities5.txt
awk 'NR==FNR {a[$1]; next} !($1 in a)' outliers.txt ethnicities5.txt > ethnicities6.txt
```

```bash
for K in 5 6 7 8 9 10; \
do admixture --cv merged7.bed -j8 $K | tee log${K}.out; done
grep -h CV log*.out

admixture --cv merged6.bed -j8 12
grep -v -f <(awk '{print $1}' outliers.txt | sort | uniq) ethnicities5.txt > ethnicities6.txt
python safe_plot_admixture.py merged6.12.Q ethnicities6.txt
```

For plotting in R
```bash
awk 'NR==FNR {ethnicity[FNR]=$2; population[FNR]=$3; sampleID[FNR]=$1; next} {for (i=1; i<=NF; i++) print sampleID[FNR], ethnicity[FNR], population[FNR], $i, i}' ethnicities6.txt merged7.5.Q > equal_samples.tsv
```

Problem: ADMIXTURE plot doesn't look like there is much different between anyone
Potential solutin: adding data from other sources

a) Downloading turkic,siberian,caucasus,jewish (has uzbek) data from Estonian Biocentre:
```bash
wget https://evolbio.ut.ee/turkic/turkic.fam
wget https://evolbio.ut.ee/turkic/turkic.bim
wget https://evolbio.ut.ee/turkic/turkic.bed
wget https://evolbio.ut.ee/caucasus/caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand.bim
wget https://evolbio.ut.ee/caucasus/caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand.fam
wget https://evolbio.ut.ee/caucasus/caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand.bed
wget https://evolbio.ut.ee/sakha/sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand.bim
wget https://evolbio.ut.ee/sakha/sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand.fam
wget https://evolbio.ut.ee/sakha/sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand.bed
wget https://evolbio.ut.ee/jew/jew_paper_data_dbSNP-b131_pos-b37_1KG_strand.bim
wget https://evolbio.ut.ee/jew/jew_paper_data_dbSNP-b131_pos-b37_1KG_strand.fam
wget https://evolbio.ut.ee/jew/jew_paper_data_dbSNP-b131_pos-b37_1KG_strand.bed
```

b) Merging estonian data together:
```bash
plink --bfile caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand --bmerge jew_paper_data_dbSNP-b131_pos-b37_1KG_strand --make-bed --out all1
plink --bfile all1 --bmerge sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand --make-bed --out all2
plink --bfile all2 --bmerge turkic --make-bed --out all3
```

c) I made metadata.txt in excel that contains 3 columns: ID, ethnicity, region
```bash
cat metadata.txt | awk '{print $1 "\t" $1}' > ethnic.txt 
plink --bfile all3 --keep ethnic.txt --make-bed --out all4
```

d) Now I want to merge it with my kazakh + HGDP data; copying 
cp ../../p3/ethnicities6.txt ./metadata_kaz.txt
cp ../../p3/*merged6* ./

e) Changing build from 37 to 38 and merging datasets and their metadata:
```bash
comm -12 <(awk '{print $2}' merged6.bim | sort) <(awk '{print $2}' all4.bim | sort) > common_snps.txt
cat merged6.bim | awk '{print $2"\t" $1}' > dictionary_chr
cat merged6.bim | cut -f 2,4 > dictionary_pos
plink --bfile all4 --extract common_snps.txt --make-bed --out all5
plink2 --bfile all5 --update-chr dictionary_chr --update-map dictionary_pos --sort-vars --make-pgen --out all6
plink2 --pfile all6 --make-bed --out all7
plink --bfile all7 --exclude all8-merge.missnp --make-bed --out all8
plink --bfile merged6 --exclude all8-merge.missnp --make-bed --out merged7
plink --bfile merged7 --bmerge all8 --make-bed --out all9
cat metadata.txt metadata_kaz.txt > ethnic1.txt
```

f) QC:
```bash
plink --bfile all9 --geno 0.02 --make-bed --out all10
plink --bfile all10 --mind 0.02 --make-bed --out all11
plink --bfile all11 --maf 0.001 --make-bed --out all12
plink --bfile all12 --genome --min 0.2 --out pihat_min0.2
plink --bfile all12 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
plink --bfile all12 --snps-only 'just-acgt' --make-bed --out all13
```

Remove individuals that don't have a counterpart in fam file from metadata
```bash
awk '{print $1}' all13.fam | grep -Fwf - ethnic1.txt > ethnic2.txt
```

g) PCA:
```bash
plink2 --bfile all13 --pca 10 --out all_pca
python plot_eigenvec.py all_pca.eigenvec ethnic2.txt
```

h) ADMIXTURE - don't fotget to run 2,3,4 and 13,14,15 and plot them
```bash
plink --bfile all13 --indep-pairwise 1000 150 0.4 --out pruned_data
plink --bfile all13 --extract pruned_data.prune.in --make-bed --out all14
for K in 13 14 15; do admixture --cv all14.bed -j8 $K | tee log${K}.out; done
for K in 2 3 4; do admixture --cv all14.bed -j8 $K | tee log${K}.out; done

python safe_plot_admixture.py all14.5.Q ethnic2.txt
python safe_plot_admixture.py all14.6.Q ethnic2.txt
python safe_plot_admixture.py all14.7.Q ethnic2.txt
python safe_plot_admixture.py all14.8.Q ethnic2.txt
python safe_plot_admixture.py all14.9.Q ethnic2.txt
python safe_plot_admixture.py all14.10.Q ethnic2.txt
python safe_plot_admixture.py all14.11.Q ethnic2.txt
python safe_plot_admixture.py all14.12.Q ethnic2.txt
```

i) Fst
```bash
cat ethnic2.txt | awk '{print $1 "\t" $1 "\t" $2}' > ethnic3.txt
plink2 --bfile all14 --fst CATPHENO --within ethnic3.txt --double-id --out fst_output
chmod +x plot_fst_heatmap.py
./plot_fst_heatmap.py fst_output.fst.summary
```

Find ALDH2 gene in kazakh and other populations and see whether we absorb alcohol better or rose than other central asians or europeans
SNPs:

Alcohol-related SNPs (only 1 present in kazakh dataset and none in merged)
rs2018417
rs28626993
rs28913916
rs28914782
rs41275697
rs41275699
rs113075608
rs13306164
rs141556759
rs150631941
rs190914158
rs201108880
rs201582342



More data to use for PCA:

**Raw files:**
Simons: can be accessed at https://www.simonsfoundation.org/simons-genome-diversity-project/ via cancer genomics cloud seven bridges metadata can be accessed at https://www.nature.com/articles/nature18964#Sec10 Supplementary Table 1

1000 genomes can be accessed at https://www.internationalgenome.org/data-portal/population

**Geno data:**
10000 ancient and present genomes for qpAdm: https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data

**PLINK data**
HGDP + 1000 genomes in plink: https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg - has also data for HapMap
Simons in plink: https://reichdata.hms.harvard.edu/pub/datasets/sgdp/ - but it only has 1-3 people for each population


Getting alzheimer rsID:
```bash
plink --bfile kaz12 --extract rsID_table --make-bed --out AD_table 
plink --bfile kaz12 --extract rsID_study --make-bed --out AD_study 
plink2 --bfile AD_table --freq --out AD_table
plink2 --bfile AD_study --freq --out AD_study
```

Visualizing with AncestryPainter 
ln -s ~/tools/AncestryPainter_v5/AncestryPainter.pl ./
awk '{print $1}' all14.fam | grep -Fwf - ethnic2.txt > ethnic4.txt
cat ethnic4.txt | awk '{print $2"\t" $1}'  > ethnic5.ind
perl AncestryPainter.pl -i ethnic5.ind -q ./all14.8.Q -t Kazakh -o Kazakh -l nolines -f png
perl AncestryPainter.pl -i ethnic5.ind -q all14.8.Q -f png
