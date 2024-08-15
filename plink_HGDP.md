HGDP: accessed at https://www.hagsc.org/hgdp/files.html accessed metadata at https://www.internationalgenome.org/data-portal/data-collection/hgdp 

1) make HGDP into plink binary file with 38 build

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
comm -12 <(awk '{print $2}' kaz12_autosomal.bim | sort) <(awk '{print $2}' HGDP.bim | sort) > common_snps.txt
cat kaz12_autosomal.bim | awk '{print $2"\t" $1}' > dictionary_chr
cat kaz12_autosomal.bim | cut -f 2,4 > dictionary_pos
plink --bfile HGDP --extract common_snps.txt --make-bed --out HGDP1
plink2 --bfile HGDP1 --update-chr dictionary_chr --update-map dictionary_pos --sort-vars --make-pgen --out HGDP2
plink2 --pfile HGDP2 --make-bed --out HGDP3
```

2) QC of HGDP data
```bash
plink --bfile HGDP3 --geno 0.02 --make-bed --out HGDP4
plink --bfile HGDP4 --mind 0.02 --make-bed --out HGDP5
plink --bfile HGDP5 --maf 0.001 --make-bed --out HGDP6
plink --bfile HGDP6 --genome --min 0.2 --out pihat_min0.2
plink --bfile HGDP6 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
```

```bash
plink --bfile HGDP6 --remove to_remove.txt --make-bed --out HGDP7
plink --bfile HGDP7 --snps-only 'just-acgt' --make-bed --out HGDP8
```

4) PCA

Only keeping selected ethicities
```bash
cut -f 1,5 metadata.txt > ethnicities1.txt
awk 'NR==FNR {ids[$1]; next} $1 in ids {print $1"\t" $2}' HGDP8.fam ethnicities.txt > ethnicities2.txt
cat ethnicities2.txt <(awk '{print $1 "\tKazakh"}' kaz12_autosomal.fam) > ethnicities3.txt
cat ethnicities3.txt | grep -e "Kazakh" -e "Uygur" -e "Hazara" -e "Russian" -e "French" -e "Basque" -e "Bergamo" -e "Pathan" -e "Sindhi" -e "Kalash" -e "Adygei" -e "Bedouin" -e "Mozabite" -e "Japanese" -e "Northern" -e "Mongolian" -e "Yakut" -e "Han" | awk '{print $1"\t" $2}' > ethnicities4.txt
cat ethnicities4.txt | awk '{print $1"\t" $1}' > selected_ethnicities.txt
plink --bfile HGDP8 --keep selected_ethnicities.txt --biallelic-only strict --make-bed --out HGDP9
```

ADMIXTURE without kazakhs: still looks weird
```bash
cp ../hgdp/ethnicities3.txt ./
cat ethnicities3.txt | grep -e "Uygur" -e "Hazara" -e "Russian" -e "French" -e "Basque" -e "Bergamo" -e "Pathan" -e "Sindhi" -e "Kalash" -e "Adygei" -e "Bedouin" -e "Mozabite" -e "Japanese" -e "Northern" -e "Mongolian" -e "Yakut" -e "Han" -e "Yoruba" -e "Biaka" -e "Bantu" | awk '{print $1"\t" $2"\t" $2}' > ethnicities3_HGDP_only.txt
cat ethnicities3_HGDP_only.txt | awk '{print $1"\t" $1}' > selected_ethnicities_african_no_kazakh.txt
plink --bfile HGDP8 --keep selected_ethnicities_african_no_kazakh.txt --biallelic-only strict --make-bed --out HGDP8_african
plink --bfile HGDP8_african --indep-pairwise 1000 150 0.4 --out pruned_data
plink --bfile HGDP8_african --extract pruned_data.prune.in --make-bed --out HGDP8_pruned

admixture --cv HGDP8_pruned.bed -j8 8
python safe_plot_admixture.py HGDP8_pruned.8.Q ethnicities3_HGDP_only.txt
admixture --cv HGDP8_pruned.bed -j8 3
python safe_plot_admixture.py HGDP8_pruned.3.Q ethnicities3_HGDP_only.txt
python average_plot_admixture.py HGDP8_pruned.3.Q ethnicities3_HGDP_only.txt

cat ethnicities3_HGDP_only.txt | awk '{print $2"\t" $1}' > ethnicities3_HGDP.ind
perl AncestryPainter.pl -i ethnicities3_HGDP.ind -q ./all14.8.Q -f png

plink2 --bfile HGDP8_pruned --pca 10 --out all_pca
python plot_eigenvec.py all_pca.eigenvec ethnicities3_HGDP_only.txt
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

awk '{print $1}' all14.fam | grep -Fwf - ethnic2.txt > ethnic4.txt

python safe_plot_admixture.py all14.5.Q ethnic4.txt
python safe_plot_admixture.py all14.6.Q ethnic4.txt
python safe_plot_admixture.py all14.7.Q ethnic4.txt
python safe_plot_admixture.py all14.8.Q ethnic4.txt
python safe_plot_admixture.py all14.9.Q ethnic4.txt
python safe_plot_admixture.py all14.10.Q ethnic4.txt
python safe_plot_admixture.py all14.11.Q ethnic4.txt
python safe_plot_admixture.py all14.12.Q ethnic4.txt

python average_plot_admixture.py all14.5.Q ethnic4.txt

```

ADMIXTURE with 5/8 populations:
```bash
cat ethnic2.txt | grep -e "Kazakh" -e "Russian" -e "Bedouin" -e "Sindhi" -e "Japanese" -e "French" | awk '{print $1"\t" $1}' > extract_5.txt
cat ethnic2.txt | grep -e "Kazakh" -e "Russian" -e "Bedouin" -e "Sindhi" -e "Japanese" -e "French" | awk '{print $1"\t" $2"\t" $3}' > ethnic_5.txt
plink --bfile all14 --keep extract_5.txt --make-bed --out all15_5
for K in 5 8; do admixture --cv all15_5.bed -j8 $K | tee log${K}.out; done

cat ethnic2.txt | grep -e "Kazakh" -e "Russian" -e "Bedouin" -e "Sindhi" -e "Japanese" -e "French" -e "Mongolian" -e "Tajik" -e "Turks" | awk '{print $1"\t" $1}' > extract_8.txt
cat ethnic2.txt | grep -e "Kazakh" -e "Russian" -e "Bedouin" -e "Sindhi" -e "Japanese" -e "French" -e "Mongolian" -e "Tajik" -e "Turks" | awk '{print $1"\t" $2"\t" $3}' > ethnic_8.txt
plink --bfile all14 --keep extract_8.txt --make-bed --out all15_8
for K in 5 8; do admixture --cv all15_8.bed -j8 $K | tee log${K}.out; done

awk '{print $1}' all15_5.fam | grep -Fwf - ethnic_5.txt > ethnic4_5.txt
cat ethnic4_5.txt | awk '{print $2"\t" $1}'  > ethnic5_5.ind
perl AncestryPainter.pl -i ethnic5_5.ind -q ./all15_5.5.Q -t Kazakh -o Kazakh -l nolines -f png
python safe_plot_admixture.py all15_5.5.Q ethnic_5.txt

awk '{print $1}' all15_5.fam | grep -Fwf - ethnic_5.txt > ethnic4_5.txt
cat ethnic4_5.txt | awk '{print $2"\t" $1}'  > ethnic5_5.ind
perl AncestryPainter.pl -i ethnic5_5.ind -q ./all15_5.8.Q -t Kazakh -o Kazakh -l nolines -f png
python safe_plot_admixture.py all15_5.8.Q ethnic_5.txt

awk '{print $1}' all15_8.fam | grep -Fwf - ethnic_8.txt > ethnic4_8.txt
cat ethnic4_8.txt | awk '{print $2"\t" $1}'  > ethnic8_5.ind
perl AncestryPainter.pl -i ethnic8_5.ind -q ./all15_8.5.Q -t Kazakh -o Kazakh -l nolines -f png
python safe_plot_admixture.py all15_8.5.Q ethnic_8.txt

awk '{print $1}' all15_8.fam | grep -Fwf - ethnic_8.txt > ethnic4_8.txt
cat ethnic4_8.txt | awk '{print $2"\t" $1}'  > ethnic8_5.ind
perl AncestryPainter.pl -i ethnic8_5.ind -q ./all15_8.8.Q -t Kazakh -o Kazakh -l nolines -f png
python safe_plot_admixture.py all15_8.8.Q ethnic_8.txt
```

Visualizing with AncestryPainter
```bash
ln -s ~/tools/AncestryPainter_v5/AncestryPainter.pl ./
awk '{print $1}' all14.fam | grep -Fwf - ethnic2.txt > ethnic4.txt
cat ethnic4.txt | awk '{print $2"\t" $1}'  > ethnic5.ind
perl AncestryPainter.pl -i ethnic5.ind -q ./all14.8.Q -t Kazakh -o Kazakh -f png
```

**Ancient genomes:**
Changing phenotype to 1 to allow convertf (didn't test yet)
plink --bfile all14 --pheno dataset.fam --make-pheno 1 '*' --make-bed --out all15

Ancient genoems available at 
https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW
https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data

Tutorial:
https://indoaryan.com/qpadm-tutorial/

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
10000 ancient and present genomes for qpAdm: 
https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
Also at: https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW

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
