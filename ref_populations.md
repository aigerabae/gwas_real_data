HGDP: accessed at https://www.hagsc.org/hgdp/files.html 

Metadata: accessed at https://www.internationalgenome.org/data-portal/data-collection/hgdp 

1) HGDP text file 36.1 build -> plink binary file with 38 build

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

3) Only keeping selected ethicities from HGDP
```bash
cut -f 1,5 metadata.txt > ethnicities1.txt
awk 'NR==FNR {ids[$1]; next} $1 in ids {print $1"\t" $2}' HGDP8.fam ethnicities.txt > ethnicities2.txt
cat ethnicities2.txt <(awk '{print $1 "\tKazakh"}' kaz12_autosomal.fam) > ethnicities3.txt
cat ethnicities3.txt | grep -e "Kazakh" -e "Uygur" -e "Hazara" -e "Russian" -e "French" -e "Basque" -e "Bergamo" -e "Pathan" -e "Sindhi" -e "Kalash" -e "Adygei" -e "Bedouin" -e "Mozabite" -e "Japanese" -e "Northern" -e "Mongolian" -e "Yakut" -e "Han" | awk '{print $1"\t" $2}' > ethnicities4.txt
cat ethnicities4.txt | awk '{print $1"\t" $1}' > selected_ethnicities.txt
plink --bfile HGDP8 --keep selected_ethnicities.txt --biallelic-only strict --make-bed --out HGDP9
```

4) Merging kazakh and HGDP data 
First - deal with multiallelic variants
```bash
plink --bfile kaz12_autosomal --bmerge HGDP9.bed HGDP9.bim HGDP9.fam --make-bed --out merged1
plink --bfile HGDP9 --exclude merged1-merge.missnp --biallelic-only strict --make-bed --out HGDP10
plink --bfile kaz12_autosomal --exclude merged1-merge.missnp --biallelic-only strict --make-bed --out kaz13_autosomal
plink --bfile kaz13_autosomal --bmerge HGDP10.bed HGDP10.bim HGDP10.fam --make-bed --out merged2
plink --bfile merged2 --geno 0.02 --make-bed --out merged3
plink --bfile merged3 --mind 0.02 --make-bed --out merged4
echo -e "rs9267522\nrs11229\nrs75412302\nrs12660719" > duplicates.txt
plink --bfile merged4 --exclude duplicates.txt --make-bed --out merged5
```

Remove outliers that show on PCA:
```bash
nano outliers_kazakh.txt
```
```bash
1210510 1210510
1302810  1302810
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
plink --bfile merged5 --remove outliers.txt --make-bed --out merged6
```

5) Downloading turkic,siberian,caucasus,jewish (has uzbek) data from Estonian Biocentre:
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

6) Merging estonian data together:
```bash
plink --bfile caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand --bmerge jew_paper_data_dbSNP-b131_pos-b37_1KG_strand --make-bed --out all1
plink --bfile all1 --bmerge sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand --make-bed --out all2
plink --bfile all2 --bmerge turkic --make-bed --out all3
```

7) Keeping only selected ethnicities in estonian dataset (not including kazakh)
I made metadata.txt in google sheets that contains 3 columns: ID, ethnicity, region
```bash
cat metadata.txt | awk '{print $1 "\t" $1}' > ethnic.txt 
plink --bfile all3 --keep ethnic.txt --make-bed --out all4
```

e) Changing build from 37 to 38 and merging datasets and their metadata:
Use ethnicities6.txt and merged6 binary files for HGDP merged dataset

```bash
comm -12 <(awk '{print $2}' merged6.bim | sort) <(awk '{print $2}' all4.bim | sort) > common_snps.txt
cat merged6.bim | awk '{print $2"\t" $1}' > dictionary_chr
cat merged6.bim | cut -f 2,4 > dictionary_pos
plink --bfile all4 --extract common_snps.txt --make-bed --out all5
plink2 --bfile all5 --update-chr dictionary_chr --update-map dictionary_pos --sort-vars --make-pgen --out all6
plink2 --pfile all6 --make-bed --out all7
plink --bfile merged6 --bmerge all7 --make-bed --out all8
# this line above is for making a missnp file

plink --bfile all7 --exclude all8-merge.missnp --make-bed --out all8
plink --bfile merged6 --exclude all8-merge.missnp --make-bed --out merged7
```

9) Merging filtered estonian data with my kazakh + HGDP data
```bash
plink --bfile merged7 --bmerge all8 --make-bed --out all9
cat metadata.txt metadata_kaz.txt > ethnic1.txt
```

10) QC of Kazakh + HGDP + Estonian data
```bash
plink --bfile all9 --geno 0.02 --make-bed --out all10
plink --bfile all10 --mind 0.02 --make-bed --out all11
plink --bfile all11 --maf 0.001 --make-bed --out all12
plink --bfile all12 --genome --min 0.2 --out pihat_min0.2
plink --bfile all12 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
plink --bfile all12 --snps-only 'just-acgt' --make-bed --out all13
```

Updating metadata to have the same individuals as in fam file:
```bash
awk '{print $1}' all13.fam | grep -Fwf - ethnic1.txt > ethnic2.txt
```

11) PCA:
```bash
plink2 --bfile all13 --pca 10 --out all_pca
python plot_eigenvec.py all_pca.eigenvec ethnic2.txt
```

12) Fst
```bash
cat ethnic2.txt | awk '{print $1 "\t" $1 "\t" $2}' > ethnic2_fst.txt
plink2 --bfile all13 --fst CATPHENO --within ethnic2_fst.txt --double-id --out fst_output
./plot_fst_heatmap.py fst_output.fst.summary
```

13) IBD for kazakhs with all populations
plink --bfile all14 --genome --out ibd_all

14) ROH
plink --bfile kaz12_autosomal --homozyg-density 60 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-window-het 0

15) ADMIXTURE 
Pruning data and adjusting metadata:
```bash
plink --bfile all13 --indep-pairwise 1000 150 0.4 --out pruned_data
plink --bfile all13 --extract pruned_data.prune.in --make-bed --out all14
awk '{print $1}' all14.fam | grep -Fwf - ethnic2.txt > ethnic3.txt
```

ADMIXTURE
```bash
for K in 8; do admixture --cv all14.bed -j8 $K | tee log${K}.out; done
for K in 10; do admixture --cv all14.bed -j8 $K | tee log${K}.out; done
for K in 3 5 15; do admixture --cv all14.bed -j8 $K | tee log${K}.out; done
```

Plotting:
Average stacked plot:
```bash
python average_plot_admixture.py all14.8.Q ethnic3.txt
```

Regular plot:
```bash
python safe_plot_admixture.py all14.8.Q ethnic3.txt
```

Round plot with kazakhs in the middle:
```bash
ln -s ~/tools/AncestryPainter_v5/AncestryPainter.pl ./
cat ethnic3.txt | awk '{print $2"\t" $1}'  > ethnic3.ind
perl AncestryPainter.pl -i ethnic3.ind -q ./all14.8.Q -t Kazakh -o Kazakh -f png
```

Making a table with average percentages for each population:
```bash
python create_admixtures_table.py all14.8.Q ethnic3.txt 
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
perl AncestryPainter.pl -i ethnic5_5.ind -q ./all15_5.5.Q -t Kazakh -o Kazakh -f png
python safe_plot_admixture.py all15_5.5.Q ethnic_5.txt

awk '{print $1}' all15_5.fam | grep -Fwf - ethnic_5.txt > ethnic4_5.txt
cat ethnic4_5.txt | awk '{print $2"\t" $1}'  > ethnic5_5.ind
perl AncestryPainter.pl -i ethnic5_5.ind -q ./all15_5.8.Q -t Kazakh -o Kazakh -f png
python safe_plot_admixture.py all15_5.8.Q ethnic_5.txt

awk '{print $1}' all15_8.fam | grep -Fwf - ethnic_8.txt > ethnic4_8.txt
cat ethnic4_8.txt | awk '{print $2"\t" $1}'  > ethnic8_5.ind
perl AncestryPainter.pl -i ethnic8_5.ind -q ./all15_8.5.Q -t Kazakh -o Kazakh_8_5  -f png
python safe_plot_admixture.py all15_8.5.Q ethnic_8.txt

awk '{print $1}' all15_8.fam | grep -Fwf - ethnic_8.txt > ethnic4_8.txt
cat ethnic4_8.txt | awk '{print $2"\t" $1}'  > ethnic8_5.ind
perl AncestryPainter.pl -i ethnic8_5.ind -q ./all15_8.8.Q -t Kazakh -o Kazakh_8_8  -f png
python safe_plot_admixture.py all15_8.8.Q ethnic_8.txt
```
