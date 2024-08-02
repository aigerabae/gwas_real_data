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
nano outliers.txt
```

```bash
1210510 1210510
HGDP00621       HGDP00621
HGDP01270       HGDP01270
HGDP01271       HGDP01271
HGDP00175       HGDP00175
HGDP00953       HGDP00953
HGDP00949       HGDP00949
HGDP00969       HGDP00969
HGDP00959       HGDP00959
HGDP00621       HGDP00621
1302810  1302810
```

```bash
plink --bfile merged5 --remove outliers.txt --make-bed --out merged6
plink2 --bfile merged6 --pca 10 --out all_pca 
```

```bash
python plot_eigenvec.py all_pca.eigenvec ethnicities4.txt
```

6) runs of homozygosity (ROH) and Fst
```bash
plink --bfile merged6 --homozyg-density 60 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-window-het 0

cut -f 1,5 metadata.txt > metadata1.txt
awk '{print $2, "Kazakh"}' kaz12_autosomal.fam | sort >> metadata1.txt
awk 'NR==FNR {ids[$1]; next} $1 in ids {print $1"\t" $2}' merged6.fam metadata1.txt > metadata2.txt
awk '{print $1"\t" $1"\t" $2}' metadata2.txt > metadata3.txt
plink2 --bfile merged6 --fst CATPHENO --within metadata3.txt --double-id --out fst_output
chmod +x plot_fst_heatmap.py
./plot_fst_heatmap.py fst_output.fst.summary
```

7) admixture
```bash
for K in 1 2 3 4 5 6 7 8 9 10; \
do admixture --cv merged6.bed $K | tee log${K}.out; done
grep -h CV log*.out
```

Find ALDH2 gene in kazakh and other populations and see whether we absorb alcohol better or rose than other central asians or europeans
More data to use for PCA:

Simons: can be accessed at https://www.simonsfoundation.org/simons-genome-diversity-project/ via cancer genomics cloud seven bridges metadata can be accessed at https://www.nature.com/articles/nature18964#Sec10 Supplementary Table 1

1000 genomes can be accessed at https://www.internationalgenome.org/data-portal/population
