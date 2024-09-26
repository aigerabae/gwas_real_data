Affymetrix data:
```bash
# Step 1: Remove the first 5 rows from ah.txt and save as ah1.txt
tail -n +6 ah.txt > ah1.txt

# Step 2: Remove rows where column 836 is empty from ah1.txt and save as ah2.txt
awk -F'\t' '$836 != ""' ah1.txt > ah2.txt

# Step 3: Remove rows where column 835 doesn't equal 2 from ah2.txt and save as ah3.txt (except row 1)
awk -F'\t' '$835 == 2' ah2.txt > ah3.txt

# Step 4: For columns 2-805, replace AA with column 819 twice merged; AB with columns 819 and 820 merged; BB with column 820 twice merged; A with column 819 merged with 0; B with column 820 merged with 0 and "NoCall", NoCall_1" and "ZeroCN" with "0", save as ah4.txt
awk -F'\t' -v OFS='\t' '{for(i=2; i<=805; i++) {
    if($i == "AA") $i = $819 $819;
    else if($i == "AB") $i = $819 $820;
    else if($i == "BB") $i = $820 $820;
    else if($i == "A") $i = "00";
    else if($i == "B") $i = "00";
    else if($i == "NoCall" || $i == "NoCall_1" || $i == "ZeroCN") $i = "00";
} print}' ah3.txt > ah4.txt

# Step 5: Removing all non binary genotypes (long stretches of nucleotides identified earelier)
awk -F'\t' -v OFS='\t' '{keep = 1; for(i=2; i<=805; i++) {if(length($i) > 2) {keep = 0; break;}} if(keep) print}' ah4.txt > ah5.txt
sed -i 's/-/0/g' ah5.txt

# Step 6: Modify the header row in ah1.txt and concatenate with processed data without header
awk -F'\t' -v OFS='\t' 'NR == 1 {for(i=2; i<=805; i++) {sub(/_\(.*$/, "", $i); sub(/AH/, "", $i);} print}' ah1.txt > ah1_header.txt
cat ah1_header.txt ah5.txt > ah6.txt

# Step 7: Save columns 807, 836, 0s, 808 in that order to ah.map
awk -F'\t' -v OFS='\t' 'NR > 1 { print $807, $836, "0", $808 }' ah6.txt > ah.map

# Step 7: Save columns 2-805 to genotypes.txt
cut -f 2-805 ah6.txt > genotypes.txt
sed -i 's/ //g' genotypes.txt 

# Step 8: Transpose genotypes.txt and save as genotypes2.txt; sort genotypes2.txt
datamash transpose <genotypes.txt > genotypes2.txt
sort -k1,1 genotypes2.txt > genotypes3.txt

# Step 9: In id_sex_pheno.tsv, remove rows where column 2 and 3 are empty and save as id_sex_pheno1.tsv
awk '$2 != "" && $3 != ""' id_sex_pheno.tsv > id_sex_pheno1.tsv

# Step 11: sort id_sex_pheno1.tsv and fix formatting
sort -k1,1 id_sex_pheno1.tsv > id_sex_pheno2.tsv
sed -i 's/ \+/\t/g' id_sex_pheno2.tsv

# Step 12: Save ah1.ped
join -1 1 -2 1 id_sex_pheno2.tsv genotypes3.txt | awk '{printf "0\t%s\t0\t0\t", $1; for (i=2; i<=NF; i++) printf "%s\t", $i; printf "\n"}' > combined1.txt
sed -i 's/\r//g' combined1.txt
cat combined1.txt | cut -f 7- > genotypes4.txt
sed 's/./&\t/g' genotypes4.txt > genotypes5.txt
sed 's/\t\+/\t/g' genotypes5.txt > genotypes6.txt
cut -f1-6 combined1.txt | paste - genotypes6.txt > ah.ped

plink --file ah --missing-code -9,0,NA,na --make-bed --out ah

It! doesn't! Work!3
Possibly irregular .ped line.  Restarting scan, assuming multichar alleles.
Rescanning .ped file... 0%
Error: Half-missing call in .ped file at variant 5119, line 1.
(base) 
```







FinalReport to plink::
in ad.txt remove the first 10 rows and save as ad1.txt
In ad1.txt merge columns 3 and 4, so that their merged content is column 3 and save as ad2.txt
in ad2.txt use pivot_table to turn long data to wide with index=column 2, columns = column 1, values = column 3; save as pivoted.txt
save ad.ped with 0s in columns 1,3,4,5,6; for column 2 use column 1 from pivoted.txt; for columns 2 to end use columns 7 to end of pivoted.txt

gwas:
```bash
cat gwas300_FinalReport.txt | cut -f 1,2,3,4 > gwas300.txt
awk 'NR > 1 { print $3, $2, 0, $4 }' SNP_Map.txt > kaz.map
# dont forget to remove snps with 0 as chromosome!

# in gwas300.txt remove the first 10 rows and save as kaz1.txt
tail -n +11 "gwas300.txt" > kaz1.txt

# In kaz1.txt merge columns 3 and 4, so that their merged content is column 3 and save as kaz2.txt
awk 'BEGIN {OFS="\t"} {$3=$3""$4; $4=""; print}' kaz1.txt > kaz2.txt

# pivoting
awk '{printf "%s%s", $3, (NR%665608 ? OFS : ORS)}' kaz2.txt > pivoted1.txt

# save ad.ped with 0 in columns 1,3,4,5,6; for column 2 use pivoted_header.txt; for columns 2 to end use all columns from pivoted1.txt
cat kaz2.txt | cut -f 2 | uniq > pivoted_header.txt
awk 'BEGIN { OFS="\t" } { print 0, $0, 0, 0, 0, 0 }' pivoted_header.txt > temp_columns.txt
paste temp_columns.txt pivoted1.txt > kaz.ped

# Merging
plink --file kaz --make-bed --out kaz
awk -F'\t' '!seen[$1]++' GSA-24v2-0_A1_b150_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > GSA-dictionary.txt
plink --bfile kaz --update-name GSA-dictionary.txt --make-bed --out kaz1
awk '$2 !~ /^rs/' kaz1.bim | sort -k2,2 > non_rs_SNP.txt
plink --bfile kaz1 --exclude non_rs_SNP.txt --make-bed --out kaz2
awk '$1 == 0 {print $2}' kaz2.bim > exclude_snps.txt
plink --bfile kaz2 --exclude exclude_snps.txt --make-bed --out kaz3

# merge kaz3 and gnomad binary data using positions on chromosomes; use names from kaz3; save merged dataset as merged binary dataset (script available in gnomad.md)
# ... Unsuccessful. Only 20000 common SNPs... admixture still looks bad
```

Alzheimer idats to plink
```bash
awk 'NR > 1 { print $3, $2, 0, $4 }' SNP_Map.txt > alz.map
tail -n +11 "alz_FinalReport_2.txt" > alz1.txt
awk 'BEGIN {OFS="\t"} {$3=$3""$4; $4=""; print}' alz1.txt > alz2.txt
sed $'s/\r/\t/g' alz2.txt > alz3.txt
awk '{printf "%s%s", $3, (NR%253702 ? OFS : ORS)}' alz3.txt > pivoted_alz.txt
cat alz3.txt | cut -f 2 | uniq > pivoted_header.txt
awk 'BEGIN { OFS="\t" } { print 0, $0, 0, 0, 0, 0 }' pivoted_header.txt > temp_columns.txt
paste temp_columns.txt pivoted_alz.txt > alz.ped
```

Changing snp names to standard rsID names, excluding SNPs on 0 chromosome, adding sex and phenotypes:
```bash
plink --file alz --make-bed --out alz
awk -F'\t' '!seen[$1]++' InfiniumImmunoArray-24v2-0_A_b138_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > IA-dictionary.txt
plink --bfile alz --update-name IA-dictionary.txt --make-bed --out alz1
awk '$2 !~ /^rs/' alz1.bim | sort -k2,2 > non_rs_SNP.txt
plink --bfile alz1 --exclude non_rs_SNP.txt --make-bed --out alz2
awk '$1 == 0 {print $2}' alz2.bim > exclude_snps.txt
plink --bfile alz2 --exclude exclude_snps.txt --make-bed --out alz3
```

I used metadata (from different zapusks) to identify samples IDs and their matching Sentrix_ID positions; 
i wanted to manually replace 207851060016 with 207859430016 since they were the only samples that faild to replace names and they had an identical numebr of samples so i figured it could be a misspelling... but then i decided to ust keep them as they are since even if i did i would have 13 more samples with no phenotype data so i just excluded 29 samples with no phenotype data (among them 24 have no sample ID and 5 had a sample ID but were not present in the metadata).

I didn't remove other ethnicities 

```bash
plink --bfile alz3 --update-ids 'alzheimer_metadata_selected - true_dict.tsv' --make-bed --out alz4
plink --bfile alz4 --update-sex 'alzheimer_metadata_selected - sex.tsv' --make-bed --out alz5
plink --bfile alz5 --pheno 'alzheimer_metadata_selected - pheno.tsv' --make-bed --out alz6
awk '$6 == -9' alz6.fam | awk '{print $1"\t" $2}' > missing_phenotype.tsv
plink --bfile alz6 --remove missing_phenotype.tsv --allow-no-sex --make-bed --out alz7
plink --bfile alz7 --geno 0.02 --make-bed --out alz8
plink --bfile alz8 --mind 0.02 --make-bed --out alz9
plink --bfile alz9 --maf 0.001 --make-bed --out alz10
plink --bfile alz10 --genome --min 0.2 --out pihat_min0.2
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
plink --bfile alz10 --missing --out missing_report
echo "0    D180
0    C020
0    AK015
0    C132
0    C124
0    C077
0    C067
0    C170" > relatives_to_remove.tsv
plink --bfile alz10 --remove relatives_to_remove.tsv --allow-no-sex --make-bed --out alz11
plink --bfile alz11 --genome --min 0.2 --out 2pihat_min0.2
awk '$10 > 0.2 {print $1, $2, $3, $4}' 2pihat_min0.2.genome > 2related_pairs.txt
plink2 --bfile alz11 --pca 10 --out alz_pca
```

```bash
plink --bfile alz11 --covar alz_pca.eigenvec --logistic --hide-covar --thread-num 8 --out simple_logistic
plink --bfile alz11 --covar alz_pca.eigenvec --logistic --dominant --hide-covar --out dominant_results
plink --bfile alz11 --covar alz_pca.eigenvec --logistic --recessive --hide-covar --out recessive_results
plink --bfile alz11 --allow-no-sex --assoc --out assoc_results
```

```bash
cat dominant_results.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat recessive_results.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat simple_logistic.assoc.logistic | awk '$9 != "NA"' | sort -gk 9,9 | head 
cat assoc_results.assoc | awk '$9 != "NA"' | sort -gk 9,9 | head 
```
No association so far! MAX = 10e-7, 10-20 10e-5
