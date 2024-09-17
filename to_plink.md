Affyemtrix data:
Remove the first 5 rows of ah.txt and save as ah1.txt
remove rows where column 836 is empty in ah1.txt and save as ah2.txt
remove rows where column 835 doesn't equal 2 and save as ah3.txt
for every row in columns 2-805 inclusive: replace "A" with content of coplumn 819 from that same row; replace "B" with content of coplumn 820 from that same row; replace "NoCall" with "0" and save as ah4.txt
in row 1 of ah4.txt columns 2 to 805 inclusive remove in each cell everything after "_(" (inclusive); also remove "AH" from cells; save as ah5.txt
save in ah5.txt columns 807, 836, 0s, 808 (use awk in that exact order) as ah.map
save columns 2 to 805 inclusive as genotypes.txt
transpose genotypes.txt and save as genotypes2.txt
save head_ah.tsv second sheet (id_sex_pheno) from personal google sheets account
in id_sex_pheno.tsv remove rows where columns 2 and 3 have empty values and save that as id_sex_pheno1.tsv
save column 1 of genotypes2.txt as samples.tsv
use join to get a table with samples.tsv joined by columns 2 and 3 from id_sex_pheno1.tsv; save as id_sex_pheno2.tsv
save ah1.ped with column 1 = all 0s; column 2 = column 1 of genotypes2.txt; column 3 all 0s; column 4 = all 0s; column 5 = column 2 of id_sex_pheno2.tsv; column 6 = column 3 of id_sex_pheno2.tsv; columsn 7 to end = columns 2 to end of genotypes2.txt

```bash
# Step 1: Remove the first 5 rows from ah.txt and save as ah1.txt
tail -n +6 ah.txt > ah1.txt

# Step 2: Remove rows where column 836 is empty from ah1.txt and save as ah2.txt
awk -F'\t' '$836 != ""' ah1.txt > ah2.txt

# Step 3: Remove rows where column 835 doesn't equal 2 from ah2.txt and save as ah3.txt (except row 1)
awk -F'\t' '$835 == 2' ah2.txt > ah3.txt

Checking what we have in cells for genotypes:
awk -F'\t' '{for(i=2; i<=805; i++) print $i}' ah3.txt | sort | uniq -c

# Step 4: For columns 2-805, replace AA with column 819 twice merged; AB with columns 819 and 820 merged; BB with column 820 twice merged; A with column 819 merged with 0; B with column 820 merged with 0 and "NoCall", NoCall_1" and "ZeroCN" with "0", save as ah4.txt
awk -F'\t' -v OFS='\t' '{for(i=2; i<=805; i++) {
    if($i == "AA") $i = $819 $819;
    else if($i == "AB") $i = $819 $820;
    else if($i == "BB") $i = $820 $820;
    else if($i == "A") $i = $819 "00";
    else if($i == "B") $i = $820 "00";
    else if($i == "NoCall" || $i == "NoCall_1" || $i == "ZeroCN") $i = "00";
} print}' ah3.txt > ah4.txt

Checking what we have in cells for genotypes:
awk -F'\t' '{for(i=2; i<=805; i++) print $i}' ah4.txt | sort | uniq -c

Removing all non binary genotypes (long stretches of nucleotides identified earelier)
awk -F'\t' -v OFS='\t' '{keep = 1; for(i=2; i<=805; i++) {if(length($i) > 2) {keep = 0; break;}} if(keep) print}' ah4.txt > ah5.txt
sed -i 's/-/0/g' ah5.txt

Checking what we have in cells for genotypes:
awk -F'\t' '{for(i=2; i<=805; i++) print $i}' ah5.txt | sort | uniq -c

# Step 5: Modify the header row in ah1.txt and concatenate with processed data without header
awk -F'\t' -v OFS='\t' 'NR == 1 {for(i=2; i<=805; i++) {sub(/_\(.*$/, "", $i); sub(/AH/, "", $i);} print}' ah1.txt > ah1_header.txt
cat ah1_header.txt ah5.txt > ah6.txt

# Step 6: Save columns 807, 836, 0s, 808 in that order to ah.map
awk -F'\t' -v OFS='\t'  '{ print $807, $836, "0", $808 }' ah6.txt > ah.map

# Step 7: Save columns 2-805 to genotypes.txt
cut -f 2-805 ah6.txt > genotypes.txt
sed -i 's/ //g' genotypes.txt 

# Step 8: Transpose genotypes.txt and save as genotypes2.txt
datamash transpose <genotypes.txt > genotypes2.txt

# Step 9: In id_sex_pheno.tsv, remove rows where columns 2 and 3 are empty and save as id_sex_pheno1.tsv
awk '$2 != "" && $3 != ""' id_sex_pheno.tsv > id_sex_pheno1.tsv

# Step 10: Save column 1 of genotypes2.txt as samples.tsv
awk '{ print $1 }' genotypes2.txt > samples.tsv

# Step 11: Use join to get a table with samples.tsv joined by columns 2 and 3 from id_sex_pheno1.tsv and save as id_sex_pheno2.tsv
sort -k1,1 samples.tsv > samples1.tsv
sort -k1,1 id_sex_pheno1.tsv > id_sex_pheno2.tsv
join -1 1 -2 1 samples1.tsv <(awk '{print $1, $2, $3}' id_sex_pheno2.tsv) > id_sex_pheno3.tsv
sed -i 's/ \+/\t/g' id_sex_pheno3.tsv

# Step 12: Save ah1.ped
paste <(awk '{print "0"}' genotypes2.txt) <(cut -f 1 genotypes2.txt) <(awk '{print "0"}' genotypes2.txt) <(awk '{print "0"}' genotypes2.txt) <(cut -f 2 id_sex_pheno3.tsv) <(cut -f 3 id_sex_pheno2.tsv) > temp_first_6_columns.txt
cut -f 2- genotypes2.txt > temp_genotypes_columns.txt
paste temp_first_6_columns.txt temp_genotypes_columns.txt > ah1.ped








FinalReport to plink::
in ad.txt remove the first 10 rows and save as ad1.txt
In ad1.txt merge columns 3 and 4, so that their merged content is column 3 and save as ad2.txt
in ad2.txt use pivot_table to turn long data to wide with index=column 2, columns = column 1, values = column 3; save as pivoted.txt
save ad.ped with 0s in columns 1,3,4,5,6; for column 2 use column 1 from pivoted.txt; for columns 2 to end use columns 7 to end of pivoted.txt

gwas:
cat gwas300_FinalReport.txt | cut -f 1,2,3,4 > gwas300.txt
awk 'NR > 1 { print $3, $2, 0, $4 }' SNP_Map.txt > kaz.map
dont forget to remove snps with 0 as chromosome!

in gwas300.txt remove the first 10 rows and save as kaz1.txt
tail -n +11 "gwas300.txt" > kaz1.txt

In kaz1.txt merge columns 3 and 4, so that their merged content is column 3 and save as kaz2.txt
awk 'BEGIN {OFS="\t"} {$3=$3""$4; $4=""; print}' kaz1.txt > kaz2.txt

pivoting
awk '{printf "%s%s", $3, (NR%665608 ? OFS : ORS)}' kaz2.txt > pivoted1.txt
cat kaz2.txt | cut -f 2 | uniq > pivoted_header.txt

save ad.ped with 0 in columns 1,3,4,5,6; for column 2 use pivoted_header.txt; for columns 2 to end use all columns from pivoted1.txt
awk 'BEGIN { OFS="\t" } { print 0, $0, 0, 0, 0, 0 }' pivoted_header.txt > temp_columns.txt
paste temp_columns.txt pivoted1.txt > kaz.ped

plink --file kaz --make-bed --out kaz
awk -F'\t' '!seen[$1]++' GSA-24v2-0_A1_b150_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > GSA-dictionary.txt
plink --bfile kaz --update-name GSA-dictionary.txt --make-bed --out kaz1
awk '$2 !~ /^rs/' kaz1.bim | sort -k2,2 > non_rs_SNP.txt
plink --bfile kaz1 --exclude non_rs_SNP.txt --make-bed --out kaz2
awk '$1 == 0 {print $2}' kaz2.bim > exclude_snps.txt
plink --bfile kaz2 --exclude exclude_snps.txt --make-bed --out kaz3

merge kaz3 and gnomad binary data using positions on chromosomes; use names from kaz3; save merged dataset as merged binary dataset
... Unsuccessful. Only 20000 common SNPs...

Alzheimer idats to plink
awk 'NR > 1 { print $3, $2, 0, $4 }' SNP_Map.txt > alz.map
tail -n +11 "alz_FinalReport_2.txt" > alz1.txt
awk 'BEGIN {OFS="\t"} {$3=$3""$4; $4=""; print}' alz1.txt > alz2.txt
awk '{printf "%s%s", $3, (NR%253702 ? OFS : ORS)}' alz2.txt > pivoted_alz.txt

For some reason it doesn't pivot it properly... it generates or at least displays only

cat kaz2.txt | cut -f 2 | uniq > pivoted_header.txt
awk 'BEGIN { OFS="\t" } { print 0, $0, 0, 0, 0, 0 }' pivoted_header.txt > temp_columns.txt
paste temp_columns.txt pivoted_alz.txt > alz.ped




plink --file alz --make-bed --out alz
awk -F'\t' '!seen[$1]++' InfiniumImmunoArray-24v2-0_A_b138_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > IA-dictionary.txt
plink --bfile alz --update-name IA-dictionary.txt --make-bed --out alz1
awk '$2 !~ /^rs/' alz1.bim | sort -k2,2 > non_rs_SNP.txt
plink --bfile alz1 --exclude non_rs_SNP.txt --make-bed --out alz2
awk '$1 == 0 {print $2}' alz2.bim > exclude_snps.txt
plink --bfile alz2 --exclude exclude_snps.txt --make-bed --out alz3
