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

# Figure out what to do with 8623 indels that you have in the hg38 dataset:
cat kaz1.bim | awk '$5 == "I "|| $5 == "D" || $6 == "I" || $6 == "D"' | wc -l

# Figure out where can you do quality control on a dataset with indels. also you might not be able to use the binary files, may have to go back to the ped map files.

# also try running genome studio on  hg19 and see how many SNPs and indels you have and if they are the same or not
```



