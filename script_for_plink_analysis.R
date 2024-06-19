# 13 june 2024
# script for command line analyiss of PLINK default data
plink --file HB00001157 --make-bed --out HB00001157 
plink --bfile HB00001157 --geno 0.02 --make-bed --out HB00001157_1
plink --bfile HB00001157_1 --mind 0.02 --make-bed --out HB00001157_2
plink --bfile HB00001157_2 --check-sex 

# stopped here because turns out no one has sex recorded! or phenotypes
# i might go with the original raw data because this is a mess...
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
plink --bfile HB00001157_2 --remove sex_discrepancy.txt --make-bed --out HB00001157_3 



# 14 june 2024
# I decided to not use .idat files because it would be computationally intensive 
# and that part seems to be done fairly well. 
#Instead, i decided to generate my own .ped. and .fam files. 
# (because their .ped file didn't have SNP names or phenotypes
# and .map file had no information of location and chromosomes)
# To do that, I found the annotation file (to get info about SNPs locations),
# and gtReport.txt file that has information on genotypes in an appropriate format.
# I am also expecting a phenotype file sent separately


#            1) generating a .map file
# i found annotation file in .txt format in Analysis_Result.html -> 
# -> II-3. GSA MG v2 information i-> "1	Annotation File(GRCh37)" 
# and reformatted it into .map file
cut -f3,2,6 GSAMG2_SNP_Info_GRCh38.txt | awk 'BEGIN{OFS="\t"} {print $1, $2, "0", $3}' > output.map






# actually this might have been completely useless!!!! 
#skip until big new line!!
#!!!
#  !!!
#           2) generating a .ped file
# removed the first 10 lines, split the grReoirt.txt file into 8 smaller files
sed '1,10d' gtReport.txt > tmp.txt
split -l 100000 tmp.txt gtReport_chunk_
rm tmp.txt
# then i manually renamed them into gtReport1 to gtReport8

# removed spaces and replaced them with tabs, then transposed it 
# Function to process a single file
process_file() {
  file=$1
  awk '{printf $1"\t"; for(i=2;i<=NF;i++){printf $i"\t"}; print ""}' "$file" > "${file}_converted.tsv"
  awk -F'\t' 'NF' "${file}_converted.tsv" > "${file}_converted_cleaned.tsv"
  datamash transpose < "${file}_converted_cleaned.tsv" > "${file}_transposed.tsv"
}
# Export the function to make it available to GNU Parallel
export -f process_file
# Process the files in parallel
parallel -j 8 process_file ::: gtReport{1..8}

# then replaced '-' with '0' in 8 files
for i in {1..8}; do sed -i '' 's/-/0/g' getReport${i}; done
# concatenated the 8 files together
paste gtReport{1..8}_transposed.tsv > concatenated_transposed.tsv
# now need FID, IID, PID, MID, Sex, Phenotype;
# out of them most needed = IID, sex, phenotype
# the ready .map file is waiting for me in /Users/aygera/biostar/gwas/real_data/raw_data/GSAMG2_SNP_Info/output.map



# here is what i did with plink ped file:
# removed unnecesary space sand rpelaced them with tabs
awk -F'\t' '{gsub(/[[:space:]]+/,"\t"); print}' HB00001157.ped > fixed_result.ped




# 17 june 2024
# got my phenotypes file
# 1) phenotypes.xlsx -> phenotypes3.tsv
# 1.1) xlsx to csv
xlsx2csv phenotypes.xlsx phenotypes1.csv

# 1.2) csv to tsv
tr ',' '\t' < phenotypes1.csv > phenotypes2.tsv

# checked file
awk '{print $1}' phenotypes2.tsv | head

#1.3) removed header line 
sed 1d phenotypes2.tsv > phenotypes3.tsv
# saved header line as header_phenotypes.tsv
head -n 1 phenotypes2.tsv > header_phenotypes.tsv

#2)  GSAMG2_SNP_Info_GRCh38.txt -> kaz2.map
#tried removing header line from map file 
# and accidentally deleted the map file; 
# redid the process of creating it
cut -f3,2,6 GSAMG2_SNP_Info_GRCh38.txt | awk 'BEGIN{OFS="\t"} {print $1, $2, "0", $3}' > kaz.map

#2.1) realized there are too many empty fields for location in kaz.map
# how many fields are empty for column 1,2,3,4 in kaz.map?
awk 'BEGIN {count1=0; count2=0; count3=0; count4=0} {if ($1=="") count1++; if ($2=="") count2++; if ($3=="") count3++; if ($4=="") count4++} END {print "Column 1:", count1, "empty fields; Column 2:", count2, "empty fields; Column 3:", count3, "empty fields; Column 4:", count4, "empty fields."}' kaz.map
# Column 1: 0 empty fields; Column 2: 0 empty fields; Column 3: 0 empty fields; Column 4: 117857 empty fields.

# how many empty fields in column 4 and 6 GSAMG2_SNP_Info_GRCh38.txt?
awk -F'\t' '$4 ~ /^ *$/ {count++} END {print "Number of missing fields in field 4:", count}' GSAMG2_SNP_Info_GRCh38.txt
awk -F'\t' '$6 ~ /^ *$/ {count++} END {print "Number of missing fields in field 6:", count}' GSAMG2_SNP_Info_GRCh38.txt

# realized i need to use column 4 not 6 for location because 6 lacks many values
# 2.2) created new kaz1.map
cut -f3,2,4 GSAMG2_SNP_Info_GRCh38.txt | awk 'BEGIN{OFS="\t"} {print $1, $2, "0", $3}' > kaz1.map

# how many fields are empty for column 1,2,3,4 in kaz1.map?
awk 'BEGIN {count1=0; count2=0; count3=0; count4=0} {if ($1=="") count1++; if ($2=="") count2++; if ($3=="") count3++; if ($4=="") count4++} END {print "Column 1:", count1, "empty fields; Column 2:", count2, "empty fields; Column 3:", count3, "empty fields; Column 4:", count4, "empty fields."}' kaz1.map

# 2.3) delete headers line in map file 
sed 1d kaz1.map > kaz2.map


# 3) HB00001157.ped -> kaz?.ped
#now back to .ped file
# 3.1) created ped file again
awk -F'\t' '{gsub(/[[:space:]]+/,"\t"); print}' HB00001157.ped > kaz.ped

# 3.2) sorting phenotypes file and kaz.ped so that both are ordered the same according to ID column (2)
sort -t$'\t' -k2n,2 kaz.ped > kaz1.ped
sort -t$'\t' -k2n,2 phenotypes3.tsv > phenotypes4.tsv

  # viewing the result
awk '{print $1}' kaz1.ped | head
awk '{print $2}' phenotypes4.tsv | head

# view if they have the same number of rows
cat kaz1.ped | wc -l
cat phenotypes4.tsv | wc -l
# they don't! i need to remove those rows from phenotypes for which we don't have data

# 3.4) only included those tows into phenotypes5.tsv that have counterparts in kaz1.ped
awk 'FNR==NR{a[$2]++; next} ($2 in a)' kaz1.ped phenotypes4.tsv > phenotypes5.tsv

# viewed resulting columns
awk '{print $1, $2}' phenotypes5.tsv | head

# checks if 1 value in 1 table matches 1 value in another
awk -F'\t' 'FNR==NR{a[$2]++; next} {b[$2]++} END{for (i in a) if (a[i]!=1) print "Value", i, "in phenotypes5.tsv does not have exactly one match in kaz1.ped"; for (j in b) if (b[j]!=1) print "Value", j, "in kaz1.ped matches more than one value in phenotypes5.tsv"}' kaz1.ped phenotypes5.tsv

# 3.5) removing duplicates in each table
# saving separately the rows that have more than 1 matching
awk '$2 == 1314010' kaz1.ped > diff_1314010_kaz1.ped
awk '$2 == 1408910' kaz1.ped > diff_1408910_kaz1.ped

# viewing
awk 'NR==1 {split($0, line1, "\t")} NR==2 {split($0, line2, "\t")} END {for (i=1; i<=length(line1); i++) {if (line1[i] != line2[i]) print "Column " i ": Line 1 -", line1[i], "Line 2 -", line2[i]}}' diff_1314010_kaz1.ped
awk 'NR==1 {split($0, line1, "\t")} NR==2 {split($0, line2, "\t")} END {for (i=1; i<=length(line1); i++) {if (line1[i] != line2[i]) print "Column " i ": Line 1 -", line1[i], "Line 2 -", line2[i]}}' diff_1408910_kaz1.ped
cat header_phenotypes.tsv 
awk '$2 == 1314010' phenotypes5.tsv 
awk '$2 == 1408910' phenotypes5.tsv 
# so we have 4 patients that are different but have the same indices; 

# 3.6) excluding all 4 of them
awk -F'\t' '$2 != "1314010" && $2 != "1408910"' kaz1.ped > kaz2.ped
awk -F'\t' '$2 != "1314010" && $2 != "1408910"' phenotypes5.tsv > phenotypes6.tsv

# check for match again
awk -F'\t' 'FNR==NR{a[$2]++; next} {b[$2]++} END{for (i in a) if (a[i]!=1) print "Value", i, "in phenotypes6.tsv does not have exactly one match in kaz2.ped"; for (j in b) if (b[j]!=1) print "Value", j, "in kaz2.ped matches more than one value in phenotypes6.tsv"}' kaz2.ped phenotypes6.tsv
cat kaz2.ped | wc -l
cat phenotypes6.tsv | wc -l
# now i can match them together
# need to choose the factor to use
# possible options: BMI (use 25 as cutting value for obesity)

# tuberulosis
cut -f 11 phenotypes6.tsv | sort | uniq -c
# too few patients! there has to be at least 100 patients to have decent power

#for some reason the formatting broke and content of some columns appeared in others
# when i opened the same file in excel it was fine
# i saved it as phenotypes7.tsv and now it works fine!
cut -f 10 phenotypes7.tsv | sort | uniq -c
cut -f 11 phenotypes7.tsv | sort | uniq -c
cut -f 12 phenotypes7.tsv | sort | uniq -c
cut -f 13 phenotypes7.tsv | sort | uniq -c
cut -f 18 phenotypes7.tsv | sort | uniq -c
cut -f 19 phenotypes7.tsv | sort | uniq -c
cut -f 20 phenotypes7.tsv | sort | uniq -c
# yep -works fine now!

# counted number of numerical values for BMI
awk -F'\t' '$17 ~ /^[0-9]+(\.[0-9]+)?$/ {count++} END {print count}' phenotypes7.tsv

# patients with BMI < 18.5 
cut -f 17 phenotypes7.tsv | awk '$1 < 18.5 && $1 != ""' | wc -l
# patients with BMI < 24.9 
cut -f 17 phenotypes7.tsv | awk '$1 > 24.9 && $1 != ""' | wc -l

# removed non-kazakh patients
awk -F'\t' '$18 == "kazakh"' phenotypes7.tsv > phenotypes8.tsv
awk 'FNR==NR{a[$2]++; next} ($2 in a)' phenotypes8.tsv kaz2.ped > kaz3.ped

# check for matches
awk -F'\t' 'FNR==NR{a[$2]++; next} {b[$2]++} END{for (i in a) if (a[i]!=1) print "Value", i, "in phenotypes8.tsv does not have exactly one match in kaz3.ped"; for (j in b) if (b[j]!=1) print "Value", j, "in kaz3.ped matches more than one value in phenotypes8.tsv"}' kaz3.ped phenotypes8.tsv
cat kaz3.ped | wc -l
cat phenotypes8.tsv | wc -l
# perfect!

# I create a table phenotypes9.tsv that has numeric representations of data in phenotype8.tsv. 
# 1) I used Sample ID from column 2
# 2) I  used sex data in column 13 where Sex (1=male; 2=female; other=unknown)
# 3) I used column 17 of phenotypes8.tsv to divide patients into
# overweight 2 (BMI>25) 
# vs normal 1 (BMI<24.99999)
# vs. underweight 3 (<18.5)
# vs not available (-9)
cat phenotypes8.tsv | cut -f 2,13,17 > phenotypes9.tsv
awk -F"\t" '{ if ($2 == "male") $2 = 1; else if ($2 == "female") $2 = 2; if (!($3 ~ /^[0-9]*(\.[0-9]*)?$/)) $3 = -9; else { if ($3 < 18.5) $3 = 3; else if ($3 >= 18.5 && $3 <= 24.9999) $3 = 1; else if ($3 > 25) $3 = 2; } print $0; }' phenotypes9.tsv > phenotypes10.tsv
sed 's/ /\t/g' phenotypes10.tsv > phenotypes11.tsv

# create file kaz4.ped that has columns 1 to 4, 7 to 1532448 from kaz3.ped
# but column 5 from column 2 of phenotypes11.tsv,
# column 6 from column 3 of phenotype11.tsv
# checking number of columns in kaz3.ped
awk '{print NF; exit}' kaz3.ped

cut -f 1-4 kaz3.ped > kaz4_1.ped
cut -f 2,3 phenotypes11.tsv > phenotypes12.tsv
paste kaz4_1.ped phenotypes12.tsv > kaz4_2.ped
cut -f 7-1532448 kaz3.ped > kaz4_3.ped
paste kaz4_2.ped kaz4_3.ped > kaz5.ped

# need to change order of columns in map file 
awk '{print NF; exit}' kaz2.map
awk -F"\t" '{print $2 "\t" $1 "\t" $3 "\t" $4}' kaz2.map > kaz3.map