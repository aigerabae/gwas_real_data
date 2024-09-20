This script described work I did on annovar file from Aset:

Adding extended (KAZ_MAF + allele count) + changing ref/alt to kaz_ref/alt and DB (database) ref/alt:
```bash
cat autosomal_ext_for_annovar.vcf | tail -n +31 | cut -f 234,235 > for_extended_vcf.tsv
paste kaz_gwas_for_annovar_224.LATEST.annovar.hg19_multianno.header.txt for_extended_vcf.tsv > final_annovared_extended.tsv
sed -i '' -e 's/REF/kaz_ref/g' -e 's/ALT/kaz_alt/g' -e 's/Ref/DB_ref/g' -e 's/Alt/DB_alt/g' final_annovared_extended.tsv
```

Print how many missing values are in columns 1 to 452:
```bash
awk -F"\t" '{for(i=1;i<=452;i++) if($i == ".") count[i]++} END{for(i=1;i<=452;i++) print "Column " i ": " count[i]}' final_annovared_extended.tsv 
```

Print mutation type and exonic function (3 databases):
```bash
cat final_annovared_extended.tsv | cut -f 6 | sort | uniq -c
cat final_annovared_extended.tsv | cut -f 11 | sort | uniq -c
cat final_annovared_extended.tsv | cut -f 16 | sort | uniq -c
cat final_annovared_extended.tsv | cut -f 9 | sort | uniq -c
cat final_annovared_extended.tsv | cut -f 14 | sort | uniq -c
cat final_annovared_extended.tsv | cut -f 19 | sort | uniq -c
```

checking if ref/alt are the same: - yes!
```bash
awk -F'\t' '$4 != $457 || $5 != $458' final_annovared_extended.tsv | cut -f 4,5,457,458 | head 
```

Print rows with non-missing MAFs: 
```bash
awk -F'\t' '{if($303 != "." && $307 != "." && $306 != "." && $27 != ".") print $0}' final_annovared_extended.tsv > rows_with_mafs.tsv
```

how many of these rows with MAFs are exonic: ~ about 14000/22000
```bash
cat rows_with_mafs.tsv | cut -f 6 | grep -w "exonic" | wc -l
cat rows_with_mafs.tsv | cut -f 11 | grep -w "exonic" | wc -l
cat rows_with_mafs.tsv | cut -f 16 | grep -w "exonic" | wc -l
```

Making a table with rsID, ref/alt from databases, ref/alt from kazakh, kazakh MAFs and other population MAFs; then renaming the MAFs:
```bash
cut -f 456,4,5,457,458,687,303,307,306,27 rows_with_mafs.tsv | awk -F'\t' '{print $7, $1, $2, $8, $9, $10, $6, $5, $3, $4}' OFS='\t' > mafs_only.tsv
sed -i -e 's/kaz_alt_frq/Kazakh_MAF/g' \
       -e 's/AF_nfe/European_MAF/g' \
       -e 's/AF_eas/EastAsian_MAF/g' \
       -e 's/SAS.sites.2015_08/SouthAsian_MAF/g' \
       -e 's/AF_afr/African_MAF/g' mafs_only.tsv
```

Then I made sure that ref/alt are the same in Kazakh and ref populations:
```bash
awk '$2 != $4 || $3 != $5' mafs_only.tsv
```

Fold change:
```python
#!/usr/bin/env python3

import csv

# Input and output file names
input_file = "mafs_only.tsv"
output_file = "fold_change.tsv"

# Function to calculate fold change or handle the 0/NA cases
def calculate_fold_change(kazakh_maf, referent_maf):
    # Check for "NA" values (represented as ".")
    if kazakh_maf == "NA" or referent_maf == "NA":
        return "NA"
    kazakh_maf = float(kazakh_maf)
    referent_maf = float(referent_maf)

    # Handle the 0 cases or calculate the fold change
    if kazakh_maf == 0 and referent_maf == 0:
        return "kz=db=0"
    elif kazakh_maf == 0:
        return "kz=0"
    elif referent_maf == 0:
        return "DB=0"
    else:
        return kazakh_maf / referent_maf if kazakh_maf > referent_maf else -(referent_maf / kazakh_maf)

# Open the input file and process it
with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')
    
    # Write the header row
    header = next(reader)
    writer.writerow(["ID", "FoldChange_Europeans", "FoldChange_EastAsians", "FoldChange_SouthAsians", "FoldChange_Africans"])
    
    # Process each row
    for row in reader:
        snp = row[0]
        kazakh_maf = row[5] if row[5] != "." else "NA"  # Kazakh MAF from column 6
        european_maf = row[6] if row[6] != "." else "NA"  # European MAF from column 7
        east_asian_maf = row[7] if row[7] != "." else "NA"  # East Asian MAF from column 8
        south_asian_maf = row[8] if row[8] != "." else "NA"  # South Asian MAF from column 9
        african_maf = row[9] if row[9] != "." else "NA"  # African MAF from column 10

        # Calculate fold changes
        fold_change_europeans = calculate_fold_change(kazakh_maf, european_maf)
        fold_change_east_asians = calculate_fold_change(kazakh_maf, east_asian_maf)
        fold_change_south_asians = calculate_fold_change(kazakh_maf, south_asian_maf)
        fold_change_africans = calculate_fold_change(kazakh_maf, african_maf)

        # Write the result for this SNP
        writer.writerow([snp, fold_change_europeans, fold_change_east_asians, fold_change_south_asians, fold_change_africans])

print(f"Fold change table saved to {output_file}")
```

Executing the code:
```bash
chmod +x calculate_fold_change.py 
./calculate_fold_change.py
```

Modifying the fold change file to include the MAFs and making separate files for SNPs different with different populations; 
removing extra return characters
removing all lines where there is 0 value:

```bash
paste fold_change.tsv <(cut -f 6,7,8,9,10 mafs_only.tsv) | awk -F'\t' '{print $1, $6, $2, $7, $3, $8, $4, $9, $5, $10}' OFS='\t' > fold_change1.tsv
sed -i '' 's/\r//' fold_change1.tsv
grep -v '=' fold_change1.tsv > fold_change2.tsv
awk -F'\t' '$2 >= 0.02232143 && $4 >=  0.001 && $6 >=  0.001 && $8 >=  0.01 && $10 >=  0.001' fold_change2.tsv > fold_change3.tsv
```

checking how many exonic MAFs are left:
```bash
awk -F'\t' '$6 == "exonic" {print $456}' rows_with_mafs.tsv > exonic_rsIDs.tsv
awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a' exonic_rsIDs.tsv fold_change1.tsv | wc -l
awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a' exonic_rsIDs.tsv fold_change2.tsv | wc -l
awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a' exonic_rsIDs.tsv fold_change3.tsv | wc -l
```

Getting specific values different from each of these populations:
```bash
awk -F'\t' '{if (($3 > 5 || $3 < -5) && ($5 > 5 || $5 < -5) && ($7 > 5 || $7 < -5) && ($9 > 5 || $9 < -5)) print $0}' fold_change3.tsv > mafs_unique.tsv
awk -F'\t' '{if ($3 > 5 || $3 < -5) print $0}' fold_change3.tsv > mafs_change_euro.tsv
awk -F'\t' '{if ($5 > 5 || $5 < -5) print $0}' fold_change3.tsv > mafs_change_eastAsia.tsv
awk -F'\t' '{if ($7 > 5 || $7 < -5) print $0}' fold_change3.tsv > mafs_change_southAsia.tsv
awk -F'\t' '{if ($9 > 5 || $9 < -5) print $0}' fold_change3.tsv > mafs_change_afro.tsv
```

Checking how many exonic SNPs are left in each file: 4 246 148 99 374
```bash
awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a' exonic_rsIDs.tsv mafs_unique.tsv | wc -l
awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a' exonic_rsIDs.tsv mafs_change_euro.tsv | wc -l
awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a' exonic_rsIDs.tsv mafs_change_eastAsia.tsv | wc -l
awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a' exonic_rsIDs.tsv mafs_change_southAsia.tsv | wc -l
awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a' exonic_rsIDs.tsv mafs_change_afro.tsv | wc -l
```

Finding genes for lactase, alcohol, and pharmacogenes:
cat final_annovared_extended.tsv | grep -e "exonic" -e "ExonicFunc.knownGene" | grep -e "LCT" -e "ExonicFunc.knownGene" | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > lactose.tsv
cat final_annovared_extended.tsv | grep -e "exonic" -e "ExonicFunc.knownGene" | grep -e "ALDH2" -e "ADH" -e "ExonicFunc.knownGene" | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > alcohol.tsv
cat final_annovared_extended.tsv | grep -e "exonic" -e "ExonicFunc.knownGene" | grep -w -e "IFNL3" -e "NUDT15" -e "SLCO1B1" -e "TPMT" -e "UGT1A1" -e "CFTR" -e "CYP2B6" -e "CYP2C19" -e "CYP2C9" -e "CYP2D6" -e "CYP3A5" -e "CYP4F2" -e "DPYD" -e "VKORC1" -e "ExonicFunc.knownGene" | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > pharm_idda.tsv

cat final_annovared_extended.tsv | grep -e "exonic" -e "ExonicFunc.knownGene" | grep -w -f druggable_genome.tsv  | awk '$14 == "ExonicFunc.knownGene" || $70 == "D" && $76 == "D"' | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > druggable_mafs.tsv
```

```bash
awk -F'\t' '$14 == "ExonicFunc.knownGene" || $70 == "D" && $76 == "D"' final_annovared_extended.tsv | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > deleterious.tsv

cat deleterious.tsv | cut -f 12 | awk 'NR > 1 {sum += $1 * 224} END {print sum / 224}' 
# There was an average of 81.2555 likely deleterious nonsynonymous single nucleotide variants (SNVs) (those which are classified as deleterious by both SIFT and PolyPhen databases) per individual

awk -F'\t' '$70 == "D" && $76 == "D" {for (i=463; i<=686; i++) if ($i == "1/1") count[i]++} END {for (i=463; i<=686; i++) print "Individual " i-462 ": " count[i]}' final_annovared_extended.tsv 

awk -F'\t' '$70 == "D" && $76 == "D" {for (i=463; i<=686; i++) if ($i == "1/1") count[i]++} END {sum=0; for (i=463; i<=686; i++) sum+=count[i]; print "Average: ", sum/(686-463+1)}' final_annovared_extended.tsv
# An average of 19 variants per individual were homozygotes, therefore representing the presence of a considerable amount of potentially deleterious variation in the Kazakh population. 

awk -F'\t' '$14 == "ExonicFunc.knownGene" || $333 == "drug_response"' final_annovared_extended.tsv | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > clinvar_drug_response.tsv
awk -F'\t' '$14 == "ExonicFunc.knownGene" || $333 == "pathogenic"|| $333 == "risk_factor"' final_annovared_extended.tsv | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > clinvar_pathogenic_riskfactor.tsv





Calculating Fst as a measure of differential MAF:

```bash
nano calculate_fst.py
```

```python
#!/usr/bin/env python3
import pandas as pd
import sys

# Check if the correct number of arguments is passed
if len(sys.argv) != 3:
    print("Usage: calculate_fst.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Load the MAF table
maf_df = pd.read_csv(input_file, sep='\t')

# Define a function to calculate Fst between two populations based on their MAFs
def calculate_fst(maf1, maf2):
    p_avg = (maf1 + maf2) / 2  # Average allele frequency between populations
    H_T = 2 * p_avg * (1 - p_avg)  # Total heterozygosity
    H_S1 = 2 * maf1 * (1 - maf1)  # Subpopulation heterozygosity for population 1
    H_S2 = 2 * maf2 * (1 - maf2)  # Subpopulation heterozygosity for population 2
    H_S = (H_S1 + H_S2) / 2  # Average subpopulation heterozygosity
    if H_T == 0:  # To avoid division by zero
        return 0
    else:
        Fst = (H_T - H_S) / H_T  # Fst formula
        return Fst

# Calculate Fst for each SNP between Kazakhs and the other populations
fst_columns = ['Fst_European', 'Fst_EastAsian', 'Fst_SouthAsian', 'Fst_African']
for pop in fst_columns:
    population_maf = pop.split('_')[1] + '_MAF'
    maf_df[pop] = maf_df.apply(lambda row: calculate_fst(row['Kazakh_MAF'], row[population_maf]), axis=1)

# Save the results to a new file
maf_df.to_csv(output_file, sep='\t', index=False)

print(f"Fst calculations completed. Results saved to {output_file}")
```

chmod +x calculate_fst.py 
./calculate_fst.py mafs_only.tsv mafs_fst.tsv

cat mafs_fst.tsv | grep -f idda_pharm_rsIDs.txt  > idda_pharm_fsts.txt
