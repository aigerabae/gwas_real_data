This file contains information about making PCA with world data:

 HGDP:
 accessed at https://www.hagsc.org/hgdp/files.html
 accessed metadata at https://www.internationalgenome.org/data-portal/data-collection/hgdp
To do:
1) extract only IDs into lists
2) plot PCA for Russian, Yakut, Japanese, Uygur, Sardinian, Adygei, Hazara, Mongolian, Northern Han + Kazakh

replaced unnecessary spaces and tabs to make sure formatting is good in original pedigree file
```bash
awk -F'\t' '{gsub(/[[:space:]]+/,"\t"); print}' HGDP.txt > HGDP1.txt
```

manually added word SNP as the first field of the first column -it was missing
transposed rows and columns so that each sample is a row and each SNP - a column

split -l 90000 HGDP1.txt file_
I manually renamed them into gtReport1 to gtReport8, then using code bwlow removed spaces and replaced them with tabs, then transposed it 
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
paste gtReport{1..8}_transposed.tsv > concatenated_transposed.tsv
mv concatenated_transposed.tsv ./HGDP2.txt
```

extracted yakut, russian, and japanese IDs into their respective files
```bash
awk -F"\t" '{if ($5 == "Yakut") print $1, $2}' metadata.txt > list_yakut.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print $1 "\t" $2 }' list_yakut.txt > list_yakut1.txt
awk -F"\t" '{if ($5 == "Japanese") print $1, $2}' metadata.txt > list_japanese.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print $1 "\t" $2 }' list_japanese.txt > list_japanese1.txt
awk -F"\t" '{if ($5 == "Russian") print $1, $2}' metadata.txt > list_russian.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print $1 "\t" $2 }' list_russian.txt > list_russian1.txt
awk -F"\t" '{if ($5 == "Uygur") print $1, $2}' metadata.txt > list_uyghur.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print $1 "\t" $2 }' list_uygur.txt > list_uygur1.txt
awk -F"\t" '{if ($5 == "Sardinian") print $1, $2}' metadata.txt > list_sardinian.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print $1 "\t" $2 }' list_sardinian.txt > list_sardinian1.txt
awk -F"\t" '{if ($5 == "Adygei") print $1, $2}' metadata.txt > list_adygei.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print $1 "\t" $2 }' list_adygei.txt > list_adygei1.txt
awk -F"\t" '{if ($5 == "Hazara") print $1, $2}' metadata.txt > list_hazara.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print $1 "\t" $2 }' list_hazara.txt > list_hazara1.txt
awk -F"\t" '{if ($5 == "Mongolian") print $1, $2}' metadata.txt > list_mongolian.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print $1 "\t" $2 }' list_mongolian.txt > list_mongolian1.txt
awk -F"\t" '{if ($5 == "Northern Han") print $1, $2}' metadata.txt > list_northernhan.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print $1 "\t" $2 }' list_northernhan.txt > list_northernhan1.txt
```

got those rows corresponding to sample IDs in lists:
```bash
awk 'NR==FNR {a[$1]; next} $1 in a' list_japanese1.txt HGDP2.txt > japanese_SNP.txt
awk 'NR==FNR {a[$1]; next} $1 in a' list_russian1.txt HGDP2.txt > russian_SNP.txt
awk 'NR==FNR {a[$1]; next} $1 in a' list_yakut1.txt HGDP2.txt > yakut_SNP.txt
awk 'NR==FNR {a[$1]; next} $1 in a' list_uygur1.txt HGDP2.txt > uygur_SNP.txt
awk 'NR==FNR {a[$1]; next} $1 in a' list_sardinian1.txt HGDP2.txt > sardinian_SNP.txt
awk 'NR==FNR {a[$1]; next} $1 in a' list_adygei1.txt HGDP2.txt > adygei_SNP.txt
awk 'NR==FNR {a[$1]; next} $1 in a' list_hazara1.txt HGDP2.txt > hazara_SNP.txt
awk 'NR==FNR {a[$1]; next} $1 in a' list_mongolian1.txt HGDP2.txt > mongolian_SNP.txt
awk 'NR==FNR {a[$1]; next} $1 in a' list_northernhan1.txt HGDP2.txt > northernhan_SNP.txt
```

making proper map/ped files from them: 
merged 3 together
```bash
cat list_japanese1.txt list_russian1.txt list_yakut1.txt list_uygur1.txt list_sardinian1.txt list_adygei1.txt list_hazara1.txt list_mongolian1.txt list_northernhan1.txt > all_1_1.ped
cat japanese_SNP.txt russian_SNP.txt yakut_SNP.txt uygur_SNP.txt sardinian_SNP.txt adygei_SNP.txt hazara_SNP.txt mongolian_SNP.txt northernhan_SNP.txt > all_1_2.ped
```

splitting genotypes into 2 columns
```bash
awk 'BEGIN {FS=OFS="\t"} {for(i=2; i<=NF; i++) gsub(/./,"\t&",$i)} 1' all_1_2.ped > all_1_2_1.ped
```

removing extra tabs
```bash
awk -F'\t' '{gsub(/[[:space:]]+/,"\t"); print}' all_1_2_1.ped > all_1_2_2.ped
```

adding family id, maternal ID, paternal ID, and phenotype
```bash
cut -f1 all_1_1.ped > temp_column.txt
paste temp_column.txt all_1_1.ped > all_1_3.ped
awk '{print $1, $2, 0, 0}' OFS="\t" all_1_3.ped > all_1_4.ped
paste all_1_4.ped <(cut -f2 all_1_1.ped) <(cut -f3 all_1_4.ped) <(cut -f2- all_1_2_2.ped) > all_2.ped
awk '{print $2, $1, 0, $3}' OFS="\t" HGDP_Map.txt > all.map
```

checking if all rows have the same number of fields
```bash
awk '{print NF}' all_2.ped
```

adding missing fields to lines that have less fields
```bash
awk 'NR==1 {num_fields=NF; header=$0; next} 
{for (i=NF+1; i<=num_fields; i++) $i=0; print} 
END {print header}' all_2.ped > all_2_1.ped
```

checking if all rows have the same number of fields
```bash
awk '{print NF}' all_2_1.ped
```

to plink binary
```bash
cp all_2_1.ped ./3all.ped
cp all.map ./3all.map
plink --file 3all --missing-code -9,0,NA,na --make-bed --out all2
```

merging my dataset with world dataset
```bash
??? awk '{print $2}' all2.bim > all2_rsids.txt
awk 'NR==FNR{a[$1];next}($2 in a){print $2,$4}' all2_rsids.txt kaz9.bim > mapped_rsids.txt
plink --bfile kaz9 --extract  mapped_rsids.txt --make-bed --out kaz10
plink --bfile all2 --extract  mapped_rsids.txt --make-bed --out all3
plink --bfile all3 --update-map mapped_rsids.txt --make-bed --out all4
```

```bash
echo merge_list.txt > kaz10.bed       kaz10.bim       kaz10.fam
```

```bash
plink --bfile all4 --merge-list merge_list.txt --make-bed --out merged_dataset 
```

PCA
```bash
plink --bfile merged_dataset --geno 0.02 --make-bed --out merged_dataset1
plink --bfile merged_dataset1 --mind 0.02 --make-bed --out merged_dataset2
plink2 --bfile merged_dataset2 --pca 10 --out all_pca 

cat metadata.txt | cut -f 1,5 > ethicities.txt
cat kaz10.fam | cut -f 1 -d ' ' >> temp_ethicities.txt 
cat temp_ethicities.txt | awk '{$2 = "\tkazakh"; print }' >> ethicities.txt

python plot_eigenvec.py all_pca.eigenvec ethicities.txt
```

where plot_eigenvec.py has script:
```python
import sys
import pandas as pd
import plotly.express as px

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python plot_eigenvec.py <eigenvec_file> <ethnicities_file>")
    sys.exit(1)

# Read the eigenvec file, skipping the first line
eigenvec_file = sys.argv[1]
with open(eigenvec_file, 'r') as f:
    lines = f.readlines()[1:]

# Read the ethnicities file
ethnicities_file = sys.argv[2]
with open(ethnicities_file, 'r') as f:
    ethnicities = {line.split()[0]: line.split()[1] for line in f.readlines()}

# Extract data from the file
data = {
    'PC1': [float(line.split()[2]) for line in lines],
    'PC2': [float(line.split()[3]) for line in lines],
    'sample_id': [line.split()[0] for line in lines],
    'ethnicity': [ethnicities[line.split()[0]] for line in lines]
}

df = pd.DataFrame(data)

# Create interactive scatter plot
fig = px.scatter(df, x='PC1', y='PC2', color='ethnicity', hover_data=['sample_id', 'ethnicity'])

# Save the plot as an HTML file
fig.write_html('interactive_plot.html')

# Show the plot in a browser
fig.show()
```

To do: include other populations for central asian PCA - Pathan, Uygur

More data to use for PCA:

Simons:
can be accessed at https://www.simonsfoundation.org/simons-genome-diversity-project/ via cancer genomics cloud seven bridges
metadata can be accessed at https://www.nature.com/articles/nature18964#Sec10 Supplementary Table 1

1000 genomes
can be accessed at https://www.internationalgenome.org/data-portal/population
