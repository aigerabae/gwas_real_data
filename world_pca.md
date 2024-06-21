This file contains information about making PCA with world data:

 HGDP:
 accessed at https://www.hagsc.org/hgdp/files.html
 accessed metadata at https://www.internationalgenome.org/data-portal/data-collection/hgdp
To do:
1) extract only IDs into lists
2) plot PCA for Russian, Yakut, Japanese

extracted yakut, russian, and japanese IDs into their respective files
awk -F"\t" '{if ($5 == "Yakut") print $1, $2}' metadata.txt > list_yakut.txt
awk -F"\t" '{if ($5 == "Japanese") print $1, $2}' metadata.txt > list_japanese.txt
awk -F"\t" '{if ($5 == "Russian") print $1, $2}' metadata.txt > list_russian.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print }' list_yakut.txt > list_yakut1.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print }' list_russian.txt > list_russian1.txt
awk '{ if ($2 == "female") $2 = 2; else if ($2 == "male") $2 = 1; print }' list_japanese.txt > list_japanese1.txt

replaced unnecessary spaces and tabs to make sure formatting is good
awk -F'\t' '{gsub(/[[:space:]]+/,"\t"); print}' HGDP.txt > HGDP1.txt

manually added word SNP as the first field of the first column -it was missing
transposed rows and columns so that each sample is a row and each SNP - a column

split -l 90000 HGDP1.txt file_
I manually renamed them into gtReport1 to gtReport8, then using code bwlow removed spaces and replaced them with tabs, then transposed it 
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

got those rows corresponding to sample IDs in lists:
awk 'NR==FNR {a[$1]; next} $1 in a' list_japanese1.txt HGDP2.txt > japanese_SNP.txt
awk 'NR==FNR {a[$1]; next} $1 in a' list_russian1.txt HGDP2.txt > russian_SNP.txt
awk 'NR==FNR {a[$1]; next} $1 in a' list_yakut1.txt HGDP2.txt > yakut_SNP.txt

making proper map/ped files from them: 
merged 3 together
cat japanese_SNP.txt russian_SNP.txt yakut_SNP.txt > all.ped
cat list_japanese1.txt list_russian1.txt list_yakut1.txt > all.ped





 Simons:
 accessed at https://www.simonsfoundation.org/simons-genome-diversity-project/ via cancer genomics cloud seven bridges
 accessed metadata at https://www.nature.com/articles/nature18964#Sec10 Supplementary Table 1

 1000 genomes
 accessed at https://www.internationalgenome.org/data-portal/population

To do:
1) figure out if my VCF file is good
