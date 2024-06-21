This file contains information about making PCA with world data:

 HGDP:
 accessed at https://www.hagsc.org/hgdp/files.html
 accessed metadata at https://www.internationalgenome.org/data-portal/data-collection/hgdp
To do:
1) extract only IDs into lists
2) plot PCA for Russian, Yakut, Japanese

extracted yakut, russian, and japanese IDs into their respective files
awk -F"\t" '{if ($5 == "Yakut") print $1}' metadata.txt > list_yakut.txt
awk -F"\t" '{if ($5 == "Japanese") print $1}' metadata.txt > list_japanese.txt
awk -F"\t" '{if ($5 == "Russian") print $1}' metadata.txt > list_russian.txt

replaced unnecessary spaces and tabs to make sure formatting is good
awk -F'\t' '{gsub(/[[:space:]]+/,"\t"); print}' HGDP.txt > HGDP1.txt

manually added word SNP as the first field of the first column -it was missing
transposed rows and columns so that each sample is a row and each SNP - a column

split -l 90000 HGDP1.txt file_
parallel 'datamash transpose < {file_{1..8}} > {.}_transposed.txt' 
paste file_{1..8}_transposed.txt > HGDP2.txt


got those rows corresponding to sample IDs in lists:




 Simons:
 accessed at https://www.simonsfoundation.org/simons-genome-diversity-project/ via cancer genomics cloud seven bridges
 accessed metadata at https://www.nature.com/articles/nature18964#Sec10 Supplementary Table 1

 1000 genomes
 accessed at https://www.internationalgenome.org/data-portal/population

To do:
1) figure out if my VCF file is good
