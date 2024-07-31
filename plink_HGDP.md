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
nano remove_relatives.sh
```

```bash
#!/bin/bash

# Check if the input file exists
if [ ! -f "related_pairs.txt" ]; then
  echo "File not found!"
  exit 1
fi

# Create a working copy of the related pairs file
cp related_pairs.txt related_pairs_working.txt

# Extract all IDs from the file and count their occurrences
awk '{print $2; print $4}' related_pairs_working.txt | sort | uniq -c | sort -nr > id_count.txt

# Initialize an empty list for removal
> to_remove.txt

# Keep removing the individual with the highest count until no pairs remain
while [ -s related_pairs_working.txt ]; do
  # Get the individual with the highest count
  to_remove=$(head -n 1 id_count.txt | awk '{print $2}')
  
  # Add this individual to the removal list with FID and IID columns
  echo "$to_remove $to_remove" >> to_remove.txt
  
  # Remove all pairs involving this individual from related_pairs_working.txt
  awk -v remove="$to_remove" '$2 != remove && $4 != remove' related_pairs_working.txt > temp.txt
  mv temp.txt related_pairs_working.txt
  
  # Recount the occurrences of each individual
  awk '{print $2; print $4}' related_pairs_working.txt | sort | uniq -c | sort -nr > id_count.txt
done

echo "Individuals to be removed (FID IID):"
cat to_remove.txt
```

```bash
plink --bfile HGDP6 --remove to_remove.txt --make-bed --out HGDP7
plink --bfile HGDP7 --snps-only 'just-acgt' --make-bed --out HGDP8
```

4) runs of homozygosity (ROH)
plink --bfile kaz12_autosomal --homozyg-density 60 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-window-het 0
plink --bfile HGDP7 --homozyg-density 60 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-window-het 0

tasks for thursday:
- understand whats Fst
- check how to visualize ROH
- check how to combine HGDP and kaz data
- redo PCA: include only chosen populations into the final dataset, see how many of them will be present, rewrite the python script for making the PCA graph

use all populations from HDGP and SGDP for ROH and Fst and only selected eurasian populations for PCA
