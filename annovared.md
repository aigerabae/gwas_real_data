This script described work I did on annovar file from Aset:

Adding extended (KAZ_MAF + allele count) + changing ref/alt to kaz_ref/alt and DB (database) ref/alt:
```bash
cat autosomal_ext_for_annovar.vcf | tail -n +31 | cut -f 234,235 > for_extended_vcf.tsv
paste kaz_gwas_for_annovar_224.FINAL.annovar.hg38_multianno.header.txt for_extended_vcf.tsv > final_annovared_extended.tsv
sed -i '' -e 's/REF/kaz_ref/g' -e 's/ALT/kaz_alt/g' -e 's/Ref/DB_ref/g' -e 's/Alt/DB_alt/g' final_annovared_extended.tsv
```

Print how many missing values are in columns 60 to 111:
```bash
awk -F"\t" '{for(i=6;i<=111;i++) if($i == ".") count[i]++} END{for(i=6;i<=111;i++) print "Column " i ": " count[i]}' final_annovared_extended.tsv 
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

Print rows with non-missing MAFs:
```bash
awk -F'\t' '{for(i=57;i<=101;i++) if($i != ".") {print $0; next}; for(i=22;i<=26;i++) if($i != ".") {print $0; next}}' final_annovared_extended.tsv > rows_with_mafs.tsv
```

how many of these rows with MAFs are exonic: ~ about 1600
```bash
cat rows_with_mafs.tsv | cut -f 6 | grep -w "exonic" | wc -l
cat rows_with_mafs.tsv | cut -f 11 | grep -w "exonic" | wc -l
cat rows_with_mafs.tsv | cut -f 16 | grep -w "exonic" | wc -l
```

table with rsID, ref/alt from databases, ref/alt from kazakh, kazakh MAFs and all other MAFs:
```bash
cut -f 4,5,118,119,117,348,22-26,57-101 rows_with_mafs.tsv > mafs_only.tsv
```

In Excel I manually moved the last 4 columns to be the first 4 columns.

Then I made sure that ref/alt are the same in Kazakh and ref populations:
```bash
awk '$2 != $5 || $3 != $6' mafs_only.tsv
```

Gene list for GO: (KnownGene)
awk '$13 !~ /dist/ && $13 != "." {print $13}' final_annovared_extended.tsv > genelist_knowngene.txt
awk '$19 !~ /dist/ && $19 != "." {print $19}' final_annovared_extended.tsv > genelist_ensemble.txt
awk '$7 !~ /dist/ && $7 != "." {print $7}' final_annovared_extended.tsv > genelist_reqseq.txt

Fold change:
awk 'NR==1 {print $0, "Kazakh_MAF", "European_MAF", "EastAsian_MAF", "SouthAsian_MAF", "African_MAF", "MiddleEast_MAF"} 
NR>1 {getline file < "mafs_only.tsv"; split(file,maf,"\t"); print $0, maf[4], maf[54], maf[51], maf[56], maf[47], maf[53]}' fold_change_table.tsv > fold_change_with_mafs.tsv

cut -f 1,4,54,51,56,47,53 mafs_only.tsv > mafs_gnomad.tsv
sed -i '' -e 's/gnomad41_genome_AF_afr/afr_maf/g' -e 's/gnomad41_genome_AF_eas/east_asian_maf/g' -e 's/gnomad41_genome_AF_mid/mid_east_maf/g' -e 's/gnomad41_genome_AF_nfe/euro_maf/g' -e 's/gnomad41_genome_AF_sas/south_asia_maf/g'  mafs_gnomad.tsv
