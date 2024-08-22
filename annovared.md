This script described work I did on annovar file from Aset:

Adding extended (KAZ_MAF + allele count):
```bash
cat autosomal_ext_for_annovar.FINAL.annovar.hg38_multianno.header.txt | cut -f 347,348 > for_extended_vcf.tsv
paste kaz_gwas_for_annovar_224.FINAL.annovar.hg38_multianno.header.txt for_extended_vcf.tsv > final_annovared_extended.tsv
```

Print how many missing values are in columns 60 to 111:
```bash
awk -F"\t" '{for(i=6;i<=111;i++) if($i == ".") count[i]++} END{for(i=6;i<=111;i++) print "Column " i ": " count[i]}' kaz_gwas_for_annovar_224.FINAL.annovar.hg38_multianno.header.txt 
```

Print mutation type and exonic function (3 databases):
```bash
cat kaz_gwas_for_annovar_224.FINAL.annovar.hg38_multianno.header.txt | cut -f 6 | sort | uniq -c
cat kaz_gwas_for_annovar_224.FINAL.annovar.hg38_multianno.header.txt | cut -f 11 | sort | uniq -c
cat kaz_gwas_for_annovar_224.FINAL.annovar.hg38_multianno.header.txt | cut -f 16 | sort | uniq -c
cat kaz_gwas_for_annovar_224.FINAL.annovar.hg38_multianno.header.txt | cut -f 9 | sort | uniq -c
cat kaz_gwas_for_annovar_224.FINAL.annovar.hg38_multianno.header.txt | cut -f 14 | sort | uniq -c
cat kaz_gwas_for_annovar_224.FINAL.annovar.hg38_multianno.header.txt | cut -f 19 | sort | uniq -c
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

table with rsID, kazakh MAFs and all other MAFs:
```bash
cat rows_with_mafs.tsv | awk '{printf "%s ", $117; for(i=22;i<=26;i++) printf "%s ", $i; for(i=57;i<=101;i++) printf "%s ", $i; print ""}' > mafs_only.tsv
```
