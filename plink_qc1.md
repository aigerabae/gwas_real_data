This file has code for command-line processing of HB00001157.ped and phenotypes.tsv files to generate the final version .ped file

1) removed header line 
```bash
sed 1d phenotypes.tsv > phenotypes1.tsv
```

2) removed extra spaces and tabs from ped file
```bash
awk -F'\t' '{gsub(/[[:space:]]+/,"\t"); print}' HB00001157.ped > kaz.ped
```

3) sorting phenotypes file and kaz.ped so that both are ordered the same according to ID column (2)
```bash
sort -t$'\t' -k2n,2 kaz.ped > kaz1.ped
sort -t$'\t' -k2n,2 phenotypes1.tsv > phenotypes2.tsv
```

view if they have the same number of rows - they don't
```bash
cat kaz1.ped | wc -l
cat phenotypes2.tsv | wc -l
```

4) only included those tows into phenotypes5.tsv that have counterparts in kaz1.ped
```bash
awk 'FNR==NR{a[$2]++; next} ($2 in a)' kaz1.ped phenotypes2.tsv > phenotypes3.tsv
```

checks if 1 value in 1 table matches 1 value in another
```bash
awk -F'\t' 'FNR==NR{a[$2]++; next} {b[$2]++} END{for (i in a) if (a[i]!=1) print "Value", i, "in phenotypes3.tsv does not have exactly one match in kaz1.ped"; for (j in b) if (b[j]!=1) print "Value", j, "in kaz1.ped matches more than one value in phenotypes3.tsv"}' kaz1.ped phenotypes3.tsv
```

5) saving separately the rows that have more than 1 matching
```bash
awk '$2 == 1314010' kaz1.ped > diff_1314010_kaz1.ped
awk '$2 == 1408910' kaz1.ped > diff_1408910_kaz1.ped
```

viewing those patients phenotypic data
```bash
awk '$2 == 1314010' phenotypes3.tsv 
awk '$2 == 1408910' phenotypes3.tsv 
```

so we have 4 patients that are different but have the same indices; 
5) excluding all 4 of them
```bash
awk -F'\t' '$2 != "1314010" && $2 != "1408910"' kaz1.ped > kaz2.ped
awk -F'\t' '$2 != "1314010" && $2 != "1408910"' phenotypes3.tsv > phenotypes4.tsv
```

check for match again
```bash
awk -F'\t' 'FNR==NR{a[$2]++; next} {b[$2]++} END{for (i in a) if (a[i]!=1) print "Value", i, "in phenotypes4.tsv does not have exactly one match in kaz2.ped"; for (j in b) if (b[j]!=1) print "Value", j, "in kaz2.ped matches more than one value in phenotypes4.tsv"}' kaz2.ped phenotypes4.tsv
cat kaz2.ped | wc -l
cat phenotypes4.tsv | wc -l
```

6) removed non-kazakh patients
```bash
awk -F'\t' '$18 == "kazakh"' phenotypes4.tsv > phenotypes5.tsv
awk 'FNR==NR{a[$2]++; next} ($2 in a)' phenotypes5.tsv kaz2.ped > kaz3.ped
```

check for matches
```bash
awk -F'\t' 'FNR==NR{a[$2]++; next} {b[$2]++} END{for (i in a) if (a[i]!=1) print "Value", i, "in phenotypes5.tsv does not have exactly one match in kaz3.ped"; for (j in b) if (b[j]!=1) print "Value", j, "in kaz3.ped matches more than one value in phenotypes5.tsv"}' kaz3.ped phenotypes5.tsv
cat kaz3.ped | wc -l
cat phenotypes5.tsv | wc -l
```
perfect!

7) I create a table phenotypes9.tsv that has numeric representations of data in phenotype8.tsv. 
1) I used Sample ID from column 2
2) I  used sex data in column 13 where Sex (1=male; 2=female; other=unknown)
```bash
cat phenotypes5.tsv | cut -f 2,13 > phenotypes6.tsv
awk -F"\t" '{ if ($2 == "male") $2 = 1; else if ($2 == "female") $2 = 2; print $0; }' phenotypes6.tsv > phenotypes7.tsv
sed 's/ /\t/g' phenotypes7.tsv > phenotypes8.tsv
```

8) checking number of columns in kaz3.ped and merging the ped file and the genders
```bash
awk '{print NF; exit}' kaz3.ped
cut -f 1-4 kaz3.ped > kaz4_1.ped
cut -f 2 phenotypes8.tsv > phenotypes9.tsv
paste kaz4_1.ped phenotypes9.tsv > kaz4_2.ped
cut -f 6-1532448 kaz3.ped > kaz4_3.ped
paste kaz4_2.ped kaz4_3.ped > kaz5.ped
```

9) Copy to another folder and rename kaz5.ped HB00001157.map as kaz.ped and kaz.map

Further analysis is avilable at plink_qc2.md
