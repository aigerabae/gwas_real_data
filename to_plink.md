Affyemtrix data:
Remove the first 5 rows of ah.txt and save as ah1.txt
remove rows where column 836 is empty in ah1.txt and save as ah2.txt
remove rows where column 835 doesn't equal 2 and save as ah3.txt
for every row in columns 2-805 inclusive: replace "A" with content of coplumn 819 from that same row; replace "B" with content of coplumn 820 from that same row; replace "NoCall" with "0" and save as ah4.txt
in row 1 of ah4.txt columns 2 to 805 inclusive remove in each cell everything after "_(" (inclusive); also remove "AH" from cells; save as ah5.txt
save in ah5.txt columns 807, 836, 0s, 808 (use awk in that exact order) as ah.map
save columns 2 to 805 inclusive as genotypes.txt
transpose genotypes.txt and save as genotypes2.txt
save head_ah.tsv second sheet (id_sex_pheno) from personal google sheets account
in id_sex_pheno.tsv remove rows where columns 2 and 3 have empty values and save that as id_sex_pheno1.tsv
save column 1 of genotypes2.txt as samples.tsv
use join to get a table with samples.tsv joined by columns 2 and 3 from id_sex_pheno1.tsv; save as id_sex_pheno2.tsv
save ah1.ped with column 1 = all 0s; column 2 = column 1 of genotypes2.txt; column 3 all 0s; column 4 = all 0s; column 5 = column 2 of id_sex_pheno2.tsv; column 6 = column 3 of id_sex_pheno2.tsv; columsn 7 to end = columns 2 to end of genotypes2.txt

```bash
#!/bin/bash

# Step 1: Remove the first 5 rows from ah.txt and save as ah1.txt
tail -n +6 ah.txt > ah1.txt

# Step 2: Remove rows where column 836 is empty from ah1.txt and save as ah2.txt
awk '$836 != ""' ah1.txt > ah2.txt

# Step 3: Remove rows where column 835 doesn't equal 2 from ah2.txt and save as ah3.txt
awk '$835 == 2' ah2.txt > ah3.txt

# Step 4: For columns 2-805, replace "A" with col 819, "B" with col 820, and "NoCall" with "0", save as ah4.txt
awk '{
  for(i=2; i<=805; i++) {
    if($i == "A") $i = $819;
    else if($i == "B") $i = $820;
    else if($i == "NoCall") $i = "0";
  }
  print
}' ah3.txt > ah4.txt

# Step 5: Modify the header row in ah4.txt, saving the result as ah5.txt
awk 'NR == 1 {
  for(i=2; i<=805; i++) {
    sub(/_\(.*$/, "", $i);
    sub(/AH/, "", $i);
  }
}
{ print }
' ah4.txt > ah5.txt

# Step 6: Save columns 807, 836, 0s, 808 in that order to ah.map
awk '{ print $807, $836, "0s", $808 }' ah5.txt > ah.map

# Step 7: Save columns 2-805 to genotypes.txt
awk '{ for(i=2; i<=805; i++) printf "%s ", $i; print "" }' ah5.txt > genotypes.txt

# Step 8: Transpose genotypes.txt and save as genotypes2.txt
awk '
{
  for (i=1; i<=NF; i++)  {
    a[NR,i] = $i
  }
}
NF>p { p = NF }
END {
  for(j=1; j<=p; j++) {
    str=a[1,j]
    for(i=2; i<=NR; i++) {
      str=str" "a[i,j];
    }
    print str
  }
}' genotypes.txt > genotypes2.txt

# Step 9: In id_sex_pheno.tsv, remove rows where columns 2 and 3 are empty and save as id_sex_pheno1.tsv
awk '$2 != "" && $3 != ""' id_sex_pheno.tsv > id_sex_pheno1.tsv

# Step 10: Save column 1 of genotypes2.txt as samples.tsv
awk '{ print $1 }' genotypes2.txt > samples.tsv

# Step 11: Use join to get a table with samples.tsv joined by columns 2 and 3 from id_sex_pheno1.tsv and save as id_sex_pheno2.tsv
join -1 1 -2 1 samples.tsv <(awk '{print $1, $2, $3}' id_sex_pheno1.tsv) > id_sex_pheno2.tsv

# Step 12: Save ah1.ped
awk '{
  if (NR > 1) print "0", $1, "0", "0", $2, $3, $0;
}' genotypes2.txt > ah1.ped
```
