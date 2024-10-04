reproducibility solving:

awk -F'\t' '!seen[$1]++' GSA-24v2-0_A1_b150_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > GSA-dictionary.txt

plink --bfile kaz --update-name GSA-dictionary.txt --make-bed --out kaz1
awk '$2 !~ /^rs/' kaz1.bim | sort -k2,2 > kaz_non_rs_SNP.txt
plink --bfile kaz1 --exclude kaz_non_rs_SNP.txt --make-bed --out kaz2

plink --bfile custom_kaz --update-name GSA-dictionary.txt --make-bed --out custom_kaz1
awk '$2 !~ /^rs/' custom_kaz1.bim | sort -k2,2 > custom_non_rs_SNP.txt
plink --bfile custom_kaz1 --exclude custom_non_rs_SNP.txt --make-bed --out custom_kaz2

diff <(cut -f2  kaz2.bim | sort) <(cut -f2  custom_kaz2.bim | sort)
echo "Exclusive to kaz2.bim: $(comm -23 <(cut -f2 kaz2.bim | sort) <(cut -f2 custom_kaz2.bim | sort) | wc -l), Exclusive to custom_kaz2.bim: $(comm -13 <(cut -f2 kaz2.bim | sort) <(cut -f2 custom_kaz2.bim | sort) | wc -l), Common: $(comm -12 <(cut -f2 kaz2.bim | sort) <(cut -f2 custom_kaz2.bim | sort) | wc -l)"
Exclusive to kaz2.bim: 34161, Exclusive to custom_kaz2.bim: 17188, Common: 641539

echo "Exclusive to kaz.bim: $(comm -23 <(cut -f2 kaz.bim | sort) <(cut -f2 custom_kaz.bim | sort) | wc -l), Exclusive to custom_kaz.bim: $(comm -13 <(cut -f2 kaz.bim | sort) <(cut -f2 custom_kaz.bim | sort) | wc -l), Common: $(comm -12 <(cut -f2 kaz.bim | sort) <(cut -f2 custom_kaz.bim | sort) | wc -l)"
Exclusive to kaz.bim: 117857, Exclusive to custom_kaz.bim: 17244, Common: 648364
