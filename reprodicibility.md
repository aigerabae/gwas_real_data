reproducibility solving:

```bash
awk -F'\t' '!seen[$1]++' GSA-24v2-0_A1_b150_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > GSA-dictionary.txt

plink --bfile kaz --update-name GSA-dictionary.txt --make-bed --out kaz1
awk '$2 !~ /^rs/' kaz1.bim | sort -k2,2 > kaz_non_rs_SNP.txt
plink --bfile kaz1 --exclude kaz_non_rs_SNP.txt --make-bed --out kaz2

plink --bfile custom_kaz --update-name GSA-dictionary.txt --make-bed --out custom_kaz1
awk '$2 !~ /^rs/' custom_kaz1.bim | sort -k2,2 > custom_non_rs_SNP.txt
plink --bfile custom_kaz1 --exclude custom_non_rs_SNP.txt --make-bed --out custom_kaz2

diff <(cut -f2  kaz2.bim | sort) <(cut -f2  custom_kaz2.bim | sort)
echo "Exclusive to kaz2.bim: $(comm -23 <(cut -f2 kaz2.bim | sort) <(cut -f2 custom_kaz2.bim | sort) | wc -l), Exclusive to custom_kaz2.bim: $(comm -13 <(cut -f2 kaz2.bim | sort) <(cut -f2 custom_kaz2.bim | sort) | wc -l), Common: $(comm -12 <(cut -f2 kaz2.bim | sort) <(cut -f2 custom_kaz2.bim | sort) | wc -l)"
# Exclusive to kaz2.bim: 34161, Exclusive to custom_kaz2.bim: 17188, Common: 641539

echo "Exclusive to kaz.bim: $(comm -23 <(cut -f2 kaz.bim | sort) <(cut -f2 custom_kaz.bim | sort) | wc -l), Exclusive to custom_kaz.bim: $(comm -13 <(cut -f2 kaz.bim | sort) <(cut -f2 custom_kaz.bim | sort) | wc -l), Common: $(comm -12 <(cut -f2 kaz.bim | sort) <(cut -f2 custom_kaz.bim | sort) | wc -l)"
# Exclusive to kaz.bim: 117857, Exclusive to custom_kaz.bim: 17244, Common: 648364

# I removed all old kaz files into separate folder and named it korean

# Make FID same as IID
$ cat custom_kaz2.fam | grep -wf list.txt | cut -d " " -f 1 

# I manually changed 1 of 122 and 1 of 123 into 122_1 and 123_1 because ther were 2 of each of them

plink --bfile custom_kaz2 --remove "reproducibility - to_remove.tsv" --make-bed --out custom_kaz3
plink --bfile custom_kaz3 --update-ids "reproducibility - changing_fid.tsv" --make-bed --out custom_kaz4

plink --bfile custom_kaz4 --impute-sex --make-bed --out custom_kaz5
plink --bfile custom_kaz5 --geno 0.02 --make-bed --out custom_kaz6
plink --bfile custom_kaz6 --mind 0.02 --make-bed --out custom_kaz7
plink --bfile custom_kaz7 --maf 0.001 --make-bed --out custom_kaz8
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' custom_kaz8.bim > snp_1_22.txt
plink --bfile custom_kaz8 --extract snp_1_22.txt --make-bed --out custom_kaz9

# Attention! At this stage I removed all non-nucleotide SNPs (indels). But weirdly they are already not here. Probably got filteredt out at some point.
plink --bfile custom_kaz9 --snps-only 'just-acgt' --make-bed --out custom_kaz10

# Merging with HGDP:
$ cp ../redo_july/working_with_ref_data/hgdp_estonian/HGDP9.bed ./
$ cp ../redo_july/working_with_ref_data/hgdp_estonian/HGDP9.bim ./
$ cp ../redo_july/working_with_ref_data/hgdp_estonian/HGDP9.fam ./

plink --bfile custom_kaz10 --bmerge HGDP9.bed HGDP9.bim HGDP9.fam --make-bed --out merged1
plink --bfile HGDP9 --exclude merged1-merge.missnp --biallelic-only strict --make-bed --out HGDP10
plink --bfile custom_kaz10 --exclude merged1-merge.missnp --biallelic-only strict --make-bed --out custom_kaz11
plink --bfile custom_kaz11 --bmerge HGDP10.bed HGDP10.bim HGDP10.fam --make-bed --out merged2
plink --bfile merged2 --geno 0.02 --make-bed --out merged3
plink --bfile merged3 --mind 0.02 --make-bed --out merged4
echo -e "rs9267522\nrs11229\nrs75412302\nrs12660719" > duplicates.txt
plink --bfile merged4 --exclude duplicates.txt --make-bed --out merged6
# merged6 instead of merged5 because in the original code i removed some outliers which i did earlier

# Strangely only 54000 SNPs remaining

wget https://evolbio.ut.ee/turkic/turkic.fam
wget https://evolbio.ut.ee/turkic/turkic.bim
wget https://evolbio.ut.ee/turkic/turkic.bed
wget https://evolbio.ut.ee/caucasus/caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand.bim
wget https://evolbio.ut.ee/caucasus/caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand.fam
wget https://evolbio.ut.ee/caucasus/caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand.bed
wget https://evolbio.ut.ee/sakha/sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand.bim
wget https://evolbio.ut.ee/sakha/sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand.fam
wget https://evolbio.ut.ee/sakha/sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand.bed
wget https://evolbio.ut.ee/jew/jew_paper_data_dbSNP-b131_pos-b37_1KG_strand.bim
wget https://evolbio.ut.ee/jew/jew_paper_data_dbSNP-b131_pos-b37_1KG_strand.fam
wget https://evolbio.ut.ee/jew/jew_paper_data_dbSNP-b131_pos-b37_1KG_strand.bed

plink --bfile caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand --bmerge jew_paper_data_dbSNP-b131_pos-b37_1KG_strand --make-bed --out all1
plink --bfile all1 --bmerge sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand --make-bed --out all2
plink --bfile all2 --bmerge turkic --make-bed --out all3

cp ../redo_july/working_with_ref_data/hgdp_estonian/metadata.txt ./

cat metadata.txt | awk '{print $1 "\t" $1}' > ethnic.txt 
plink --bfile all3 --keep ethnic.txt --make-bed --out all4

comm -12 <(awk '{print $2}' merged6.bim | sort) <(awk '{print $2}' all4.bim | sort) > common_snps.txt
cat merged6.bim | awk '{print $2"\t" $1}' > dictionary_chr
cat merged6.bim | cut -f 2,4 > dictionary_pos
plink --bfile all4 --extract common_snps.txt --make-bed --out all5
plink2 --bfile all5 --update-chr dictionary_chr --update-map dictionary_pos --sort-vars --make-pgen --out all6
plink2 --pfile all6 --make-bed --out all7
plink --bfile merged6 --bmerge all7 --make-bed --out all8
# this line above is for making a missnp file

plink --bfile all7 --exclude all8-merge.missnp --make-bed --out all8
plink --bfile merged6 --exclude all8-merge.missnp --make-bed --out merged7

plink --bfile merged7 --bmerge all8 --make-bed --out all9
cp ../redo_july/working_with_ref_data/hgdp_estonian/metadata_kaz.txt ./
cat metadata.txt metadata_kaz.txt > ethnic1.txt

plink --bfile all9 --geno 0.02 --make-bed --out all10
plink --bfile all10 --mind 0.02 --make-bed --out all11
plink --bfile all11 --maf 0.001 --make-bed --out all12
plink --bfile all12 --genome --min 0.2 --out pihat_min0.2
plink --bfile all12 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
plink --bfile all12 --snps-only 'just-acgt' --make-bed --out all13
awk '{print $1}' all13.fam | grep -Fwf - ethnic1.txt > ethnic2.txt

# I downloaded metadata_hgdp_estonian_kz_final - original_order.tsv and saved as ethnic_final.tsv
# also downlaoded pca and fst plotting scripts from git
# i also removed some samples that are present in all13 but not in metadata somehow
cat ethnic_final.tsv | awk '{print $1"\t"$1}' > to_keep_from_metadata.tsv
plink --bfile all13 --keep to_keep_from_metadata.tsv --make-bed --out all14

plink2 --bfile all14 --pca 10 --out all_pca
python plot_eigenvec.py all_pca.eigenvec ethnic_final.tsv
cat ethnic_final.tsv | awk '{print $1 "\t" $1 "\t" $2 "\t" $3}' > ethnic2_fst.txt
plink2 --bfile all13 --fst CATPHENO --within ethnic2_fst.txt --double-id --out fst_output
./plot_fst_heatmap.py fst_output.fst.summary sorting_order.tsv


# Making a table of MAFs and checking manually if i have the pharmacogenes there in the same amounts:
plink2 --bfile custom_kaz10  --freq --out maf_custom_kaz10

# comparing how many SNPs I have in common at custom_kaz10 stage and kaz12_autosomal:
cp ../redo_july/working_with_ref_data/hgdp_estonian/kaz12_autosomal.bed ./
cp ../redo_july/working_with_ref_data/hgdp_estonian/kaz12_autosomal.bim ./
cp ../redo_july/working_with_ref_data/hgdp_estonian/kaz12_autosomal.fam ./

echo "Exclusive to kaz12_autosomal.bim: $(comm -23 <(cut -f2 kaz12_autosomal.bim | sort) <(cut -f2 custom_kaz10.bim | sort) | wc -l), Exclusive to custom_kaz10.bim: $(comm -13 <(cut -f2 kaz12_autosomal.bim | sort) <(cut -f2 custom_kaz10.bim | sort) | wc -l), Common: $(comm -12 <(cut -f2 kaz12_autosomal.bim | sort) <(cut -f2 custom_kaz10.bim | sort) | wc -l)"

# Exclusive to kaz12_autosomal.bim:    24014, Exclusive to custom_kaz10.bim:     9689, Common:   513941

grep maf_custom_kaz10.afreq -w -e rs2108622 -e rs3745274 -e rs3745274 -e rs4148323 -e rs2070959 -e rs4988235 -e rs1573496 -e rs671 -e rs4148323 -e rs2056900 -e rs2076740 -e rs189261858 -e rs12484684 
# for some reason some SNPs didn't show in this grep search some I searched a few that were missing manually; 1 was missing (atopic dermatitis)

# i made the same table for kaz12_autosomal to make MAF comparison easier
plink2 --bfile kaz12_autosomal  --freq --out maf_kaz12_autosomal
grep maf_kaz12_autosomal.afreq -w -e rs2108622 -e rs3745274 -e rs3745274 -e rs4148323 -e rs2070959 -e rs4988235 -e rs1573496 -e rs671 -e rs4148323 -e rs2056900 -e rs2076740 -e rs189261858 -e rs12484684 
# same here - some had to be searched manually

# Thank god... all but one SNP from ones I described are present in the re-done dataset with almost identical MAFs (I suppose the tiny difference in 1 of them comes from different calling algrithms but I'd say its negligible). The one that isn't present is for atopic dermatitis but I'm sure I can come up with something if I take a look at insertions and deletions instead. Yay! I have re-done PCA and FST (although now I have considerably less SNPs) but it still looks oretty much the same. Now I need to redo the figure with the mutations and potentially recalculate the numbers with different mutations ("average of 79 non-synymous mutations per person" or something alogn the lines). I have uploaded all work I did in that redo_october folder in Google DRive to continue on Monday on the workstation.
```
