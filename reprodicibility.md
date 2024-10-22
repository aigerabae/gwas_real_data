reproducibility solving:

# RAW intensity fiels to VCF:
```bash
# I removed 72 samples from my SampleSheet in google sheets to just have 224 filtered ones. I also changed path to idats be on the current disk
# Downloaded reference file for hg19 from https://support.illumina.com/downloads/genome-fasta-files.html

array-analysis-cli genotype call \
    --bpm-manifest /home/user/biostar/gwas/redo_october/making_vcf/hg19_manifest/GSA-24v2-0_A1_hg19.bpm \
    --cluster-file /home/user/biostar/gwas/redo_october/making_vcf/clustering_file/GSA-24v2-0_A1_ClusterFile.egt \
    --idat-sample-sheet /home/user/biostar/gwas/redo_october/making_vcf/sample_sheet/SampleSheet_224_aap_cli.csv \
    --output-folder /home/user/biostar/gwas/redo_october/making_vcf/gtc

dragena genotype gtc-to-vcf \
    --bpm-manifest /home/user/biostar/gwas/redo_october/making_vcf/hg19_manifest/GSA-24v2-0_A1_hg19.bpm \
    --genome-fasta-file /home/user/biostar/gwas/redo_october/making_vcf/GRCh37_genome/GRCh37_genome.fa \
    --gtc-sample-sheet /home/user/biostar/gwas/redo_october/making_vcf/sample_sheet/SampleSheet_224_aap_cli_with_gtc.csv \
    --csv-manifest /home/user/biostar/gwas/redo_october/making_vcf/hg19_manifest/GSA-24v2-0_A1.csv \
    --output-folder /home/user/biostar/gwas/redo_october/making_vcf/output_dragena_vcf

bcftools merge ./*.vcf.gz -o merged_output.vcf
# I renamed sentrix names to regular sample names
```

# Quality control of kaz data:
```bash
plink --vcf merged_output.vcf --make-bed --out custom_kaz
awk -F'\t' '!seen[$1]++' GSA-24v2-0_A1_b150_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > GSA-dictionary.txt

plink --bfile custom_kaz --update-name GSA-dictionary.txt --make-bed --out custom_kaz1
awk '$2 !~ /^rs/' custom_kaz1.bim | sort -k2,2 > custom_non_rs_SNP.txt
plink --bfile custom_kaz1 --exclude custom_non_rs_SNP.txt --make-bed --out custom_kaz2

# Make FID same as IID
cat custom_kaz2.fam | grep -wf list.txt | cut -d " " -f 1
# I manually changed 1 of 122 and 1 of 123 into 122_1 and 123_1 because ther were 2 of each of them

plink --bfile custom_kaz2 --remove "reproducibility - to_remove.tsv" --make-bed --out custom_kaz3
plink --bfile custom_kaz3 --update-ids "reproducibility - changing_fid.tsv" --make-bed --out custom_kaz4

plink --bfile custom_kaz4 --impute-sex --make-bed --out custom_kaz5
plink --bfile custom_kaz5 --geno 0.02 --make-bed --out custom_kaz6
plink --bfile custom_kaz6 --mind 0.02 --make-bed --out custom_kaz7
plink --bfile custom_kaz7 --maf 0.001 --make-bed --out custom_kaz8

# Attention! At this stage I removed all non-nucleotide SNPs (indels). But weirdly they are already not here. Probably got filteredt out at some point.
plink --bfile custom_kaz8 --snps-only 'just-acgt' --make-bed --out custom_kaz9
```

# Working with reference data:
```bash
# Merging with HGDP:
cp ../redo_july/working_with_ref_data/hgdp_estonian/HGDP9.bed ./
cp ../redo_july/working_with_ref_data/hgdp_estonian/HGDP9.bim ./
cp ../redo_july/working_with_ref_data/hgdp_estonian/HGDP9.fam ./

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

# admixture
plink --bfile all14 --indep-pairwise 1000 150 0.4 --out pruned_data
plink --bfile all14 --extract pruned_data.prune.in --make-bed --out all15
awk '{print $1}' all15.fam | grep -Fwf - ethnic2.txt > ethnic3.txt
for K in 8; do admixture --cv all15.bed -j8 $K | tee log${K}.out; done
ln -s ~/tools/AncestryPainter_v5/AncestryPainter.pl ./
cat ethnic3.txt | awk '{print $2"\t" $1}'  > ethnic3.ind 
perl AncestryPainter.pl -i ethnic3.ind -q ./all15.8.Q -t Kazakh -o Kazakh -f png
```

# Getting docs and graphs for publication:
```bash
# autosomal, mitochnodrial, y-chr
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' custom_kaz9.bim > snp_1_22.txt
plink --bfile custom_kaz9 --extract snp_1_22.txt --make-bed --out custom_kaz9_autosomal

awk '{ if ($1 == 26) print $2 }' custom_kaz9.bim > snp_mitoch.txt
plink --bfile custom_kaz9 --extract snp_mitoch.txt --make-bed --out custom_kaz9_mitoch

awk '{ if ($1 == 24) print $2 }' custom_kaz9.bim > snp_y.txt
plink --bfile custom_kaz9 --extract snp_y.txt --make-bed --out custom_kaz9_y_chr

#mafs
plink2 --bfile custom_kaz9  --freq --out maf_custom_kaz9
plink2 --bfile custom_kaz9_autosomal  --freq --out maf_custom_kaz9_autosomal
plink2 --bfile custom_kaz9_mitoch  --freq --out maf_custom_kaz9_mitoch
plink2 --bfile custom_kaz9_y_chr --freq --out maf_custom_kaz9_y_chr
sed -i 's/ALT_FREQS/kaz_alt_frq/g' maf_custom_kaz9_autosomal.afreq  maf_custom_kaz9_mitoch.afreq maf_custom_kaz9_y_chr.afreq maf_custom_kaz9.afreq
sed -i 's/OBS_CT/kaz_allele_count/g' maf_custom_kaz9_autosomal.afreq  maf_custom_kaz9_mitoch.afreq maf_custom_kaz9_y_chr.afreq maf_custom_kaz9.afreq

#making vcfs
plink --bfile custom_kaz9_autosomal --recode vcf --out custom_kaz9_autosomal
plink --bfile custom_kaz9_mitoch --recode vcf --out custom_kaz9_mitoch
plink --bfile custom_kaz9_y_chr --recode vcf --out custom_kaz9_y_chr
plink --bfile custom_kaz9 --recode vcf --out custom_kaz9_all

# making extended vcfs
bcftools view -h custom_kaz9_autosomal.vcf > kaz_a1.vcf && bcftools view -H custom_kaz9_autosomal.vcf > kaz_a2.vcf && tail -n +2 maf_custom_kaz9_autosomal.afreq | cut -f 6,7 > added_info.txt && head -n 1 maf_custom_kaz9_autosomal.afreq | cut -f 6,7 > added_header.txt && (sed '$d' kaz_a1.vcf; paste <(tail -n 1 kaz_a1.vcf) added_header.txt) > kaz_a4.vcf && paste kaz_a2.vcf added_info.txt > kaz_a3.vcf && cat kaz_a4.vcf kaz_a3.vcf > kaz_a5.vcf && mv kaz_a5.vcf ./autosomal_ext_for_annovar.vcf

# ROH
plink --bfile custom_kaz9_autosomal --homozyg-density 60 --homozyg-gap 500 --homozyg-window-snp 100 --homozyg-window-het 0

# IBD between kazakhs only and between different ethnicities
plink --bfile custom_kaz9_autosomal --genome --out ibd_kaz
plink --bfile all14 --genome --out ibd_all
./average_ibd.py ibd_all.genome ethnic_final.tsv
cat ibd_all_average.txt | grep Kazakh | sort -gk 3,3 > average_ibd_kazakh.tsv
```

# IGNORE! Work in progress notes 
```bash
# all but 1 mutation is present in the new dataset; the PCA and Fst graphs I updated; the claculations in the phram part kept intact due to no annotated vcf

# Let's work with indels and CNVs:
cat custom_kaz.bim | cut -f 2 | grep "CNV" | wc -l
# 2288 CNVs

cat custom_kaz.bim | cut -f 5-6 | grep -e "I" -e "D" | wc -l
# 8628 indels

# Same result, no quality score in dragena but it creates indexed gz file so its a bit more cinvenient. Can't do CNV because don't have CN model. Might not really need it since in the ped/map file i got from genomestudio all CNV variants are labeled as 0 amount

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

# comparing korean kaz and my custom kaz:
diff <(cut -f2  kaz2.bim | sort) <(cut -f2  custom_kaz2.bim | sort)
echo "Exclusive to kaz2.bim: $(comm -23 <(cut -f2 kaz2.bim | sort) <(cut -f2 custom_kaz2.bim | sort) | wc -l), Exclusive to custom_kaz2.bim: $(comm -13 <(cut -f2 kaz2.bim | sort) <(cut -f2 custom_kaz2.bim | sort) | wc -l), Common: $(comm -12 <(cut -f2 kaz2.bim | sort) <(cut -f2 custom_kaz2.bim | sort) | wc -l)"
# Exclusive to kaz2.bim: 34161, Exclusive to custom_kaz2.bim: 17188, Common: 641539

echo "Exclusive to kaz.bim: $(comm -23 <(cut -f2 kaz.bim | sort) <(cut -f2 custom_kaz.bim | sort) | wc -l), Exclusive to custom_kaz.bim: $(comm -13 <(cut -f2 kaz.bim | sort) <(cut -f2 custom_kaz.bim | sort) | wc -l), Common: $(comm -12 <(cut -f2 kaz.bim | sort) <(cut -f2 custom_kaz.bim | sort) | wc -l)"
# Exclusive to kaz.bim: 117857, Exclusive to custom_kaz.bim: 17244, Common: 648364

# I removed all old kaz files into separate folder and named it korean
