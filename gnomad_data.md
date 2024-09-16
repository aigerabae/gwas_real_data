Downloading gnomad data (HGDP + 1000 Genomes):

```bash
gsutil cp gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/f2_fst/hgdp_tgp.bim ./
gsutil cp gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/f2_fst/hgdp_tgp.bed ./
gsutil cp gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/f2_fst/hgdp_tgp.fam ./
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README_populations.md
```

getting metadata, formatting it, and finding samples from eurasian populations (from google sheets manually created list):
```bash
cat gnomad_meta_updated.tsv | cut -f 1,176 > pops.tsv
sed -i 's/[[:space:]]*$//' eurasian_to_keep.tsv
sed -i 's/[[:space:]]*$//' pops.tsv
grep -Ff eurasian_to_keep.tsv pops.tsv > samples_to_keep.tsv
```

Updating SNP names in gnomad data using kazakh dataset as a reference:
Removing duplicate SNPs in kazakh and gnomad datasets:
```bash
awk '!seen[$1,$4]++' kaz3.bim > kaz4.bim
awk '!seen[$1,$4]++' gnomad.bim > gnomad1.bim
```

Creating a dictionary of common SNPs with gnoamd names and their rsID namesd based on position:
```bash
comm -12 <(awk '{print $1, $4}' kaz3.bim | sort) <(awk '{print $1, $4}' gnomad.bim | sort) | awk 'NR==FNR{pos[$1 FS $2]; next} ($1 FS $4) in pos {print $2}' - kaz4.bim > common_rsIDs.txt
comm -12 <(awk '{print $1, $4}' kaz3.bim | sort) <(awk '{print $1, $4}' gnomad.bim | sort) | awk 'NR==FNR{pos[$1 FS $2]; next} ($1 FS $4) in pos {print $2}' - gnomad1.bim > common_gnomad_names.txt
paste common_gnomad_names.txt common_rsIDs.txt > dict_names.txt
```

```bash
plink --bfile gnomad --extract common_gnomad_names.txt --make-bed --out gnomad1
plink2 --bfile gnomad1 --update-name dict_names.txt 2 1 --make-bed --out gnomad2
```

Obtaining a list of samples with IID + FID to keep in gnomad dataset:
```bash
cat samples_to_keep_metadata.tsv | cut -f 1 > samples_to_keep.tsv
awk 'NR==FNR {fam[$2] = $1; next} {if ($1 in fam) print fam[$1], $1}' gnomad2.fam samples_to_keep.tsv > samples_to_keep2.tsv
plink --bfile gnomad2 --keep samples_to_keep2.tsv --make-bed --out gnomad3
```

Running admixture on gnomad data without kazakhs: ~ still not very nice; even worse than just estonian + hgdp
```bash
for K in 8; do admixture --cv gnomad3.bed -j8 $K | tee log${K}.out; done
ln -s ~/tools/AncestryPainter_v5/AncestryPainter.pl ./
cat samples_to_keep_metadata.tsv | awk '{print $3"\t" $1}'  > samples.ind
grep -Ff <(cut -d " " -f2 samples_to_keep2.tsv) samples.ind > samples1.ind
perl AncestryPainter.pl -i samples1.ind -q ./gnomad3.8.Q -t Uygur -o Uygur_adm  -f png
```
