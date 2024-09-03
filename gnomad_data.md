Downloading gnomad data (HGDP + 1000 Genomes):

gsutil cp gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/f2_fst/hgdp_tgp.bim ./
gsutil cp gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/f2_fst/hgdp_tgp.bed ./
gsutil cp gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/f2_fst/hgdp_tgp.fam ./

gnomad data:
Getting rsIDs:
plink --bfile hgdp_tgp --keep <(echo -e "Surui HGDP00843") --recode vcf --out small
cp ~/tools/annovar/convert2annovar.pl ./
cp ~/tools/annovar/annotate_variation.pl ./
perl convert2annovar.pl -format vcf4 small.vcf > input.avinput
perl annotate_variation.pl --downdb avsnp151 -buildver hg38 -webfrom annovar humandb/

