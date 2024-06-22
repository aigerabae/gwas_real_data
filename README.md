This is a methods section for GWAS study using 300 Kazakh genomes.

**I) Information about sequencing**
1) 296 Samples, 766221 SNPs were genotyped 
2) Platform: GSA MG v2
3) Scan protocol: Standard Illumina procedures using Illumina iScan scanner
4) Reference genome: GRCh38 Genome Reference Consortium Human Build 38

**II) Raw data processing**	

Normalized signal intensity and genotype were computed using Illumina’s GenomeStudio v.2 software.
1) Clustering
  - Make genotype calls across all samples using a standard Infinium Bead Chip cluster file.
  - The standard cluster file (*egt file) supplied by Illumina for each Infinium BeadChip type is generated using a diverse set of more than 200 HapMap1DNA samples in an Illumina laboratory. Some SNP probes (include custom probes) were clustered manually for reduce the number of spurious region calls, and increase the accuracy of the results.
custom cluster tech note

2) Callrate check
  - Quality checking by plotting p10 GC and sample call rate.
<img width="600" alt="Screenshot 2024-06-19 at 15 13 25" src="https://github.com/aigerabae/gwas_real_data/assets/155903885/956ffb0b-db4f-4041-94c8-ca7f010bd87e">
 
  - If a sample passes the intact DNA sample QC criteria but when the callrate is significantly lower than other samples, then the sample can be re-experiment is some cases. Call rates were consistently high in the experiment; no samples were excluded at this stage

* Call rate, p10 GC, GenCall Score are available in phenotypes.tsv 
* _Call Rate_: Percentage of SNPs (expressed as a decimal) whose GenCall score is greater than the specified threshold.
* 	_p10 GC_: 10th percentile GenCall score over all SNPs for this sample.
* 	 GenCall Score_: quality metric that indicates the reliability of each genotype call.

3) Genotype matrix export
Make a text file that contains the genotype of entire samples and probes.

4) Make input files for third party tools: *.ped & *.map files to execute PLINK.
- HB00001157.ped, HB00001157.map were made available by sequencing provider
- phenotypes.tsv file that contained patient information was provided separately


**III) Data processing in PLINK** (command-line code in plink_qc)
1) prepared .map file was used for further analysis
2) prepared .ped file needed further processing to include phenotypes (analysis shown in plink_qc1.md)
3) prepared .ped and .map files were processed for quality control (analysis shown in plink_qc2.md)
4) VCF file was generated (analysis shown in plink_qc3.md)
5) PCA was performed using other populations from Human Genome Diversity Project (HGDP) (analysis shown in world_pca.md)
