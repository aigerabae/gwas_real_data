genomestudio:

SampleSheet needs to have this format: 
[Header]
INVESTIGATOR NAME;Aygerim
PROJECT NAME;gwas300
EXPERIMENT NAME;gwas300
DATE;4.10.2024 5:57:45
[Manifests]
A;GSA-24v2-0_A1_hg19
[Data]
Sample_ID;SentrixBarcode_A;SentrixPosition_A;Path;Aux
1;203960450010;R01C01;C:\biostar\gwas300\new_gwas300\idats300\203960450010;0
16;203960450010;R01C02;C:\biostar\gwas300\new_gwas300\idats300\203960450010;0
2;203960450010;R02C01;C:\biostar\gwas300\new_gwas300\idats300\203960450010;0
17;203960450010;R02C02;C:\biostar\gwas300\new_gwas300\idats300\203960450010;0
4;203960450010;R03C01;C:\biostar\gwas300\new_gwas300\idats300\203960450010;0
18;203960450010;R03C02;C:\biostar\gwas300\new_gwas300\idats300\203960450010;0

(keep in mind that it's not exactly what google sheets gives when you downlaod as csv. instead of commas you have semicolons)

To get proper sample Names you need to make the sample sheet and load intensities from it. for data directory set up the general directiory or the idat directory (i dont rememebr which one is correct)

You can use clustering file from the genotyping company

Hopefully I can get plink plugin to work with proper Sample IDs
