Work I do on Affymetrix data; unrelated to GWAS 300 Kazakh analysis

headers description = rows 1-5 
column names = row 6

columns 2-805 = patients genotypes
column 807 - chromosome id
808 - start
809 - stop
810 - strand
- how many have 1 or another strand? = all strands are positive (+)

811 - rsID
819 - allele A
820 - allele B
821 - Ref Allele
822 - Alt Allele
836 - extended rsID
- whats the differemt with dbSNP?

864 - nAA
865 - nAB
866 - nBB

867 - nA
868 - nB
869 - n_CN0
870 - n_NC
- how many nA are present? = 2029 non-0 A counts
- how many nB are present? = 2053
- how many n_CN0 are present? = 1120
- how many n_NC are present? = weird result, only a small amount of 0s
- In general, how many non AA, AB,BB genotypes are present?

Ped/map files ~ what I am going to extract:
Map file: 807 (chromosome), 811 (rsID), 808 (start), 809 (stop)
- extract only those rows where 808=809
Ped file: 0s (FID), row 6 columns 2-805 (IID, transpose), 0s (PID), 0s (MID), Sex (from phenotype data), Phenotype (from phenotype data), Calls row 7 to 790.125 columns 2-805 (transpose)
What columns do I need to save?
2-805,807,811,808,809,

1) transpose
Remove header and select certain columns; Split,transpose (fastest way - with csvtk!),put back together:
```bash
cat ah.txt | tail -n +6 | cut -f 2-805,807,808,809,836 > ah1.tsv
split -l 50000 -d ah1.tsv ah2_chunks
seq -w 00 15 | parallel -j 8 'csvtk -t transpose ah2_chunks{} > ah3_chunks{}'
paste ah3_chunks00 ah3_chunks01 ah3_chunks02 ah3_chunks03 ah3_chunks04 ah3_chunks05 ah3_chunks06 ah3_chunks07 ah3_chunks08 ah3_chunks09 ah3_chunks10 ah3_chunks11 ah3_chunks12 ah3_chunks13 ah3_chunks14 ah3_chunks15 > ah4.tsv
```

Result: 
rows 1 to 804 - patients
805 - chrmosome
806 - start
807 - end
808 - rsID

cat ah4.tsv | tail -n 4 > ah5_1.map
cat ah4.tsv | tail -n 1 > ah5_1.ped
cat ah4.tsv | wc -l
cat ah4.tsv | head -n 804 > ah5_2.ped
cat ah5_1.ped ah5_2.ped > ah5_3.ped

Map file:
805,808,0s,806,807(?)
Ped file:
0s, 1(?), 0s, 0s, sex, phenotype, 2-804
