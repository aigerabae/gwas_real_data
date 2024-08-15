Directory: ~/biostar/gwas/redo_july/working_with_ref_data/ancient/

```bash
sed 's/-9/1/g' kaz12_autosomal.fam > kaz12_autosomal_1.fam
```

```bash
nano convertf_param.par
```
```bash
genotypename:   kaz12_autosomal.bed
snpname:        kaz12_autosomal.bim
indivname:      kaz12_autosomal_1.fam
outputformat:   EIGENSTRAT
genotypeoutname:        kaz12.geno
snpoutname:     kaz12.snp
indivoutname:   kaz12.ind
```
```bash
convertf -p convertf_param.par
```

```bash
nano merge_param.par
```
```bash
geno1:  dataset/v54.1.p1_1240K_public.geno
snp1:  dataset/v54.1.p1_1240K_public.snp
ind1:  dataset/v54.1.p1_1240K_public.ind
geno2:  kaz12.geno
snp2:  kaz12.snp
ind2:  kaz12.ind
genooutfilename:  merged.geno
snpoutfilename:  merged.snp
indoutfilename:  merged.ind
outputformat:  EIGENSTRAT
```

```bash
mergeit -p merge_param.par
```

Now need to construct left and right files and run the analysis!
