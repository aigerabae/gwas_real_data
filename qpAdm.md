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

```bash
nano parqpfstat.txt
```
```bash
DIR: ~/biostar/gwas/redo_july_working_with_ref_data/ancient/
S1: 10Oct21
S1X: 10Oct21
indivname: merged.ind
snpname: merged.snp
genotypename: merged.geno
poplistname: ./bin/fstat_23andme/lista.txt
fstatsoutname: ./bin/fstat_23andme/fstatsa.txt
allsnps: YES
inbreed: NO
scale: NO
```

```bash
nano ./bin/fstat_23andme/lista.txt
```
```bash
Sohi
Mbuti.DG
Irula.DG
Turkey_N
Laos_LN_BA.SG
China_Tianyuan
Ami.DG
Karitiana.DG
Iran_GanjDareh_N
Iran_C_SehGabi
Iran_ShahrISokhta_BA1
Iran_ShahrISokhta_BA2
Turkmenistan_Gonur_BA_1
Dai.DG
Russia_Ust_Ishim_HG.DG
Chukchi.DG
Saami.DG
Georgia_Kotias.SG
Russia_Kostenki14.SG
Russia_Tyumen_HG
Russia_MLBA_Sintashta
Russia_DevilsCave_N.SG
Luxembourg_Loschbour.DG
Czech_BellBeaker
Kazakhstan_Central_Saka.SG
Portugal_MN.SG
Kazakhstan_Andronovo.SG
```

```bash
qpfstats -p parqpfstat.txt > bin/fstat_23andme/qpfstatlog.txt
```

I stopped here!
```bash
qpfstats -p parqpfstat.txt > bin/fstat_ancestry/qpfstatlog.txt
```

```bash
nano parqpadm.txt
```
```bash
fstatsname: ~/biostar/gwas/redo_july_working_with_ref_data/ancient/bin/fstat_ancestry/fstatsa.txt
popleft: ~/biostar/gwas/redo_july_working_with_ref_data/ancient/bin/fstat_ancestry/left.txt
popright: ~/biostar/gwas/redo_july_working_with_ref_data/ancient/bin/fstat_ancestry/right.txt
details: YES
```

```bash
nano left.txt
```
```bash
Sohi
Iran_ShahrISokhta_BA2
Kazakhstan_Andronovo.SG
Turkmenistan_Gonur_BA_1
```

```bash
nano right.txt
```
```bash
Mbuti.DG
China_Tianyuan
Karitiana.DG
Russia_Ust_Ishim_HG.DG
Ami.DG
Dai.DG
Turkey_N
Georgia_Kotias.SG
Russia_Kostenki14.SG
Iran_GanjDareh_N
```

```bash
qpAdm -p parqpadm.txt > sohi_qpadm_output.txt
```
