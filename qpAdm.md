Directory: ~/biostar/gwas/redo_july/working_with_ref_data/ancient/

sed 's/-9/1/g' kaz12_autosomal.fam > kaz12_autosomal_1.fam

nano convertf_param.par 
genotypename:   kaz12_autosomal.bed
snpname:        kaz12_autosomal.bim
indivname:      kaz12_autosomal_1.fam
outputformat:   EIGENSTRAT
genotypeoutname:        kaz12.geno
snpoutname:     kaz12.snp
indivoutname:   kaz12.ind

convertf -p convertf_param.par

nano merge_param.par
geno1:  v54.1.p1_1240K_public.geno
snp1:  v54.1.p1_1240K_public.snp
ind1:  v54.1.p1_1240K_public.ind
geno2:  kaz12.geno
snp2:  kaz12.snp
ind2:  kaz12.ind
genooutfilename:  merged.geno
snpoutfilename:  merged.snp
indoutfilename:  merged.ind
outputformat:  EIGENSTRAT

mergeit -p merge_param.par

I stopped here! (takes time to run)

nano parqpfstat.txt
DIR: ~biostar/gwas/redo_july_working_with_ref_data/ancient/bin
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

