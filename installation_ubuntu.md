
```bash
curl http://data.biostarhandbook.com/install/bash_profile.txt >> ~/.bash_profile
curl http://data.biostarhandbook.com/install/bashrc.txt >> ~/.bashrc
source ~/.bash_profile

URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
curl $URL > miniconda-installer.sh
bash miniconda-installer.sh -b
~/miniconda3/condabin/conda init

conda -V
conda update -y -n base conda
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -y --name bioinfo python=3.6
conda activate bioinfo
curl http://data.biostarhandbook.com/install/conda.txt | xargs conda install -y

efetch -db nuccore -id 2 -format gb
mkdir -p ~/bin
curl http://data.biostarhandbook.com/install/doctor.py > ~/bin/doctor.py
chmod +x ~/bin/doctor.py
~/bin/doctor.py
```

I downloaded ubuntu version of rstudio and then ran
```bash
sudo dpkg -i rstudio-2024.04.2-764-amd64.deb 
root /usr/lib/rstudio/chrome-sandbox
sudo chmod 4755 /usr/lib/rstudio/chrome-sandbox
```

sudo apt-get update
sudo apt-get install -y liblzma-dev libcurl4-openssl-dev libssl-dev libxml2-dev build-essential
sudo apt-get install libbz2-dev

$ sudo cp plink2 /usr/local/bin/
$ sudo cp plink /usr/local/bin/

sudo mkdir -p /mnt/windows/biostar/alzheimer/idats/
sudo cp -r /home/user/biostar/gwas/alzheimer/idats/ /mnt/windows/biostar/alzheimer/idats/

sudo mkdir -p /mnt/windows/biostar/alzheimer/manifest/
sudo cp -r /home/user/biostar/gwas/alzheimer/manifest/ /mnt/windows/biostar/alzheimer/manifest/

sudo mkdir -p /mnt/windows/biostar/tools/genome_studio/
sudo cp -r /home/user/tools/genome_studio/ /mnt/windows/biostar/tools/genome_studio/ 

sudo mkdir -p /mnt/windows/biostar/alzheimer/idats_clean/
sudo find /mnt/windows/biostar/alzheimer/idats/ -type f -name "*.idat" -exec cp {} /mnt/windows/biostar/alzheimer/idats_clean/ \;

Checking how many unique samples i have:
ls /mnt/windows/biostar/alzheimer/idats_clean/*_*.idat | sed 's/_Red.idat//;s/_Grn.idat//' | sort | uniq | wc -l
Answer = 288
