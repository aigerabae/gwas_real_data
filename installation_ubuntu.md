
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
```
