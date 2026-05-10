#!/bin/bash

############配置tfcomb环境 ##########
##### python = 3.10
#mamba create -n tfcomb --file required_packages.txt   ##在github下载

#mamba install pandas=1.5.2   #pandas 2.1版本remove append

eval "$(conda shell.bash hook)"

conda activate seq


bedtools instersect 



eval "$(conda shell.bash hook)"
conda activate tfcomb

jupyter notebook


