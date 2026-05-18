#!/bin/bash
#BSUB -J run-liana
#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n 16
#BSUB -R "rusage[mem=8192]"       # 16 cores × 8192 MB = ~128G
#BSUB -R "span[hosts=1]"
#BSUB -W 4:00
#BSUB -o logs/run-liana_%J.out
#BSUB -e logs/run-liana_%J.err
#BSUB -L /bin/bash

module load R/4.4.1
module load glpk/4.55

proj_dir=$(pwd)
mkdir -p logs

cd ${proj_dir}/lrp

export RENV_PATHS_PREFIX=$(cat .renv_platform | cut -d'/' -f1)

bash_cmd="Rscript ${proj_dir}/lrp/scripts/main.R"

echo $bash_cmd
($bash_cmd)
