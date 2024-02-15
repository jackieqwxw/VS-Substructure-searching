#!/bin/bash

#SBATCH --account kireevlab
#SBATCH --job-name 29p2_Substructure
#SBATCH --time 4:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=50GB
#SBATCH --array 0-1010
#SBATCH --exclude=lewis4-r630-hpc4-node233

module load miniconda3
eval "$(conda shell.bash hook)"
conda activate rdkit-tools

start=$(date +%s)

path=/home/xwpnp/group/database/ENR_5.5B/filtering/CBDD/remove_duplicates/ENR29p2

arr=(`ls $path/*`)
dbinp=${arr[$SLURM_ARRAY_TASK_ID]}
file=`(echo $dbinp | awk -F '/' '{print $NF}' | awk -F '.' '{print$1}')`
echo $dbinp

python ./substructure_search.py PLCg1-fragments-ATP.sdf $dbinp ${file}.smi ${file}.csv

end=$(date +%s)
take=$(( end - start ))
min=`echo "scale=2; $take / 60" | bc`
echo Time taken to execute commands is ${min} minites.

