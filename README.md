# VS-Substructure-searching
This tool is used for substructure searching for virtural screening

## Workflow
Install RDKit toolkit on Linux using conda:
```
conda create --name rdkit-tools python=3.8
```
```
conda activate rdkit-tools
```
```
conda install -c conda-forge rdkit
```

## Execute the substructure searching for one file
Run the script:
```
python ./substructure_search.py fragments.sdf $dbinp ${file}.smi ${file}.csv
```
*`fragments.sdf` stores the substructure fragments used for searching. The variable `$dbinp` is SMILES file name of database used for substructure searching. The `$file` denotes the outputfile name. This script is considered as an example to search the database file with the fragments of adenine + sulfone or adenine + tetrazole combinations. It will generate two types of files. One file is ${file}.smi with the format `SMILES ID` and the other file is ${file}.csv with the format `SMILES ID NAMEOFMATCHSUBSTRUCTURE`. 


## Execute the substructure searching in batch
If you perfer to perform substructure searching for a large ligand collection. First, split database into multiple small files. Then use the bash script `substructure_search.sh` to perform the calculations in bath.
```
split -l 10000 -d -e --additionl-suffix=.smi input.smi outprefixhere
```
```
ls | wc -l
```
```
sbatch substructure_search.sh
```

* Note that change the `sbatch parameter setting`, `$path`, and `fragments.sdf` variables accordingly.  
