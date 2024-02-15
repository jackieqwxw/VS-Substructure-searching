#!/usr/bin/env python3
# coding: utf-8

# In[1]:

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools

def RemoveDuplicatesFromLists(li):
    import itertools
    li.sort()
    res = list(k for k, _ in itertools.groupby(li))
    return res


# In[2]:


import sys

inputfile = sys.argv[1]
databasefile = sys.argv[2]
outputfile1 = sys.argv[3]
outputfile2 = sys.argv[4]


# In[3]:


# read fragments
suppl = Chem.SDMolSupplier(inputfile)
mol_fm = [x for x in suppl if x is not None]
mols_fm = [Chem.MolToSmiles(m) for m in mol_fm]
legends_frag = [x.GetProp("_Name") for x in suppl]
print(f'Numbers of fragments:{len(mol_fm)}, they are {legends_frag}')


# In[4]:


# read database
read_smiles = Chem.SmilesMolSupplier(databasefile,titleLine=False)
mol_db = [x for x in read_smiles if x is not None]
print('Numbers of Smiles in database:',len(mol_db))
mols_db = [Chem.MolToSmiles(m) for m in mol_db]
legends_db = [x.GetProp("_Name") for x in read_smiles]
dbzip = list(zip(mols_db,legends_db))


# In[5]:


db_smi = []
db_id = []
match_id = []

for smiles, legend in dbzip:
    molecule = Chem.MolFromSmiles(smiles)

    substructure_1 = Chem.MolFromSmiles(mols_fm[0])
    substructure_2 = Chem.MolFromSmiles(mols_fm[1])
    substructure_3 = Chem.MolFromSmiles(mols_fm[2])


    has_match_1_2 = molecule.HasSubstructMatch(substructure_1, useChirality=True) and molecule.HasSubstructMatch(substructure_2, useChirality=True)
    has_match_1_3 = molecule.HasSubstructMatch(substructure_1, useChirality=True) and molecule.HasSubstructMatch(substructure_3, useChirality=True)

    if has_match_1_2:
        db_smi.append(smiles)
        db_id.append(legend)
        match_id.append(legends_frag[0]+'-'+legends_frag[1])       

    if has_match_1_3:
        db_smi.append(smiles)
        db_id.append(legend)
        match_id.append(legends_frag[0]+'-'+legends_frag[2])

match_smi = list(zip(db_smi,db_id))
match_smi_rmdup = RemoveDuplicatesFromLists(match_smi)

match_smi_frag = list(zip(db_smi,db_id,match_id))
match_smi_frag_rmdup = RemoveDuplicatesFromLists(match_smi_frag)


# In[6]:

## Headers
print("##############################################################")
print("#  Substructure searching for virtual screening\n#")
print("#  Computational Biophysics & Drug Design (Kireev lab)\n#")
print("#  Developed by Xiaowen Wang")
print("##############################################################\n")

writer = Chem.SmilesWriter(outputfile1)
for i,j in enumerate(match_smi_rmdup[:]):
    smi, label = match_smi_rmdup[i]
    mol = Chem.MolFromSmiles(smi)
    mol.SetProp('_Name',label)
    writer.write(mol)
writer.close()
print('Hit numbers of substructure matches after removing duplicates):', writer.NumMols())

writer = Chem.SmilesWriter(outputfile2)
writer.SetProps(['NAME_Match'])
for i,j in enumerate(match_smi_frag_rmdup[:]):
    smi, label, frag = match_smi_frag_rmdup[i]
    mol = Chem.MolFromSmiles(smi)
    mol.SetProp('_Name',label)
    mol.SetProp('NAME_Match',frag)
    writer.write(mol)
writer.close()
