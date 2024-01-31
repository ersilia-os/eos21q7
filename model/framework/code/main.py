# imports
import os
import csv
import joblib
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import numpy as np
import pandas as pd
from rdkit.Chem import QED


# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# checkpoints directory
checkpoints_dir = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))

# read checkpoints (here, simply an integer number: 42)
model = joblib.load(os.path.join(checkpoints_dir, "random_forest_model.pkl"))


# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# Convert SMILES to fingerprints
mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024) for m in mols]

# Convert fingerprints to numpy array
arr_list = []
for fp in fps:
    arr = np.zeros((0,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    arr_list.append(arr)

test_x = np.stack([arr.tolist() for arr in arr_list])
test_x = test_x.astype(np.float32)

# run model
outputs = model.predict_proba(test_x)[:, 1]  
# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"])  # header
    for o in outputs:
        writer.writerow([o])