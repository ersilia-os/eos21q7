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
model = joblib.load(os.path.join(checkpoints_dir, "Random_forest_model.pkl"))

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# read SMILES from input .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# Convert SMILES to molecules
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

# If smiles don't transform to mol, add to none_list
none_list = []
for i in range(len(mols)):
    if mols[i] is None:
        none_list.append(i)
        print('add to none_list')

reg_idx = 0
for i in none_list:
    del mols[i - reg_idx]
    reg_idx += 1

# Modify index
if len(none_list) != 0:
    for i in none_list:
        del smiles_list[i - reg_idx]
        reg_idx += 1

# Create fingerprint
bit_info_list = []  # Bit vector
bit_info = {}  # Bit vector
fps = []
b = 0

# Mol to fingerprint Bit Vector
for a in mols:
    fps.append(AllChem.GetMorganFingerprintAsBitVect(a, 3, nBits=1024, bitInfo=bit_info))
    bit_info_list.append(bit_info.copy())

# To array
arr_list = list()
for i in range(len(fps)):
    array = np.zeros((0,), dtype=np.int8)
    arr_list.append(array)

for i in range(len(fps)):
    bit = fps[i]
    DataStructs.ConvertToNumpyArray(bit, arr_list[i])

test_x = np.stack([i.tolist() for i in arr_list])
test_x = test_x.astype(np.float32)
test_finprt = pd.DataFrame(test_x)

# Create physicochemical properties
qe = [QED.properties(mol) for mol in mols]
qe = pd.DataFrame(qe)

# Merge QED properties dataframe with train_finprt
test_finprt = pd.concat([test_finprt, qe], axis=1)

# Convert the test_finprt dataframe to numpy array
test_x_combined = test_finprt.values.astype(np.float32)

# run model
outputs = model.predict_proba(test_x_combined)[:, 1]

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"])  # header
    for o in outputs:
        writer.writerow([o])
