import os
import csv
import joblib
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem import QED
from sklearn.preprocessing import StandardScaler

# Define ECFP calculation function
def ecfp_calculation(smiles):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    bit_info_list = []
    bit_info = {}
    fps = []
    for mol in mols:
        fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=1024, bitInfo=bit_info))
        bit_info_list.append(bit_info.copy())
    arr_list = [np.zeros((0,), dtype=np.int8) for _ in range(len(fps))]
    for i, bit in enumerate(fps):
        DataStructs.ConvertToNumpyArray(bit, arr_list[i])
    ecfps = np.stack([i.tolist() for i in arr_list])
    df = pd.DataFrame(ecfps.astype(np.float32))
    return df

# Define QED calculation function
def qed_calculation(smiles):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    qe = [QED.properties(mol) for mol in mols]
    qe_df = pd.DataFrame(qe)
    return qe_df

# Define function to scale QED descriptors
def scale_qed(df):
    scaler = StandardScaler()
    scaled_array = scaler.fit_transform(df)
    scaled_df = pd.DataFrame(scaled_array, columns=df.columns)
    return scaled_df

# Define function to preprocess the data
def preprocess_data(df):
    # Detect the column containing SMILES strings
    smiles_column = df.columns[df.apply(lambda col: col.apply(lambda x: Chem.MolFromSmiles(str(x)) is not None)).all()][0]
    smiles = df[smiles_column].tolist()
    
    # Calculate ECFP descriptors
    ecfp_df = ecfp_calculation(smiles)
    
    # Calculate QED descriptors
    qed_df = qed_calculation(smiles)
    
    # Scale the QED descriptors
    scaled_qed_df = scale_qed(qed_df)
    
    # Concatenate the descriptor dataframes
    processed_df = pd.concat([ecfp_df, scaled_qed_df], axis=1)
    
    return processed_df

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# checkpoints directory
checkpoints_dir = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))

# read checkpoints (here, simply an integer number: 42)
model = joblib.load(os.path.join(checkpoints_dir, "random_forest_model.pkl"))

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# read input data
input_data = pd.read_csv(input_file)

# preprocess data
processed_data = preprocess_data(input_data)

# convert the processed data to numpy array
test_x_combined = processed_data.values.astype(np.float32)

# run model
outputs = model.predict_proba(test_x_combined)[:, 1]

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"])  # header
    for o in outputs:
        writer.writerow([o])
