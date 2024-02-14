import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem import QED
from sklearn.preprocessing import StandardScaler
import pickle
import os

# Define paths
root = os.path.dirname(os.path.abspath(__file__))
checkpoints_dir = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))

class DataProcessor:
    def __init__(self):
        pass

    def ecfp_calculation(self, smiles):
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

    def qed_calculation(self, smiles):
        mols = [Chem.MolFromSmiles(smi) for smi in smiles]
        qe = [QED.properties(mol) for mol in mols]
        qe_df = pd.DataFrame(qe)
        return qe_df

    # Function to load scaler from a pickle file
    def load_scaler(self, scaler_filepath):
        with open(scaler_filepath, 'rb') as f:
            scaler = pickle.load(f)
        return scaler

        
    def preprocess_data(self, df):
        # Detect the column containing SMILES strings
        smiles_column = df.columns[df.apply(lambda col: col.apply(lambda x: Chem.MolFromSmiles(str(x)) is not None)).all()][0]
        smiles = df[smiles_column].tolist()
        
        # Calculate ECFP descriptors
        ecfp_df = self.ecfp_calculation(smiles)
        
        # Calculate QED descriptors
        qed_df = self.qed_calculation(smiles)
        
        # Load scaler from pickle file
        scaler = self.load_scaler(scaler_filepath)

        # Scale QED descriptors using the loaded scaler
        scaled_array = scaler.transform(qed_df)
        scaled_qed_df = pd.DataFrame(scaled_array, columns=qed_df.columns)
        
        # Concatenate the descriptor dataframes
        processed_df = pd.concat([ecfp_df, scaled_qed_df], axis=1)
        
        return processed_df



# Define scaler file path
scaler_filepath = os.path.join(checkpoints_dir, "scaler.pkl")
