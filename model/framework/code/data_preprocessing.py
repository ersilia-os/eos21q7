import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem import QED
from sklearn.preprocessing import StandardScaler

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

    def scale_qed(self, df):
        scaler = StandardScaler()
        scaled_array = scaler.fit_transform(df)
        scaled_df = pd.DataFrame(scaled_array, columns=df.columns)
        return scaled_df

    def preprocess_data(self, df):
        # Detect the column containing SMILES strings
        smiles_column = df.columns[df.apply(lambda col: col.apply(lambda x: Chem.MolFromSmiles(str(x)) is not None)).all()][0]
        smiles = df[smiles_column].tolist()
        
        # Calculate ECFP descriptors
        ecfp_df = self.ecfp_calculation(smiles)
        
        # Calculate QED descriptors
        qed_df = self.qed_calculation(smiles)
        
        # Scale the QED descriptors
        scaled_qed_df = self.scale_qed(qed_df)
        
        # Concatenate the descriptor dataframes
        processed_df = pd.concat([ecfp_df, scaled_qed_df], axis=1)
        
        return processed_df
