import os
import csv
import joblib
import sys
from data_preprocessing import DataProcessor
import pandas as pd
import numpy as np

# Define paths
root = os.path.dirname(os.path.abspath(__file__))
checkpoints_dir = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))

# Load model
model = joblib.load(os.path.join(checkpoints_dir, "random_forest_model.pkl"))

# Read input file
input_file = sys.argv[1]
output_file = sys.argv[2]
input_data = pd.read_csv(input_file)

# Preprocess data
processor = DataProcessor()
processed_data = processor.preprocess_data(input_data)

# Convert processed data to numpy array
test_x_combined = processed_data.values.astype(np.float32)

# Run model
outputs = model.predict_proba(test_x_combined)[:, 1]

#check input and output have the same length
input_len = len(input_data)
output_len = len(outputs)
assert input_len == output_len

# Write output to CSV
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"])  # header
    for o in outputs:
        writer.writerow([o])
