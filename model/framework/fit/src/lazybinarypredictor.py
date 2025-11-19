import os
import sys
import pandas as pd
import numpy as np
import json
from lazyqsar.qsar import LazyBinaryQSAR
from lazyqsar.utils.logging import logger
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
import collections
import shutil

root = os.path.dirname(os.path.abspath(__file__))

test_file = os.path.join(root, "..","data", "test.csv")
model_folder = os.path.join(root, "..", "results", "validation_model")

test = pd.read_csv(test_file)
x_test = test["st_smiles"].astype(str).tolist()
y_test = test["bin"].astype(int).tolist()


def evaluate(x_test, y_test, output_dir):
    model = LazyBinaryQSAR.load(output_dir)
    y_pred = model.predict_proba(smiles_list=x_test)[:, 1]
    logger.info("ROC-AUC: {0}".format(roc_auc_score(y_test, y_pred)))
    return roc_auc_score(y_test, y_pred), y_pred

roc_score, y_hat = evaluate(x_test, y_test)

report = {}
report["y_true"] = list(y_test)
report["y_hat"] = list(y_hat)
report["roc_auc"] = roc_score

with open(os.path.join(root, "..", "results",f"report_validation.json"), "w") as f:
    json.dump(report, f, indent=2)