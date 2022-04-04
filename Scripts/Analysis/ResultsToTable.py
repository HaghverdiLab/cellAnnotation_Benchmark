#!/usr/bin/python
import sys
import os
import numpy as np
import pandas as pd 

def getBase(file, method, length):
    file = file.rstrip('/')
    base=os.path.basename(file)
    print(base)
    base = os.path.splitext(base)[0]
    if method == "Seurat":
        base= base.replace("txt","")
    base= base.split(sep="_")
    
    #if len(base)==length:
    #    base = base + ["None"]
    print(base)
    return base
    
args = sys.argv
file = args[1]
output = args[2]
method = args[3]
print("Starting....")
if method in ["Seurat", "SCN", "CellID", "SingleR"]:
    if "prediction" in file:
        sys.exit()
    base = getBase(file, method,2)
    X = pd.read_csv(file, sep="\t", header=0)
    X = X.dropna()
    print(X)
    y_true = X["class_"]
    y_pred = X["predicted"]
    out = [method] + base
    print(out)
elif method == "ItClust":
    base = getBase(file, method,2)
    file = file + "/results.txt"
    X = pd.read_csv(file, sep="\t", header=0)
    y_true = X["class_"]
    y_pred = X["predicted_celltype"]
    out = ["ItClust"] + base
elif method == "MLP":
    base = getBase(file, method,2)
    file = file + "/predictions.csv"
    X = pd.read_csv(file, sep=",", header=None)
    y_true = X[0]
    y_pred = X[1]
    out = ["MLP"] + base
else:
    sys.exit("The method has to be either scPotter, Seurat or ItClust")
print("Get Accuracies")
class_accuracies = []
print(np.unique(y_true))
print("-----------")
print(np.unique(y_pred))
print("-----------")

for class_ in np.unique(y_true):
    print(class_)
    class_acc = np.mean(y_pred[y_true == class_] == class_)
    class_accuracies.append(class_acc)
print("------------------------------------------")
accuracy=np.mean(y_pred == y_true)
print(accuracy)

out_types = out + class_accuracies + [accuracy]

out_types = ",".join(str(e) for e in out_types)
print(out_types)   
print(output)
with open(output, 'a') as file:
    file.write(out_types +'\n')

