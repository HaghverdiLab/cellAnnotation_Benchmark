We usimport ItClust as ic
import scanpy as sc
import os
from numpy.random import seed
from tensorflow.compat.v1 import set_random_seed
import pandas as pd
import numpy as np
import warnings
import sys

print("Starting...")

os.environ["CUDA_VISIBLE_DEVICES"]="1"
warnings.filterwarnings("ignore")

args = sys.argv
trainFolder = args[1]
testFolder = args[2] # + args[1] + "/"
os.chdir(args[3])

#Set seeds
seed(20180806)
np.random.seed(10)
set_random_seed (20180806) # on GPU may be some other default

# Source data
adata_train = sc.read_csv((trainFolder+"data_train.txt"),
                           delimiter=",", first_column_names=True)
celltypes = pd.read_csv((trainFolder +"meta_train.txt"), sep=',')
adata_train.obs["celltype"] = celltypes.class_.values
adata_train.obs_names_make_unique(join="-")
adata_train.var_names_make_unique(join="-")
adata_test = sc.read_csv((testFolder + "data_test.txt"),
                          delimiter=",",first_column_names=True)
adata_test.var_names = adata_train.var_names
adata_test.obs_names_make_unique(join="-")
adata_test.var_names_make_unique(join="-")

# Fit ItClust model
clf=ic.transfer_learning_clf()
clf.fit(adata_train, adata_test)

# Prediction
pred, prob, cell_type_pred=clf.predict()
pred.head()
