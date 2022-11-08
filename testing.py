#from load_data import data_preprocessing
import networkx as nx
import pandas as pd
import numpy as np
import os

from results_processing import results_analysis
from libAP import LSOprimizer

#true_case = open("data/true_HT.txt", "r").read().split("\n")[:-1]
#true_control = open("data/true_control.txt", "r").read().split("\n")[:-1]

#path_expr, path_net = 'Sample.txt', 'IID_LUNG.graphml'
#GE, GX, labels_ids, rev = data_preprocessing(path_expr, path_net, log2=True, size=3000)

cwd= os.getcwd()

# load Gene expression 
df = pd.read_csv('~/clustering/Network-Clustering-timeseries/Sample.txt', sep='\t')
#df.head()
GE = df.set_index('Geneid')
#pat = np.array(list(df.columns))
#GE= df.transpose()
#GE= df.to_numpy().reshape(71,2,58051)

#print(pat,pat.ndim)

#n, m = GE.shape

G = '/home/rna/clustering/Network-Clustering-timeseries//Network.graphml'
T = 20
alpha = 0.01  # speed of temperature decrease
L_min = 150
L_max = 200

optimizer = LSOprimizer(GE, G, T, L_min, L_max, max_iter=20)
#print(optimizer)
nodes, labels, sc = optimizer.run_ls()
print(nodes,labels,sc)