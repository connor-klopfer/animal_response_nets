import pandas
import networkx
import numpy as np
from matplotlib import pyplot
import os

PREFIX = '../analysis_results/'

df = pandas.read_csv(PREFIX + 'edgelist_treelizards_0weightexcl.csv')

#df = df[['trt', 'id1', 'id2', 'distance', 'weight']]

elist = df[['id1','id2']].values
#eweights = df['distance'].values
eweights = df['weight'].values

elist2 = [ (e[0],e[1],w) for e,w in zip(elist,eweights) ]

G = networkx.DiGraph()

G.add_weighted_edges_from(elist2)
#G = networkx.from_edgelist(elist2)

fig,ax = pyplot.subplots(1,1)
networkx.draw_networkx(G)
