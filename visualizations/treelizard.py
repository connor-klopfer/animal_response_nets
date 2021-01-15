import pandas
import networkx
import numpy as np
from matplotlib import pyplot
import os

#from vis_params import rcParams
import vis_params

#PREFIX = '../analysis_results/'

df = pandas.read_csv(vis_params.ANALYSIS_FOLDER + 'edgelist_treelizards_0weightexcl.csv')

#df = df[['trt', 'id1', 'id2', 'distance', 'weight']]

elist = df[['id1','id2']].values
#eweights = df['distance'].values
eweights = df['weight'].values

elist2 = [ (e[0],e[1],w) for e,w in zip(elist,eweights) ]

G = networkx.DiGraph()

G.add_weighted_edges_from(elist2)
#G = networkx.from_edgelist(elist2)

fig,ax = pyplot.subplots(1,1, constrained_layout=True)

pos_xy = networkx.spring_layout(G)
networkx.draw_networkx(G,
#        arrowsize=eweights,
        pos = pos_xy,
        font_color=[1,1,1],
        font_size=8,
        edge_color=[0,0,0,0.4],
        node_shape='s',
        node_size=400
)

fig.savefig(vis_params.IMAGES_FOLDER + 'treelizard.png')
fig.show()