import pandas
import matplotlib
from matplotlib import pyplot
import numpy as np
import networkx

import vis_params  # local parameters file.

import reorder  # node reordering, for e.g. circular layout.
#

treatment_palette = {
                    'sick': [0,0.5,0],
                    'control': [0.5,0,0.5]
                }

#

df = pandas.read_csv(vis_params.ANALYSIS_FOLDER + 'bats_pre_post_edgelist.csv')

all_bats = np.unique( df[['bat1','bat2']].values )
n = len(all_bats)

cond_dict = dict( list( df.groupby('condition') ) )
ncond = len(cond_dict)

# assign a node color for every bat by treatment.
treatment_col_dict = {b: treatment_palette['control'] for b in all_bats}
for row in df.iloc:
    treatment_col_dict[ row['bat1'] ] = treatment_palette[ row['treatment_bat1'] ]
#


####
# make the figure.
fig,ax = pyplot.subplots( 1,ncond, figsize=(10,5), constrained_layout=True )

for i,c in enumerate(['pre', 'post']):
    df_s = cond_dict[c]
    
    elist = df_s[['bat1','bat2','weight']].values
    G = networkx.Graph()
    G.add_weighted_edges_from(elist)
    
    # reorder for some sort of clustering
    G = reorder.reorder_nx_Graph(G)
    
    edge_widths = [G[u][v]['weight']/1000 for u,v in G.edges()]
    
    if c=='pre':
        pos_xy = networkx.spectral_layout(G)
        pos_xy = networkx.spring_layout(G, pos=pos_xy, iterations=1, k=0.4)
    else:
        # try manual adjustment
        # many iterations will get to a steady state - but 
        # no real "message" about how network adjusts after treatment.
#        pos_xy = networkx.spectral_layout(G)
        pos_xy = networkx.spring_layout(G, pos=pos_xy, iterations=1, k=0.4)
    #
    
    node_colors = [treatment_col_dict[b] for b in G.nodes()]
    
    pyplot.sca( ax[i] )
    networkx.draw_networkx( G, 
        node_color = node_colors,
        pos = pos_xy,
        arrowsize=edge_widths,
        font_color='w',
        font_size=12,
        edge_color=[0,0,0,0.2],
        node_shape='s')
    
    ax[i].set_title(c)
#

# make a legend to indicate color meaning
for k,v in treatment_palette.items():
    ax[0].scatter([],[], c=[v], label=k, s=200)
ax[0].legend(loc='best')


if __name__=="__main__":
    fig.show()
    fig.savefig(vis_params.IMAGES_FOLDER + 'bats.png')