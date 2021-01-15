from scipy import spatial,sparse,cluster
import numpy as np
import networkx

def get_compressed(mat):
    n = len(mat)
    return np.array( [mat[i,j] for i in range(1,n) for j in range(i+1,n)] )
#

def get_permutation(mat):
    uppertri = get_compressed(mat)
#    Z = cluster.hierarchy.median(uppertri)
    Z = cluster.hierarchy.centroid(uppertri)
    o = cluster.hierarchy.leaves_list(Z)
    return o
#

def nx_graph_to_adj(nx_G):
    nodelist_orig = nx_G.nodes()
    edgelist = nx_G.edges()
    edgelist = np.array(edgelist)
    n2i = {n:i for i,n in enumerate(nodelist_orig)}
    n = len(nodelist_orig)
    A = np.zeros( (n,n) )
    
    weights = [nx_G[u][v].get('weight',1) for u,v in edgelist]
    weights = np.reshape(weights, (len(edgelist),1) )
    edgelist = np.hstack( [edgelist, weights] )
    for u,v,w in edgelist:
        A[n2i[u],n2i[v]] = w
        A[n2i[v],n2i[u]] = w    # do we want this?
    return A
#
    

def reorder_nx_Graph(nx_G):
    nodelist_orig = np.array( nx_G.nodes() )
    A_sym = nx_graph_to_adj(nx_G)   # beware: symmetric even if G is directed.
    o = get_permutation(A_sym)
    nodelist_new = nodelist_orig[o]
    
    G2 = networkx.Graph()
    G2.add_nodes_from(nodelist_new)

    edgelist = nx_G.edges()
    edgelist = np.array(edgelist)
    weights = [nx_G[u][v].get('weight',1) for u,v in edgelist]
    weights = np.reshape(weights, (len(edgelist),1) )
#    edgelist = np.hstack( [edgelist, weights] )
    edgelist = [(e[0],e[1],{'weight':w}) for e,w in zip(edgelist,weights)]
    G2.add_edges_from( edgelist )
    return G2
#