import vis_params
from matplotlib import pyplot
from sklearn import cluster, linear_model
import numpy as np
import pandas

from mpl_toolkits import mplot3d

def summarize(metadata,cluster_labels):
    animal_list = np.unique(metadata['animal'].values)
    animal_counts = {k:sum(metadata['animal']==k) for k in animal_list}
    
    clusters = np.array(list(cluster_labels.keys()))
    cluster_sizes = [len(cluster_labels[c]) for c in clusters]
    presentation_order = np.argsort(cluster_sizes)[::-1]
    clusters = clusters[ presentation_order ]
    
    for k in clusters:
        v = cluster_labels[k]
        #
        
        print('Cluster %s (n=%i)'%(str(k),len(v)))
        print('================')
        
        subs = metadata['animal'].iloc[v].values
        subs_unique = np.unique(subs)
        
        cluster_belong = np.zeros(len(subs_unique))
        rel_belong = np.zeros(len(subs_unique))
#        print(len(subs_unique))
        
        for j,s in enumerate(subs_unique):
            cluster_belong[j] = sum(subs==s)/len(v)
            rel_belong[j] = sum(subs==s)/animal_counts[s]
        #
        ord = np.argsort(cluster_belong)
        ord = ord[::-1] # decreasing order
        for o in ord:
#            print(subs_unique[o],type(subs_unique[o]))
#            print(cluster_belong[o],type(cluster_belong[o]))
#            print(rel_belong[o],type(rel_belong[o]))
            print('\t%s: %.2f%% of cluster \t %.2f%% of animal.'%(str(subs_unique[o]).ljust(15),100*cluster_belong[o],100*rel_belong[o]))
        #
        print('')
    #
    
#
        
    

df = pandas.read_csv( vis_params.ANALYSIS_FOLDER + 'all_individual_results.csv')

data = df[['eigenvector_centrality', 'betweeness', 'degree']].values

meta = df[['animal','condition']]

data[data[:,1]==0,1] = 0.1*np.min(data[data[:,1]!=0,1])
data[:,1] = np.log10(data[:,1])
data[:,2] = np.log10(data[:,2])
#data = np.log10(data)

#

dbs = cluster.DBSCAN(eps=0.15)

cluster_labels = dbs.fit_predict(data)

classes = {k:np.where(cluster_labels==k)[0] for k in np.unique(cluster_labels)}

#

fig = pyplot.figure(constrained_layout=True, figsize=(8,8))
ax = fig.add_subplot( projection='3d' )

clusters = np.array(list(classes.keys()))
cluster_sizes = [len(classes[c]) for c in clusters]
presentation_order = np.argsort(cluster_sizes)[::-1]
clusters = clusters[ presentation_order ]


#for i,(k,v) in enumerate( classes.items() ):
for i,k in enumerate(clusters):
    v = classes[k]
    
    if k==-1:
        continue
    ax.scatter(data[v,0],data[v,1],data[v,2], c=[pyplot.cm.tab10(i%20)], label='Cluster %s'%str(k))
    if i==9:
        break
#

v = classes[-1]
ax.scatter(data[v,0],data[v,1],data[v,2], c=[[0.8,0.8,0.8,0.5]], label='Noise cluster')

ax.set_xlabel('e-vec centrality')
ax.set_ylabel(r'$\log_{10}$(betweenness)')
ax.set_zlabel(r'$\log_{10}$(degree)')
ax.legend(loc='upper left')

fig.savefig('individual_clustering_dbscan.png')
fig.show()
pyplot.ion()

summarize(meta,classes)
