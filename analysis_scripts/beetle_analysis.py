#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""
Analysis of a nested dictionary of NetworkX graphs to measure network metrics: 

    1.Average density, Number of nodes, Mean Betweenness Centrality, and Mean Eigenvector Centrality 
    2.Node Degree centrality, Node Eigenvector Centrality, Node Closeness Centrality, and Node Betweenness Centrality
    
Metrics are compiled into groups (1,2) then exported into a .csv file
"""


def get_graph_metrics(network):
    k3_dict=[]
    append_network_mean=[]
    append_network_cent=[]
    for k1,v1 in network.items():
        for k2,v2 in v1.items():
            v2=v2[0]
            #Average network density values
            density   = nx.density(v2)

            #Number of nodes
            tot_nodes = nx.number_of_nodes(v2)

            #Total average betweenness values for the entire network
            ##Listing all individual betweenness values, adding them up and dividing by the full 
            tot_betweenness=0
            betweenness=(nx.betweenness_centrality(v2))
            for k, v in list(betweenness.items()): 
                tot_betweenness=v+tot_betweenness
            mean_betweenness = (tot_betweenness/tot_nodes)

            #Total Eigenvector centrality values for the entire network
            tot_eigenvector=0
            eigenvector=(nx.eigenvector_centrality(v2))
            for k, v in list(eigenvector.items()): 
                tot_eigenvector=v+tot_eigenvector
            mean_eigengenvector = (tot_eigenvector/tot_nodes)


            #take measurements for centrality values and insert into a dictionary under Pandas
            df2 = pd.DataFrame(dict(
            DEGREE_CENTRALITY      = nx.degree_centrality(v2),
            EIGENVECTOR            = nx.eigenvector_centrality(v2),
            CLOSENESS_CENTRALITY   = nx.closeness_centrality(v2),
            BETWEENNESS_CENTRALITY = nx.betweenness_centrality(v2)))

            #compile mean values into dictionary
            d={
            'TOT_NODES':tot_nodes,
            'MEAN_EIGENVECTOR'       : mean_eigengenvector,
            'MEAN_BETWEENNESS'       : mean_betweenness}
            k3=k1,k2
            k3_dict.append(k3)
            df3=pd.DataFrame(data=d,index=[k3])
            append_network_mean.append(df3)
            append_network_cent.append(df2)
    append_network_mean=pd.concat(append_network_mean)
    append_network_cent=pd.concat(append_network_cent,keys=k3_dict)
    #print(append_network_mean)
    print(append_network_cent)
            #print(k1,k2,df2)
            #print(k1,k2,df3)
            
            #move to a .csv file
    append_network_cent.to_csv('centrality_measures.csv')
    append_network_mean.to_csv('mean_measures.csv')


# In[ ]:




