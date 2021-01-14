#!/usr/bin/env python
# coding: utf-8

"""
Taking the association networks from the beetle data and saving to a friendlier format type, for easier use.

Used VThorsethief's code as a template
"""


from os.path import join
import pandas as pd
import networkx as nx

def get_beetle_networks():
    """Import the beetle encounters as an edgelist, partitioned into the components. 
    """
    filename = "SRI_edgelist.csv"
    df = pd.read_csv(join(filename))

    all_nets = {}
    # For all treatement levels 
    for treatment in df['Group.ID.Populations'].unique():
        all_nets[treatment] = {}
        
        # For 'before' and 'after' injection
        for treatment_period in df['Treatment.Period'].unique():
            all_nets[treatment][treatment_period] = []
            
            treatment_subset = df.loc[(df['Treatment.Period'] == treatment_period) & (df['Group.ID.Populations'] == treatment)]
            # Append network from pandas dataframe into master dictionary 
            all_nets[treatment][treatment_period].append(nx.from_pandas_edgelist(treatment_subset,"Focal.Letter","Social.Partner", edge_attr=True))
 
    return all_nets






