"""
import the mice networks into a dictionary that allows for easier manipulation. 

Final Structure (dict):

    {
        'LPS':{
            'before' : [net1, net2, ...]
            'after' : [net1, net2, ...]
        }, 
        'Control':{
            'before' : [net1, net2, ...]
            'after' : [net1, net2, ...]
        }, 
        'no_injection: {
            'before' : [net1, net2, ...]
            'after' : [net1, net2, ...]
        }
    }

NOTE: The 'before injection' network is matched with it's corresponding 'after injection network by index. 
"""
from os.path import join
import pandas as pd
import networkx as nx

from mice_cleaning import print_partitioned_list

def get_mice_networks():
    """Import the mice encounters as an edgelist, partitioned into the components. 

    """
    filename = "partitioned_dryad_encounters.csv"
    parent_dir = join("data", "mice")
    df = pd.read_csv(join(parent_dir, filename))

    all_nets = {}
    # For all tratement levels 
    for treatment in df['treatment'].unique():
        all_nets[treatment] = {}
        
        # For 'before' and 'after' injection
        for timepoint in df['timepoint'].unique():
            all_nets[treatment][timepoint] = []
            condition_subset = df.loc[(df['treatment'] == treatment) & (df['timepoint'] == timepoint)]
            
            # For all the components in that condition.
            for c in condition_subset['component_num'].unique():
                # Subset down to an edglist for that component
                component_subset = condition_subset.loc[
                    condition_subset['component_num'] == c].drop(
                        ['treatment', 'timepoint', 'component_num'], axis = 1)
                
                # Append network from pandas dataframe into master dictionary 
                all_nets[treatment][timepoint].append(nx.from_pandas_edgelist(component_subset, 
                source= "first_mouse", target="second_mouse", edge_attr=True))
    
    print_partitioned_list(all_nets)
    

if __name__ == "__main__":
    get_mice_networks()