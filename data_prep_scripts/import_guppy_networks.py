"""
Import guppy association networks that are cleaned into a newtowkxobject, and the node attributes 
associated with the other files. 

By: Connor Klopfer 09Jan21

There are 16 groups of high risk and low risk each.

Imported file filetype is a dictionary, one key for high risk, low risk, and a list for all the groups within each risk
group.

all_networks = {
    low_risk: [net, net, net, ...],
    high_risk: [net, net, net, ...]
}
"""
import os
import pandas as pd
import networkx as nx

def import_individual_data():
    """Import individual data from 'Individual Level Data.csv' from the 
    repository. 

    Returns:
        Pandas DataFreame: Pandas Dataframe, in the same format as the original file. 
    """
    filepath = os.path.join("data", "guppies", "Individual Level Data.csv")
    guppy_ind = pd.read_csv(filepath)
    return guppy_ind

def import_assoc_networks():
    """Import clean association networks as networkx graph object. All graphs
    are in a dictionary. 

    Returns:
        dict: Each key is 'high_risk' and 'low_risk' which maps to a list of networkx 
        Graph objects. Each list should be 16 graphs long. 
    """
    # Parent dir of the clean network edgelist files. 
    parent_dir = os.path.join("data", "guppies", "reformatted_networks")
    all_nets = os.listdir(parent_dir)
    # Dictionary has two groups 
    guppy_networks = {
        'low_risk': [None] * 16,
        'high_risk' : [None] * 16
    }
    # Iterate though all the edgelist files in the parent directory 
    for n in all_nets:
        n_split = n.split("_")
        # Group Number (1-16), Index is Group Number - 1, (because of 0 indexing)
        group_n = int(n_split[1]) - 1
        group_risk = "{}_{}".format(n_split[2], n_split[3])
        # Read in file
        guppy_edgelist = pd.read_csv(os.path.join(parent_dir, n))
        # Convert to Networkx from directory. 
        guppy_networks[group_risk][group_n] = nx.convert_matrix.from_pandas_edgelist(
            guppy_edgelist, source = "guppy1", target = 'guppy2', edge_attr=True)
        
    return guppy_networks

def import_assoc_networks_w_attrs():
    """Import the associations from the clean files, and add on the node attributes. 

    Returns:
        dict: Dictionary as returned by import_assoc_networks() but wth node attributes. 
    """
    nets_no_attr = import_assoc_networks()
    nets_w_attr = append_guppy_ind_attr(nets_no_attr)
    return nets_w_attr


def append_guppy_ind_attr(all_networks):
    """Add node attributes for the guppy association networks, from the 
    'Individual Level Data.csv' file. 
    NOTE: Only include IndividualID and body length, but easy to add more. 

    Args:
        all_networks (dict): Dictionary containing list of networks as returned by import_assoc_networks()

    Returns:
        dict: Same format as passed parameter except nodes now have attributes. 
    """
    ind_data = import_individual_data()
    for risk in all_networks:
        risk_label = "H" if risk == "high_risk" else "L"
        for n_idx, n in enumerate(all_networks[risk]):
            n_attrs = {}
            for node in n.nodes():
                ind_id = "{}{}{}".format(risk_label, n_idx + 1, node)
                ind_length = ind_data.loc[ind_data['IndividualID'] == ind_id, 'Length_cm'].values[0]
                n_attrs[node] = {'ID' : ind_id, 'length': ind_length}
                nx.set_node_attributes(n, n_attrs)
        
    return all_networks


if __name__ == "__main__":
    # import_individual_data()
    # import_assoc_networks()
    import_assoc_networks_w_attrs()