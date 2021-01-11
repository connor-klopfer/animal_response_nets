"""
Taking the association networks from the guppy data and saving to a friendlier format type, for easier use.
"""
import os
import csv

def reformat_raw_association_networks(directory):
    """Import the raw association matrices from the original file 
    on the repository. The current file is not very easy to import. 

    Args:
        directory (string): Path file as returned by os.path.join() where 
        original file is stored. 
    """
    filename = os.path.join(directory, "Association Networks.csv")
    # Store in a dictionary
    all_networks = {}
    # Open original file 
    with open(filename, "r") as datafile:
        fin = csv.reader(datafile, delimiter = ",")
        for row in fin:
            # IF "Risk" is in the row, it's a header containing node labels
            if "Risk" in row[0]:
                # Group Number (1-16)
                group_id = row[0].split("Risk")[-1]
                # Get Node labels for this network
                low_risk_header, high_risk_header = process_header(row)
                # In the original file, like-group numbers are on the same row. Index 
                # forst by group number, then partition by risk. 
                all_networks[group_id] = {
                    'low_risk':[],
                    'high_risk': []
                }
            # IF the first element in the row is '', it's a spacer row. 
            elif row[0] != '':  
                # Data Row
                temp_lr_edgelist, temp_hr_edgelist = process_non_header(low_risk_header, high_risk_header, row)  
                # Add to edge list of the network.
                all_networks[group_id]['low_risk'].extend(temp_lr_edgelist)
                all_networks[group_id]['high_risk'].extend(temp_hr_edgelist)

    # Archive to file
    archive_edgelist(all_networks)
               

def process_header(header_row):
    """Get the NODE IDs from header rows in the association networks. 

    Args:
        header_row (list): List of strings, where the first element is the group (Risk + Number) and 
        the following elements are the node IDs. 

    Returns:
        list, list: list of node labels (strings of numbers) for the low risk and high risk groups respectively. 
    """
    # Theres and empty cell between the high risk and low risk groups 
    break_idx = ["HighRisk" in x for x in header_row].index(True) - 1
    # Low risk nodes 
    low_risk = [x for x in header_row[1:break_idx] if x != '']
    # High Risk
    high_risk = [x for x in header_row[break_idx + 2:] if x != '']
    return low_risk, high_risk

def process_non_header(low_risk_header, high_risk_header, row):
    """Process row containing the edge weights in the association network. 

    Args:
        low_risk_header (list): list of the node labels in the low risk group, each 
            item is a string of a number 
        high_risk_header (list): list of the node labels in the high risk group, each 
            item is a string of a number
        row (list): List of association counts, the first element in the list is the node label. 

    Returns:
        list[list]: A list of lists, where each element in the list is a list with the source, target and edge weight. 
    """
    
    low_risk_weights = row[1:(1+len(low_risk_header))]
    high_risk_weights = row[11:11 + (len(high_risk_header))]
    # The first item in the matrix is the node label of the individual 
    lr_ind = row[0]
    hr_ind = row[10]

    low_risk_edgelist = []
    high_risk_edgelist = []
    # Append the edge list for low risk groups 
    if lr_ind != '':
        for idx, lr_h in enumerate(low_risk_header):
            low_risk_edgelist.append([lr_ind, lr_h, low_risk_weights[idx]])
    # For high risk groups 
    if hr_ind != '':
        for idx, hr_h in enumerate(high_risk_header):
            high_risk_edgelist.append([hr_ind, hr_h, high_risk_weights[idx]])
    # NOTE: Had to do this seperately, because there might be didfferent populations sizes 
    # between high risk and low risk groups 
    return low_risk_edgelist, high_risk_edgelist
    
def archive_edgelist(networks):
    """Write the edgelist to file. Every network for every risk groups gets 
    it's own file. The format is group_(group number)_(risk level)_risk_group.csv

    Args:
        networks (dict): a dictionary of dictionaries, the keys map to another 
        dictionary, which then has a 'high_risk' and 'low_risk' key, which maps to a list 
        of lists, which represents the edgelist. 
    """
    filename = os.path.join("data", "guppies", "reformatted_guppy_networks.csv")
    with open(filename, 'w', newline="\n") as fout:
        net_tofile = csv.writer(fout, delimiter = ",")
        net_tofile.writerow(['guppy1', "guppy2", 'weight', 'group', 'risk']) 
        
        # For each group. 
        for e in networks:
            # For high and low risk
            for risk_group in networks[e]:
                    # Header for file
                    # Edgelist to file. 
                    for r in networks[e][risk_group]:
                        net_tofile.writerow([*r, e, risk_group])
                    # net_tofile.writerows(networks[e][risk_group])

if __name__ == "__main__":
    reformat_raw_association_networks(os.path.join("data", "guppies"))
