"""
The mice data is in a weird format, so this reformats the networks into a cleaner version
"""
from os.path import join
import datetime
import pandas as pd
import networkx as nx


def import_mouse_encounters(parent_dir = None):
    dataset_directory = join("data", "mice") if parent_dir is None else parent_dir
    encounters = pd.read_csv(join(dataset_directory, "dryad_encounters.csv"))
    encounters['date_obs_dt'] = pd.to_datetime(encounters['date_obs'], format = "%d-%b")
    return encounters

def networks_by_date():
    encounters = import_mouse_encounters()
    date_levels = get_date_levels(encounters, "date_obs_dt")

    all_nets = {}
    for date in date_levels:
        encounters_subset = encounters.loc[encounters['date_obs'] == date]
        all_nets[date] = nx.from_pandas_edgelist(encounters_subset,
            source = "IDA[inside_first]", 
            target = "IDB[inside_after_IDA]", edge_attr=True, 
            create_using=nx.DiGraph)
    return all_nets

def get_date_levels(df, date_col):
    dates = df[date_col].unique()
    dates.sort()
    return list(pd.to_datetime(dates).strftime("%d-%b"))


def get_injection_data(parent_dir = None):
    dataset_directory = join("data", "mice") if parent_dir is None else parent_dir
    injections = pd.read_csv(join(dataset_directory, "dryad_injected_animal_info.csv"))
    injections['date_injected_dt'] = pd.to_datetime(injections['date_injected'], format = "%d-%b")
    return injections


def partition_network_components():
    net_list = networks_by_date()
    injections = get_injection_data()
    injection_dates = get_date_levels(injections, "date_injected_dt")
    encounters = import_mouse_encounters()
    sample_dates = get_date_levels(encounters, "date_obs_dt")
    partitioned_nets = {}

    treatment_levels = list(injections['injection'].unique())
    treatment_levels.append('no_injection')

    for t in treatment_levels:
        partitioned_nets[t] = {
            'before': [],
            'after' : []
        }
    
    subgraphs = get_network_subgraphs(net_list)
    for date in injection_dates:
        
        previous_date = sample_dates[sample_dates.index(date) - 1]
        injection_subset = injections.loc[injections['date_injected'] == date]
        for _, mouse in injection_subset.iterrows():
            injected_mouse_net, subgraphs[date] = get_subgraph_w_mouse(mouse['ID'], subgraphs[date])
            before_injection_mouse_net, subgraphs[previous_date] = get_subgraph_w_mouse(mouse['ID'], subgraphs[previous_date])
            
            if injected_mouse_net is None or before_injection_mouse_net is None:
                continue
            
            if len(injected_mouse_net.nodes()) < 10 or len(before_injection_mouse_net.nodes()) < 10:
                continue

            # Append to the injection lists 
            partitioned_nets[mouse['injection']]['after'].append(injected_mouse_net)
            partitioned_nets[mouse['injection']]['before'].append(before_injection_mouse_net) 
        
        # partition remaining mice for the rest of the networks. 
        for net in subgraphs[date]:
            if len(net.nodes()) >= 10:
                for prev_net in subgraphs[previous_date]:
                    if len(set(net.nodes()).intersection(set(prev_net.nodes()))) > 8:
                        partitioned_nets['no_injection']['after'].append(net.copy())
                        partitioned_nets['no_injection']['before'].append(prev_net.copy())
        
    print_partitioned_list(partitioned_nets)
    write_partitions_to_file(partitioned_nets)


def get_subgraph_w_mouse(mouse_id, subgraphs):
    mouse_bool_array = [mouse_id in sub.nodes() for sub in subgraphs]
    # mouse is not in network. 
    if sum(mouse_bool_array) == 0:
        return None, subgraphs
    mouse_net_index = mouse_bool_array.index(True)
    mouse_net = subgraphs[mouse_net_index].copy()
    return mouse_net, [s.copy() for idx, s in enumerate(subgraphs) if idx != mouse_net_index]


def print_partitioned_list(partition_list):
    for treatment in partition_list:
        for time in partition_list[treatment]:
            print("Treatment {}:Time {}: N items: {}".format(treatment, time, len(partition_list[treatment][time])))
            for idx, n in enumerate(partition_list[treatment][time]):
                print("{}: | N Edges: {} | N Nodes: {}: | N Intersection of Nodes: {}".format(
                    idx, len(n.edges()), len(n.nodes()), 
                    len(set(n.nodes()).intersection(set(partition_list[treatment]['before' if time == 'after' else 'after'][idx].nodes())))))
    

def get_network_subgraphs(parent_networks):
    sub_graphs = {}
    for net_date, net in parent_networks.items():
        sub_graphs[net_date] = [net.subgraph(c).copy() for c in nx.strongly_connected_components(net)]
    return sub_graphs

def write_partitions_to_file(partitioned_networks):
    all_datasets = []
    # Write the parititoned to file. 
    parent_dir = join("data", "mice")
    for treatment in partitioned_networks:
        for timepoint in partitioned_networks[treatment]:
            for idx, n in enumerate(partitioned_networks[treatment][timepoint]):
                # filename = "{}_{}_{}.csv".format(treatment, timepoint, idx)
                df = nx.to_pandas_edgelist(n, source = 'first_mouse', target = "second_mouse")

                df['treatment'] = treatment
                df['timepoint'] = timepoint
                df['component_num'] = idx
                all_datasets.append(df)
            
                # df.to_csv(join(parent_dir, filename))

    final = pd.concat(all_datasets)
    final.to_csv(join(parent_dir, "partitioned_dryad_encounters.csv"), index = False)

if __name__ == "__main__":
    partition_network_components()
