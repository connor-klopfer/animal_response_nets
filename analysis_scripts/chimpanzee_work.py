import pandas
import numpy as np
import networkx
from matplotlib import pyplot
import os

np.random.seed(0)

PREFIX = '../data/Chimpanzees_HobaiterEtAl2014/'
PREFIX = os.path.abspath(PREFIX)

########

ef = pandas.ExcelFile(PREFIX + '/' + 'Data for online storage.xlsx')
df = ef.parse('Less Strict')

def get_datetime(row):
    import datetime
    import re
    pattern = '[A-Z]{2,}\:\ ?([0-9]{1,})\:([0-9]{1,})'
    hour,minute = re.match( pattern, row['Time action'] ).groups()
    day,month,year = row['Date'].split('.')
    dt = datetime.datetime(year=int(year)+2000, month=int(month), day=int(day), hour=int(hour), minute=int(minute))
    return dt
#

df['datetime'] = [get_datetime(row) for row in df.iloc]

df = df.sort_values('datetime')

events_flat = []

non_unique_actors = []
for i,row in enumerate( df.iloc ):
    m = row['MD']
    os = []
    non_unique_actors.append( m )
    
    # extract all chimpanzees that observed the event.
    for j in range(8,df.shape[1]):
        if isinstance(row.iloc[j], str):
            os.append( row.iloc[j][:2] )
            events_flat.append( [m,row.iloc[j][:2], row['Day'], row['Novel technique']] )
        else:
            # missing string value means it's the end of the list of chimpanzees observing.
            # We also want to record the chimpanzee exhibiting the event to "itself"
            # for the purposes of tracking non-observed events.
            events_flat.append( [m,m, row['Day'], row['Novel technique']] )
            break
    # add to the record of all chimpanzees.
    non_unique_actors += os
#

df_events = pandas.DataFrame(data=events_flat, columns=['Mediator','Observer','Day','Novel technique'])

# remove duplicates of chimpanzee list.
actors = np.unique( non_unique_actors )

# try 1: create graph across all points in time.
# also construct dictionaries assigning chimpanzee names to integers.
A = np.zeros( (len(actors), len(actors)) )
adict = {a:i for i,a in enumerate(actors)}
invdict = {i:a for a,i in adict.items()}

for row in df_events.values:
    A[adict[row[0]],adict[row[1]]] += 1

# build networkx directed graph object.
G = networkx.DiGraph(A)
G = networkx.relabel.relabel_nodes(G, invdict)
edge_weights = np.array( [G[u][v]['weight'] for u,v in G.edges()] )

# pick a layout... none here are 'perfect'.
#
#nodes_xy = networkx.kamada_kawai_layout(G)
#nodes_xy = networkx.circular_layout(G)
nodes_xy = networkx.shell_layout(G)
#nodes_xy = networkx.spectral_layout(G)

if __name__ == "__main__":
    # visualize.
    fig = pyplot.figure()
    networkx.draw_networkx(G, with_labels=True, font_color='w' , width=edge_weights, node_size=200, pos=nodes_xy)

    pyplot.ion()
    pyplot.show(block=False)