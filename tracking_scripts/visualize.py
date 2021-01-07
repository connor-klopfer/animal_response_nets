import ant_tracking as a_t
from matplotlib import pyplot
import numpy as np

################
#
# Goal: Make an animation of the ants over time.
#
# Manuchehr Aminian (maminian@cpp.edu)
#

timegroups = a_t.df.groupby('t')

times_dict = timegroups.groups
times = list(times_dict.keys())	# list of time stamps

freq = 1/0.5    # readme says sampling rate is half of a second
plot_rate = 60  # seconds

######
all_ant_ids = np.unique(a_t.df['id'].values)   # how many ants?
colors = pyplot.cm.rainbow(np.linspace(0,1,len(all_ant_ids)))

ant_cmap = {aid: colors[i] for i,aid in enumerate(all_ant_ids) }

box = [-100,3200,-100,4500]   # x0,x1,y0,y1

fig,ax = pyplot.subplots(1,1, 
    constrained_layout=True,
    figsize=(9.2,6.6)
)

t = -999
pic_idx = 0
for ti in times:
    if ti >= t + plot_rate*freq:
        t = ti
        
        dfs = timegroups.get_group(t)
        x,y = dfs['x'].values, dfs['y'].values
        u,v = np.cos(dfs['theta'].values), np.sin(dfs['theta'].values)
        colors = [ant_cmap[ai] for ai in dfs['id'].values]
        ant_ids = dfs['id'].values
        
        ax.cla()
        ax.quiver( x,y,u,v , color=colors)

        ax.axis('square')
        ax.set_xlim(box[0],box[1])
        ax.set_ylim(box[2],box[3])
        
        ax.set_title('Timestamp: %i'%t)
        fig.savefig('frames/ants_colony_%s.png'%str(pic_idx).zfill(5))
        pic_idx += 1
#