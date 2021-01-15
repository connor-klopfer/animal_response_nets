import numpy as np
import networkx
from matplotlib import pyplot
#,rcParams
#rcParams['axes.titlesize'] = 18
#rcParams['axes.labelsize'] = 18


import vis_params  # local parameters file.


import chimpanzee_work as cw


c_involved = np.array([0.2,0.4,0.6,1.0])
c_not_involved = np.array([0.2,0.4,0.6,0.2])

c_learned_dict = {'moss':[0.2,0.6,0.4,1.0],
                            're-use1':[0.6,0.2,0.6,1.0]}


# arbitrary edge-case cast as moss.
# Only one event has this label.
cw.df_events.loc[cw.df_events['Novel technique']=='moss/re-use1' , 'Novel technique'] = 'moss'

techniques = np.unique( cw.df_events['Novel technique'].values )
days = np.unique( cw.df_events['Day'].values )


fig,ax = pyplot.subplots(len(techniques),len(days), 
    figsize=(3*len(days),3*len(techniques)), 
    constrained_layout=True,
    sharex=True,
    sharey=True)


for kk,c in enumerate( techniques ):
        actives = cw.df_events.loc[cw.df_events['Novel technique']==c, ['Mediator','Observer']].values
        actives = np.unique(actives)
        learned = []
        for ll,d in enumerate( days ):
            eventmask = np.logical_and( cw.df_events['Novel technique'].values==c, cw.df_events['Day'].values==d )
            eventmask = np.where(eventmask)[0]
            
            edge_sub = cw.df_events.iloc[eventmask,:2].values
#            print(c, d)
#            print(edge_sub)
            learned = np.unique( np.concatenate( [ learned, edge_sub[:,0] ] ) )
            
            # assign color...
            colors = []
            for node in cw.G.nodes():
                # has the chimpanzee exhibited the task?
                if node in learned: 
                    colors.append( c_learned_dict[c] )
                else:
                   # will the chimpanzee ever observe/be observed?
                    if node in actives:
                        colors.append( c_involved )
                    else: 
                        colors.append( c_not_involved )
            #
            
            pyplot.sca( ax[kk,ll] )
            networkx.draw_networkx( cw.G, with_labels=True, font_color='w' , 
                width=cw.edge_weights, 
                node_size=200, 
                pos = cw.nodes_xy,
                edgelist = edge_sub,
                font_size=8,
                node_color=colors
                )

            if kk==0:
                ax[kk,ll].set_title('Day %i'%(ll+1))
            if ll==0:
                ax[kk,ll].set_ylabel(c)
#

fig.savefig(vis_params.IMAGES_FOLDER + 'chimpanzee_bytechnique_byday.png')
fig.show()