#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 11:53:32 2018

@author: floriankreten
"""


import plotly.offline as py
import plotly.graph_objs as go

from pathlib import Path

import EVAL_data_output

""" Preferably use this module via EVAL_single_graph_evaluation

    A graph G can be plotted in 3D using those functions
    Graphs can be plotted in detail with up to ~100.000 elements
    
    for larger graphs, see subset_plot
    
    Produces htmls that can be opened in any (?) browser
"""


def offline_plot(G, file_name, fullpathname = False, show = True):
    
    """A 3D-Plot of a given graph G"""
    
    size_factor = 5
    #positions
    xV = [ v.position[0] for v in G.V[1:] ]
    yV = [ v.position[1] for v in G.V[1:] ]
    zV = [ v.position[2] for v in G.V[1:] ]
    
    # vertex-labels, below are shown 2 examples
    #labels = [ "T_10 " + str ( v.trait [-10:]) + " | F " + str( round(v.fitness,2) )
    #         for v in G.V[1:] ]
    labels = [ "F " + str( round(v.fitness, 2) ) for v in G.V[1:] ]    
    sizes = [v.radial_size * size_factor * 1.8 for v in G.V[1:]]
    colours = colour_by_trait(G)
    
    #DO NOT TOUCH THE FOLLOWING LINES!!!
    vertices = go.Scatter3d(x=xV,
        y=yV,
        z=zV,
        text=labels,
        hoverinfo = 'text',     #remove this line to show also the coordinates
        mode = 'markers',
        marker=dict(
            sizemode='diameter',
            size = sizes,
            color = colours,
            colorscale = 'Jet',
            cmin = 0,
            cmax = 1,
            opacity = 0.1
        )
    )
        
    dots = go.Scatter3d(x=xV,   #is a debug for opacity =1, plotly error(!)
        y=yV,
        z=zV,
        hoverinfo = 'none',
        mode = 'markers',
        marker=dict(
            sizemode='diameter',
            size = sizes,
            color = colours,
            colorscale = 'Jet',
            cmin = 0,
            cmax = 1,
            opacity = 1
        )
    )
    
    xE = []
    yE = []
    zE = []
    
    for v in G.V[2:]:
        xE += [ v.position[0], G[v.in_edge].position[0], None   ]
        yE += [ v.position[1], G[v.in_edge].position[1], None   ]
        zE += [ v.position[2], G[v.in_edge].position[2], None   ]
    
    
    edges=go.Scatter3d(x=xE,
            y=yE,
            z=zE,
            mode='lines',
            line=dict(color='rgb(125,125,125)', width = min ( min(3, size_factor/3.0), 10 ),
                      ),
            hoverinfo='none'
            )
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.5, y=1.5, z=1.5) )
        
    noaxis=dict(visible = False,
                showbackground=False,
              showline=False,
              zeroline=False,
              showgrid=False,
              showticklabels=False
              )

    layout = dict(height=900,
                       width=900,
                       showlegend=False,
                 scene=dict(
                 xaxis=noaxis,
                 yaxis=noaxis,
                 zaxis=noaxis,
                 camera=camera
                     ),
                 hovermode = 'closest',
                 shapes=[
                     dict(
                    type="circle",
                    xref="x",
                    yref="y",
                    x0=10,
                    y0=20,
                    x1=10,
                    y1=20,
                    fillcolor="blue",
                    line_color="blue"
                ) ]
            )
    
    
    data = [vertices, edges, dots]
    
    figure = go.Figure(data=data, layout=layout)
    
    # Creates directory for storage
    if not fullpathname:
        dirName = EVAL_data_output.create_folder(file_name)
    else:
        dirName = file_name
        Path(dirName).mkdir(parents=True, exist_ok=True)
    #print("Spatial plot in:", dirName)
    print("SPATIAL PLOT filename", dirName)
    name = str(dirName + "/spatial_plot.html")
    py.plot(figure, filename = name, auto_open=show )
    return
    
def subset_plot(G, subset_in, file_name, show = True, color_key = False,
                title = None, fullpathname = False, biopsy_colour_code = False):
    
    """ Plots a subset of G (without edges)
        the rest of the graph if sketched by grey dots """
    
    
    """ Plot of the given subset """
    size = 5
    subset = list ( subset_in )
    good_size = 55000
    if len(subset) > good_size:
        print("Subset for plot rescaled")
        K = int( len(subset) / good_size ) +1
        subset = subset[::K]
    
    #positions
    xV = [ v.position[0] for v in subset ]
    yV = [ v.position[1] for v in subset ]
    zV = [ v.position[2] for v in subset ]
    
    if biopsy_colour_code:
        #print("plotting biopsy with colour_coding")
        num_biopsies = len ( set ( v.in_biopsy for v in subset_in)    )
        color_vals = [ 1.0 / num_biopsies * i for i in range(num_biopsies) ]
        colors = [ color_vals[ v.in_biopsy -1]  for v in subset ]
    elif not color_key:
        colors = "blue"
    else:
        colors = [ color_key[v.cluster_group] for v in subset ]
    
    sam = go.Scatter3d(x=xV,
        y=yV,
        z=zV,
        hoverinfo = 'none',
        mode = 'markers',
        marker=dict(
            sizemode='diameter',
            size = size,
            colorscale = 'Jet',
            color = colors,
            cmin = 0,
            cmax = 1,
            opacity = 1
        )
    )
    
    
    """ Plot of the grey mass """
    good_size = 10000
    S = len(G.V[1:])
    if S > good_size:
        K = int( S / good_size ) +1
        V = G.V[1:][::K]
      
    xV = [ v.position[0] for v in V ]
    yV = [ v.position[1] for v in V ]
    zV = [ v.position[2] for v in V ]
    
    grey_mass = go.Scatter3d(x=xV,
        y=yV,
        z=zV,
        hoverinfo = 'none',
        mode = 'markers',
        marker=dict(
            sizemode='diameter',
            size = size,
            #color = colours,
            #colorscale = 'Jet',
            color = "grey",
            cmin = 0,
            cmax = 1,
            opacity = 0.02
        )
    )
    
    """ Setup of the plot """
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.5, y=1.5, z=1.5) )
        
    noaxis=dict(visible = False,
                showbackground=False,
              showline=False,
              zeroline=False,
              showgrid=False,
              showticklabels=False
              )

    layout = dict(height=900,
                       width=900,
                       showlegend=False,
                 scene=dict(
                 xaxis=noaxis,
                 yaxis=noaxis,
                 zaxis=noaxis,
                 camera=camera
                     ),
                 hovermode = 'closest',
                 title=title
                 )
    
    
    data = [sam, grey_mass]
    
    figure = go.Figure(data=data, layout=layout)
    
    # Creates directory for storage
    if not fullpathname:
        dirName = EVAL_data_output.create_folder(file_name)
    else:
        dirName = file_name
        Path(dirName).mkdir(parents=True, exist_ok=True)
    
    if title:
        name = str(dirName + "/" + title + ".html")
    else:
        name = str(dirName + "/subset_spatial_plot.html")
    py.plot(figure, filename = name, auto_open=show )
    
    return

#%%

def colour_by_trait(G):
    """returns colours for the vertices of a graph"""
    actives =  dict() #stores all living traits
    for v in G.V[1:]:
        actives[v.trait] = False
    colour_map = [ float(1/len(actives)) * i for i in range(len(actives)) ]
    
    #map colours to traits and vertices
    colours = []
    k=0
    for t in actives:
        actives[t] = colour_map[k]
        k += 1
    for v in G.V[1:]:
        colours.append(actives[v.trait])
    return colours