#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 13:13:24 2018

@author: floriankreten
"""

"""" EXPERIMENTAL """

""" The animate_plot - function takes a result from 
    PROCESS_process.over_time_simulation (see there for the usage)
    and produces a plot with timeslider for visualization over time 
    Only use this one for small fun proof-of-concept-plots as it is very memory-intensive"""  

import plotly.offline as py
import plotly.graph_objs as go

import EVAL_data_output
import PROCESS_process

import random

#%%

def sim_and_plot_with_size_frames(frame_sizes="default", filename = "predefined_3D_run"):
    """ Simulate and 3D plot at certain defined sizes """
    if frame_sizes == "default":
        frame_sizes = [i*1000 for i in range(11) ]

    assert ( max(frame_sizes)<= 100000  )
    G, grid = PROCESS_process.over_time_simulation(
        frame_sizes, spatial_size=300, key="size")
    animate_plot(G, frame_sizes, grid,
                    file_name = filename )
    
    print("Sim and time-plot successful")


#%%

def define_camera_and_layout(fig, axis_range, sliders):
    
    noaxis=dict(visible = False,
              showline=False,
              zeroline=False,
              showticklabels=False,
              range = [-axis_range,axis_range]
              )
    
    fig.update_layout(sliders=sliders)
    fig.update_layout(scene = dict(
                        xaxis = noaxis,
                        yaxis = noaxis,
                        zaxis = noaxis
                      ))
    
    return fig

def add_data_to_grid( G, frame_id, grid ):
    """ In: Graph G, timepoint time, pandaFrame grid
        Adds data from the graph to the grid
        is done during the simulation"""

    
    time = frame_id
    name_template = '{time}+{header}'

    # trait and fitness    
    labels = []
    for v in G.V[1:]:
        labels.append ( "Trait: " + str ( v.trait ) + " | Fitness " + str( round(v.fitness,2) ) )
    
    traits = [ v.trait for v in G.V[1:] ]
    
    key_labels = name_template.format ( time= time, header = "labels" )
    key_traits = name_template.format ( time= time, header = "traits" )
    
    grid = grid.append ( {'value': labels, 'key': key_labels}, ignore_index = True )
    grid = grid.append ( {'value': traits, 'key': key_traits}, ignore_index = True )
    
    # positions of vertices
    xV = [ v.position[0] for v in G.V[1:] ]
    yV = [ v.position[1] for v in G.V[1:] ]
    zV = [ v.position[2] for v in G.V[1:] ]
    
    key_xV = name_template.format ( time = time, header = "xV" )
    key_yV = name_template.format ( time = time, header = "yV" )
    key_zV = name_template.format ( time = time, header = "zV" )
    
    grid = grid.append ( {'value': xV, 'key': key_xV}, ignore_index = True )
    grid = grid.append ( {'value': yV, 'key': key_yV}, ignore_index = True )
    grid = grid.append ( {'value': zV, 'key': key_zV}, ignore_index = True )
    
    # positions of edges
    xE = []
    yE = []
    zE = []
    
    for v in G.V[2:]:
        xE += [ v.position[0], G[v.in_edge].position[0], None   ]
        yE += [ v.position[1], G[v.in_edge].position[1], None   ]
        zE += [ v.position[2], G[v.in_edge].position[2], None   ]
    
    
    key_xE = name_template.format ( time = time, header = "xE" )
    key_yE = name_template.format ( time = time, header = "yE" )
    key_zE = name_template.format ( time = time, header = "zE" )
    
    grid = grid.append ( {'value': xE, 'key': key_xE}, ignore_index = True )
    grid = grid.append ( {'value': yE, 'key': key_yE}, ignore_index = True )
    grid = grid.append ( {'value': zE, 'key': key_zE}, ignore_index = True )
    
    
    return grid


def insert_colours(G, times, grid):
    """ inserts the colours to the plot-grid according to the traits
        can only be done after the simulation is finished"""
    
    #set up colout-mapping
    length = len( G.fitness_book)
    colour_map = [ 1/length * i for i in range(length) ]
    random.shuffle(colour_map)

    
    Iterator = list ( G.fitness_book.keys() ) #list of all types
    Colour_map = {}     # injective colour mapping via a dictionary
    for i in range( len(Iterator) ):
        key = Iterator[i]
        Colour_map[key] = colour_map[i]
        
    
    name_template = '{time}+{header}'
    
    #insert colours into the grid
    def update_grid(grid, time):
        key_traits = name_template.format ( time= t, header = "traits" )
        try:
            traits = grid.loc[grid['key']== key_traits, 'value'].values[0]
        except IndexError:
            raise IndexError("Boss, I dont know what happened...")
        colours = []
        for trait in traits:
            colours.append( Colour_map[trait] )
            
        key_colours = name_template.format ( time= t, header = "colours" )
        grid = grid.append ( {'value': colours, 'key': key_colours}, ignore_index = True )
        return grid
    
    for t in times:
        grid = update_grid(grid, t)
    
    print("All colors defined")
    
    return grid
    

def extract_plot_data_from_grid(grid, time, point_size):
    """ In: Grid,  time-point
        Out: Plot-Data
        Edges are deactivated at the moment"""
    
    name_template = '{time}+{header}'
    
    key_xV = name_template.format ( time = time, header = "xV" )
    key_yV = name_template.format ( time = time, header = "yV" )
    key_zV = name_template.format ( time = time, header = "zV" )
    

    key_xE = name_template.format ( time = time, header = "xE" )
    key_yE = name_template.format ( time = time, header = "yE" )
    key_zE = name_template.format ( time = time, header = "zE" )

    key_labels  = name_template.format ( time = time, header = "labels"  )
    key_colours = name_template.format ( time = time, header = "colours" )
    
    vertices = go.Scatter3d(
            x=grid.loc[grid['key']== key_xV, 'value'].values[0],
            y=grid.loc[grid['key']== key_yV, 'value'].values[0],
            z=grid.loc[grid['key']== key_zV, 'value'].values[0],
            
        text=grid.loc[grid['key']== key_labels, 'value'].values[0],
        hoverinfo = 'text',
        mode = 'markers',
        marker=dict(
            sizemode='diameter',
            size = point_size,
            color = grid.loc[grid['key']== key_colours, 'value'].values[0],
            colorscale = 'Jet',
            cmin = 0,
            cmax = 1,
            opacity = 1
        ),
        showlegend = False
    )
    
    edges=go.Scatter3d(
            x=grid.loc[grid['key']== key_xE, 'value'].values[0],
            y=grid.loc[grid['key']== key_yE, 'value'].values[0],
            z=grid.loc[grid['key']== key_zE, 'value'].values[0],
            mode='lines',
            line=dict( color='grey', width = min (point_size/3, 3)  ),
            hoverinfo='none',
            showlegend = False
            )
    
    return vertices, edges


def animate_plot(G, times, grid, file_name = "sim_example" ):
    """ In: Graph G, time-array times, information-grid 
        Prints given Graph/grid over the timescale
        
        axis_range and point_size can be adjusted
        WARNING: SOME COMBINATIONS MIGHT CAUSE TROUBLE 
                
        Hovering: Shows trait[:-10] as default, this can be changed in 
                                "add time data to grid" """
    
    ### PRELIMINARY STEPS
    print("Preparing plot")
    point_size = 6
    # Define colors for the traits (after simulation is finished)
    grid = insert_colours(G, times, grid)
    
    # CREATE FIGURE & INSERT FRAMES
    fig = go.Figure()
    
    for step in times:
        vertices, edges = extract_plot_data_from_grid(grid, step, point_size=point_size)
        fig.add_trace(vertices)
        fig.add_trace(edges)
    
    # DEFINE SLIDER
    fig.data[len(fig.data)-1].visible = True
    fig.data[len(fig.data)-2].visible = True

    num_frames = int (len(fig.data)/2)

    steps = []
    for i in range( num_frames ):
        step = dict(
            method="update",
            args=[  {"visible": [False] * len(fig.data)}   ],
            label = "Size: " + str(max(1,times[i]))
        )
        step["args"][0]["visible"][2*i] = True  # Toggle this trace to "visible"
        step["args"][0]["visible"][2*i+1] = True  # Toggle this trace to "visible"
        steps.append(step)

    sliders = [ dict( active=len(fig.data)-1,  steps=steps) ]
    
    # DEFINE LAYOUT
    axis_range = 1.5* max( v.position[0] for v in G.V[1:]  )
    fig = define_camera_and_layout(fig, axis_range, sliders) 
    
    # CREATE FIGURE AND STORE IT
    EVAL_data_output.create_folder(file_name)
    py.plot(fig, filename= "output/" + file_name + "/spatial_growth.html")
    
    return