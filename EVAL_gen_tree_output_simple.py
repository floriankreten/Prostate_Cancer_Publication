#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 17:25:37 2023

@author: floriankreten
"""

from pathlib import Path

import numpy as np
import plotly.offline as py
import plotly.graph_objs as go

import pandas as pd

import EVAL_data_output
import EVAL_gen_tree_graph


#%%

def print_tree_simple(tree, filename, fullpathname = False,
        minimal_size = 0.01, size_mode= "relative", 
        biopsy_extra = False, use_colors = False):
    """ Simple plot of a gen-tree. Shows all mutations above a threshold """
    
    
    if size_mode == "relative":
        min_abs_size = int(float(tree.total_weight) * minimal_size)
    elif size_mode == "absolute":
        min_abs_size =  minimal_size
    
    def weight_condition(v):
        if v.branch_weight >= min_abs_size:
            return True
        else:
            return False
    
    # 1) set coordinates of nodes
    subtree_size = {}
    
    # size of subtree that is adjacent to a node, needed for placement in 2D plane
    for i in range ( tree.height -1 , -1 , -1) :
        gen = tree.generations[i]
        for v in gen:
            subtree_size[v] = max (
                    sum ( subtree_size[child] for child in v.children
                         if weight_condition(child) ),
                    1 )
            
            
    #assign (x,y) values to nodes
    coordinates = {}
    root_position = 0   # common ancestor is placed in the center
    root = tree.root


    #Finds position of the nodes, recursive node-placement
    #Defines coordinates-dictionary
    scale = 0.8
            
    def place_subtree(self, v, position, factor):
        """recursively sets position of nodes"""
            
        coordinates [v] = [ position, - len(v.name) ]
            
        if len(v.children) > 1:
            newfactor = factor * scale
        else:
            newfactor = factor
        left_border = position - newfactor * 0.5 * subtree_size[v]
        position = left_border
            
        for child in v.children:
            if weight_condition(child):
                position += (subtree_size[child] * 0.5 * newfactor)
                place_subtree(self, child, position, factor)
                position += (subtree_size[child] * 0.5 * newfactor)
            
    place_subtree(tree, root, root_position, 1)
        

    # 2) define edges of the tree
    x_edges = []
    y_edges = []
        
        
    for v in tree.V.values():
        
        if weight_condition(v):
            
            start = coordinates[v]
            for child in v.children:
                
                if weight_condition(child):
                    
                    end = coordinates[child]
            
                    x_edges.append(start[0])
                    x_edges.append(end[0])
                            
                    y_edges.append(start[1])
                    y_edges.append(end[1])
                            
                    y_edges.append(None)
                    x_edges.append(None)

    
    edges = go.Scatter(     x=x_edges,
                            y=y_edges,
                            mode='lines',
                            name='coalescence',
                            line=dict( color='black', width = 3),
                            hoverinfo='none'
                            )
        
    # 3) traits have green dots: one node with info for each mutation
    def text_node(v):
        rel_size = round ( 100 * v.branch_weight /  tree.total_weight , 1)
        liv_size = round ( 100 * v.weight /  tree.total_weight , 1)
        ret = str(v.name ) + " | branch " + str(rel_size) + "%"  + " | living  " + str(liv_size) + "%" 
        return ret
    
    
    if use_colors:
        x_vertices = [ coordinates[v][0]  for v in tree.group_points  ]
        y_vertices = [ coordinates[v][1]  for v in tree.group_points  ]
        trait_markers = [ text_node(v) for v in tree.group_points  ]
    else:
        x_vertices = [ coordinates[v][0]  for v in tree.V.values()  if weight_condition(v)  ]
        y_vertices = [ coordinates[v][1]  for v in tree.V.values()  if weight_condition(v)  ]
        trait_markers = [ text_node(v) for v in tree.V.values()  if weight_condition(v) ]
    

    max_size = 200
    min_size = 30
    
    def size_of_point(v): # linear or sqrt size-scaling, try it out
        # return max(1, max_size * v.weight / tree.total_weight )
        return max(min_size, max_size * np.sqrt(v.weight / tree.total_weight) )
        
    if use_colors:
        colours = list ( v.colour for v in tree.group_points  )
        print("For tree-plot: using pre-defined colors", colours)
        regular_sizes = [ size_of_point(v) for v in tree.group_points  for v in tree.group_points  ]
    else:
        colours = "green"
        regular_sizes = [ size_of_point(v) for v in tree.V.values()  if weight_condition(v) ]
    
    
    types = go.Scatter (
            x = x_vertices,
            y = y_vertices,
            hoverinfo = 'text',
            text = trait_markers,
            opacity = 0.99,
            mode = 'markers',
            name='Genealogical tree of mutations',
            marker = dict(size = regular_sizes, color = colours, colorscale = 'Jet', cmin = 0,
            cmax = 1 )
    )
    
    # add information about biopsies to the tree
    found_bound = 0.02
    
    if biopsy_extra:
        def biopsy_text(v):
            w = [ round( v.is_found[i] / tree.root.is_found[i] * 100, 1) for i in range(len(v.is_found)) ]
            return str(w) + " % in biopsies"

        def biopsy_color(v):
            w = [ v.is_found[i] / tree.root.is_found[i] for i in range(len(v.is_found)) ]
            for i in range(len(v.is_found)):
                if w[i] >= found_bound:
                    return "black"
            else:
                return "red"
                
        biopsy_infos = [ biopsy_text(v) for v in tree.V.values()  if weight_condition(v)  ]
        biopsy_colors =  [ biopsy_color(v) for v in tree.V.values()  if weight_condition(v)  ]
        
        biopsy_nodes = go.Scatter (
                x = x_vertices,
                y = y_vertices,
                hoverinfo = 'text',
                text = biopsy_infos,
                opacity = 1,
                mode = 'markers',
                name='Percentage within biopsies',
                marker = dict(size = min_size*0.7, color = biopsy_colors )
        )

        
    # put things together to a plotly figure and print it
    figure = {'data': [], 'layout': {} }
    noaxis = dict (autorange=True,
                   showgrid=False,
                   zeroline=False,
                   showline=False,
                   ticks='',
                   showticklabels=False )
    
    title = "Ancestral graph V3"

    layout = go.Layout( xaxis = noaxis, yaxis = noaxis, title=title,
                        hovermode = 'closest')
        
    if biopsy_extra:
        figure['data']=[edges, types, biopsy_nodes]
    else:
        figure['data']=[edges, types]
            
    figure['layout']= layout

    # Creates directory for storage
    if not fullpathname:
        dirName = EVAL_data_output.create_folder(filename)
    else:
        dirName = filename
        Path(dirName).mkdir(parents=True, exist_ok=True)
            
    name = str ( dirName + "/" + "ancestral_tree_informative.html" )

    py.plot(figure, filename = name, auto_open=False ) 
    return

#%%

def all_mut_excel_to_small_dict(filename, min_CCF = 0.1):
    df = pd.read_excel(filename)
    mutations_found = {}
    
    for row in df.iterrows():
        values = row[1]
        
        # filter for empty rows
        CCF = values[1]
        if not np.isnan(CCF):
            
            if CCF < min_CCF:
                break
            
            name = tuple(eval(values[8]))
            mutations_found[name] = round(CCF,1)
        
    return mutations_found

def biopsy_value(mutation_name, biopsy_dict):
    K = 0.0
    L = len(mutation_name)
    for key in biopsy_dict.keys():
        if L <= len(key):
            comp = tuple ( key[:L] )
            if comp == mutation_name:
                K = max( K, biopsy_dict[key]  )
    return K


#%%

def modify_root_of_tree(tree):
    # roots the tree at (0,) ("adds" vertices above the root)
    while tree.root.name != (0,):
        new_root_name = tuple(tree.root.name[:-1])
        r = EVAL_gen_tree_graph.node(new_root_name)
        r.children = [tree.root]
        r.weight = 0
        r.branch_weight = tree.root.branch_weight
        r.parent = r
        
        tree.root.parent = r
        new_gen = set()
        new_gen.add(r)
        tree.V[new_root_name] = r
        
        tree.generations.insert( 0, new_gen )
        tree.root = r
    
    return tree

#%%

def print_tree_with_four_biopsies(data_path,
        minimal_size = 0.03,
        size_mode= "relative",
        num_biopsies = 4,
        biopsy_mode = "4_quad_on_different_levels"):

    """ Simple plot of a gen-tree. Shows all mutations above a threshold,
        and if they are found in the biopsies
        Biopsies are represented via 4 dots """
    
    # at the moment only defined for 4 biopsies
    assert(num_biopsies == 4)
    
    # load data from directory (depends on the type of biopsy)
    if not data_path.endswith("/"):
        data_path = data_path + "/"
    
    tree = EVAL_data_output.get_tree_sql(data_path, fullpathname=True)
    tree = modify_root_of_tree(tree)
    
    biopsy_data = []
    
    for i in range(1, num_biopsies+1):
        filename = data_path + biopsy_mode + "/biopsy_" + str(i) + "/all_mutations.xlsx"
        biopsy = all_mut_excel_to_small_dict(filename) 
        biopsy_data.append ( biopsy )
        
    print("All data loaded from " + data_path + biopsy_mode)
    
    # define some helper variables
    if size_mode == "relative":
        min_abs_size = int(float(tree.total_weight) * minimal_size)
    elif size_mode == "absolute":
        min_abs_size =  minimal_size
    
    def weight_condition(v):
        if v.branch_weight >= min_abs_size:
            return True
        else:
            return False
    
    
    # 1) set coordinates of nodes
    subtree_size = {}
    
    # size of subtree that is adjacent to a node, needed for placement in 2D plane
    for i in range ( tree.height -1 , -1 , -1) :
        gen = tree.generations[i]
        for v in gen:
            subtree_size[v] = max (
                    sum ( subtree_size[child] for child in v.children
                         if weight_condition(child) ),
                    1 )
            
            
    #assign (x,y) values to nodes
    coordinates = {}
    root_position = 0   # common ancestor is placed in the center
    root = tree.root


    #Finds position of the nodes, recursive node-placement
    #Defines coordinates-dictionary
    scale = 0.8
            
    def place_subtree(self, v, position, factor):
        """recursively sets position of nodes"""
            
        coordinates [v] = [ position, - len(v.name) ]
            
        if len(v.children) > 1:
            newfactor = factor * scale
        else:
            newfactor = factor
        left_border = position - newfactor * 0.5 * subtree_size[v]
        position = left_border
            
        for child in v.children:
            if weight_condition(child):
                position += (subtree_size[child] * 0.5 * newfactor)
                place_subtree(self, child, position, factor)
                position += (subtree_size[child] * 0.5 * newfactor)
            
    place_subtree(tree, root, root_position, 1)
        

    # 2) define edges of the tree
    x_edges = []
    y_edges = []
        
        
    for v in tree.V.values():
        
        if weight_condition(v):
            
            start = coordinates[v]
            for child in v.children:
                
                if weight_condition(child):
                    
                    end = coordinates[child]
            
                    x_edges.append(start[0])
                    x_edges.append(end[0])
                            
                    y_edges.append(start[1])
                    y_edges.append(end[1])
                            
                    y_edges.append(None)
                    x_edges.append(None)

    
    edges = go.Scatter(     x=x_edges,
                            y=y_edges,
                            mode='lines',
                            name='coalescence',
                            line=dict( color='black', width = 3),
                            hoverinfo='none'
                            )
    

    # 3) get information about mutations in biopsies
    def biopsy_perc(v):
        biopsy_info = []
        for i in range(num_biopsies):
            b_perc = biopsy_value(v.name, biopsy_data[i])
            biopsy_info.append(b_perc)
        return biopsy_info
    
    biopsy_perc_values = [ biopsy_perc(v) for v in tree.V.values()  if weight_condition(v) ]
        
    # 4) traits have green dots: one node with info for each mutation
    
    def text_node(v):
        rel_size = round ( 100 * v.branch_weight /  tree.total_weight , 1)
        node_info = str(v.name ) + " | " + str(rel_size) + "% <br>"
        return node_info
    
    x_vertices = [ coordinates[v][0]  for v in tree.V.values()  if weight_condition(v)  ]
    y_vertices = [ coordinates[v][1]  for v in tree.V.values()  if weight_condition(v)  ]
    
    # info text for hovering_
    trait_markers = [ text_node(v) for v in tree.V.values()  if weight_condition(v) ]
    for i in range(len(trait_markers)):
        biopsy_info = biopsy_perc_values[i]
        biopsy_info = str(biopsy_info) + " % in biopsies"
        trait_markers[i] = trait_markers[i] + biopsy_info
    
    min_size = 30
    max_size = 100
    
    def size_of_point(v):
        return min_size + ( max_size - min_size) * np.sqrt(v.branch_weight / tree.total_weight)
        
    regular_sizes = [ size_of_point(v) for v in tree.V.values()  if weight_condition(v) ]
    
    types = go.Scatter (
            x = x_vertices,
            y = y_vertices,
            hoverinfo = 'text',
            text = trait_markers,
            opacity = 0.99,
            mode = 'markers',
            name='Genealogical tree of mutations',
            marker = dict(size = regular_sizes, color = 'green' )
    )
    
    # parameters for biopsy dot placement
    biopsy_dot_size = 20
    biopsy_dot_spread = 0.02
    
    H = max ( x_vertices  ) - min (x_vertices)
    V = max ( y_vertices  ) - min (y_vertices)
    
    biopsy_dot_spread_x = H * biopsy_dot_spread
    biopsy_dot_spread_y = V * biopsy_dot_spread
    
    def biopsy_shift(x,y,biopsy_ticker):
        if biopsy_ticker == 0:
            xn = x - biopsy_dot_spread_x
            yn = y + biopsy_dot_spread_y
        elif biopsy_ticker == 1:
            xn = x + biopsy_dot_spread_x
            yn = y + biopsy_dot_spread_y
        elif biopsy_ticker == 2:
            xn = x + biopsy_dot_spread_x
            yn = y - biopsy_dot_spread_y
        elif biopsy_ticker == 3:
            xn = x - biopsy_dot_spread_x
            yn = y - biopsy_dot_spread_y
        return xn, yn
    
    biopsy_colors = ["red", "black", "blue", "gold"]
    
    # add information about biopsies
    
    found_in_biopsy_bound = 5 # detected in a biopsy if above threshold
    
    biopsy_dots_x = [ [] for i in range(num_biopsies)]
    biopsy_dots_y = [ [] for i in range(num_biopsies)]
    biopsy_hover = [ [] for i in range(num_biopsies)]
    
    # place dots for biopsies
    i = 0
    for v in tree.V.values():
        if weight_condition(v):
            x = x_vertices[i]
            y = y_vertices[i]
            for j in range(num_biopsies):
                b = biopsy_perc_values[i][j]
                if b >= found_in_biopsy_bound:
                    xn, yn = biopsy_shift(x,y,j)
                    biopsy_dots_x[j].append(xn)
                    biopsy_dots_y[j].append(yn)
                    info = str(b) + "%"
                    biopsy_hover[j].append(info)
            i += 1
            
    biopsy_nodes = []
    
    for j in range(num_biopsies):
        B = go.Scatter (
            x = biopsy_dots_x[j],
            y = biopsy_dots_y[j],
            opacity = 1,
            hoverinfo = 'text',
            text = biopsy_hover[j],
            mode = 'markers',
            name='Biopsy ' + str(j+1),
            marker = dict(size = biopsy_dot_size, color = biopsy_colors[j] )
        )
        biopsy_nodes.append(B)

        
    # put things together to a plotly figure and print it
    figure = {'data': [], 'layout': {} }
    noaxis = dict (autorange=True,
                   showgrid=False,
                   zeroline=False,
                   showline=False,
                   ticks='',
                   showticklabels=False )
    
    title = "Gen graph with biopsy info"

    layout = go.Layout( xaxis = noaxis, yaxis = noaxis, title=title,
                        hovermode = 'closest')
        
    figure['data']=[edges, types]
    
    for i in range(num_biopsies):
        figure['data'].append(biopsy_nodes[i])
            
    figure['layout']= layout

    # Creates directory for storage

    dirName = data_path + biopsy_mode
    Path(dirName).mkdir(parents=True, exist_ok=True)
    
    name = str ( dirName + "/" + "ancestral_tree_informative.html" )
    
    py.plot(figure, filename = name, auto_open=False ) 
    return