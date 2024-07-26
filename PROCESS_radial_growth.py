#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 16:56:38 2020

@author: florian
"""

""" WARNING: NOT USED/TESTED AT THE MOMENT: Everything you need to implement radial growth """

""" the "ignore"-parameter is used when initializing the graph
    in this case, some update-steps are skipped
"""

import INPUT_rates as rates
from PROCESS_event_manager import UPDATE_EVENT_manager_radial_growth
import PROCESS_growth

def UPDATE_rates_radial_growth(G,x,*ignore):
    
    v = G[x]
    old_rate = v.rate_radial_growth
    
    if v.status_radial_growing == True:
        new_rate = rates.radial_cell_rate(G,x)
        v.rate_radial_growth = new_rate
        if new_rate == 0:
            G[x].status_radial_growing = False
            
    elif v.status_radial_growing == False:
        new_rate = 0
        v.rate_radial_growth = new_rate
    
    diff_rate = new_rate - old_rate
    G.graph_rate += diff_rate
    G.graph_rate = max( 0.0, G.graph_rate )    

    # needed for skipping some false updates when setting up the graph 
    if ignore:
        return

    G = UPDATE_EVENT_manager_radial_growth(G, x, diff_rate)
    
    return
    
    
def EVENT_radial_growth(G,x):
    """ Vertex G[x] tries to grow one step in radial direction """
    v = G[x]
    if v.status_radial_growing == False:
        raise TypeError("radial inactive vertex got activated")
        
    new_cells = v.radial_cells + 1    
    new_size = rates.radial_size(new_cells)
    
    # create whitelist
    SEARCH_DEPTH = int( 1.0 + 3.0 * new_size)
    # First look backwards in time
    whitelist = []
    whitelist = PROCESS_growth.get_backward_whitelist(
        G=G, key=x, steps_left = SEARCH_DEPTH,
        already_branched = False, whitelist=whitelist)
           
    # If radial growth event, look also forward in time
    edge1 = v.out_edge_1
    edge2 = v.out_edge_2
    for edge in [edge1, edge2]:
        if edge:
            whitelist = PROCESS_growth.get_forward_whitelist(
                G=G, key=edge, steps_left =SEARCH_DEPTH,
                whitelist=whitelist)
    pos = v.position
    
    if G.is_free(pos, new_size, x, whitelist):
        v.radial_cells = new_cells
        v.radial_size = new_size
    else:
        v.status_radial_growing == False
    
    UPDATE_rates_radial_growth(G,x)    
    
    return
