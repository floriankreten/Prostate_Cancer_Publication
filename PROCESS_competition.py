#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:01:14 2018

@author: floriankreten
"""

""" All you need for executing and calculating competition events
    Calles from INPUT_rates if something needs to be done """
    
""" the "ignore"-parameter is used when initializing the graph
    in this case, some update-steps are skipped
"""

import INPUT_rates as rates

import PROCESS_mutation as mutation
import PROCESS_growth as growth
#import PROCESS_radial_growth as radial_growth

from PROCESS_event_manager import UPDATE_EVENT_manager_competition


def EVENT_competition_at_v(G,x,p):
    """ A competition Event at x was triggered
        the sampling from the previous step can be re-used
        chooses an outgoing edge x->y and leads to EVENT_competition """
    v = G[x]
    
    if p <= v.rate_competition_in_edge:
        EVENT_competition(G,x,v.in_edge)
    elif p <= v.rate_competition_in_edge +  v.rate_competition_out_edge_1:
        EVENT_competition(G,x,v.out_edge_1)
    else:
        EVENT_competition(G,x,v.out_edge_2)
    return


def EVENT_competition(G,x,y):
    """ competition from x -> y
        genotype at x expands over to y
        then the rates are updated
    """
    
    if G[x].status_competing == False:
        print(G[x])
        print(G[y])
        raise ValueError ("Competition error, point not competing")
    
    # the trait in x takes over y
    G[y].trait = G[x].trait
    G[y].fitness = G[x].fitness
    
    # possible changes induced:
    mutation.UPDATE_rate_mutation(G,y)
    growth.UPDATE_rates_growth_and_branch(G,y)
    #radial_growth.UPDATE_rates_radial_growth(G,y)
    
    # competition possibly changes rates on edges adjacent to y
    # also updates x
    UPDATE_rates_competition_environment(G,y)
    return
    
def UPDATE_rates_competition_environment(G,x,*ignore):
    """ updates the competition rates of x
        and all (maximal 3) adjacent points """
    UPDATE_rates_competition(G,x,*ignore)
    
    v=G[x]
    UPDATE_rates_competition(G,v.in_edge,*ignore)
    
    if v.out_edge_1:
        UPDATE_rates_competition(G,v.out_edge_1,*ignore)
    if v.out_edge_2:
        UPDATE_rates_competition(G,v.out_edge_2,*ignore)
    return


def UPDATE_rates_competition(G,x,*ignore):
    """ updates the competition rates for all adjacent edges
        also updates the graph-rate and the event-manager"""
    v = G[x]
    old_rate =  v.rate_competition
    
    # default: deactivate competition. might get activated below
    v.status_competing = False
    
    # competition on 3 edges
    # incoming
    if v.in_edge == None:
        print(v)
        raise ValueError("In-edge not existent")
    
    else:
        new_comp = rates.competition(G,x,v.in_edge)
        v.rate_competition_in_edge = new_comp
        
        if new_comp != 0:
            v.status_competing = True
        

    
    # first outgoing        
    if v.out_edge_1:
        new_comp = rates.competition(G,x,v.out_edge_1)
        v.rate_competition_out_edge_1 = new_comp
        
        if new_comp != 0:
            v.status_competing = True
        

    
    # second outgoing        
    if v.out_edge_2:
        new_comp = rates.competition(G,x,v.out_edge_2)
        v.rate_competition_out_edge_2 = new_comp
        
        if new_comp != 0:
            v.status_competing = True
        
    # update graph-rate and vertex-rate
    v.rate_competition = v.rate_competition_in_edge \
                        + v.rate_competition_out_edge_1 \
                        + v.rate_competition_out_edge_2
    
                  
    diff_total = v.rate_competition - old_rate
    G.graph_rate += diff_total                

    # needed for skipping some false updates when setting up the graph
    if ignore:
        return
    
    G = UPDATE_EVENT_manager_competition(G,x,diff_total)
        
    return


