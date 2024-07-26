#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 10:05:09 2018

@author: floriankreten
"""

""" All you need for calculating and executing growth-events
    Geometrical calculations are shifted to PROCESS_growth_mechanism
"""

""" the "ignore"-parameter is used when initializing the graph
    in this case, some update-steps are skipped
"""

import INPUT_rates as rates
import PROCESS_structures as structures
import PROCESS_mutation as mutation
import PROCESS_radial_growth as radial_growth

from PROCESS_event_manager import UPDATE_EVENT_manager_growth

def UPDATE_rates_growth_and_branch(G,x,*ignore):
    """ IN: Graph G, key x
        updates and returns growth-rate
        needs to have correct "state"=active/inactive as implicit part of the input"""
    
    #updates growth rate
    old_growth=G[x].rate_growth
    if G[x].status_growing == True :
        new_growth=rates.growth(G,x)
        G[x].rate_growth=new_growth
        if new_growth == 0:
            G[x].status_growing = False
        
    elif G[x].status_growing == False:
        #updates growth rate to 0
        new_growth= False
        G[x].rate_growth=new_growth
        

    # updates branch-rate
    old_branch=G[x].rate_branch
    if G[x].status_branching == True:
        new_branch = rates.branch(G,x)
        G[x].rate_branch=new_branch
        if new_branch == 0:
            G[x].status_branching = False
        
    elif G[x].status_branching == False:
        #updates branch rate to 0
        new_branch= False
        G[x].rate_branch = new_branch
        
    if G[x].status_branching == G[x].status_growing == False:
        G[x].growth_direction = False
        
    #updates total rate of vertex and graph
    diff_growth = new_growth - old_growth
    diff_branch = new_branch - old_branch
    
    diff_total  = diff_growth + diff_branch
    G.graph_rate += diff_total
    
    G.graph_rate = max( 0.0, G.graph_rate )
    
    # needed for skipping some false updates when setting up the graph
    if ignore:
        return
    
    G = UPDATE_EVENT_manager_growth(G, x, diff_growth, diff_branch)
    
    return
        
def EVENT_growth(G,x):
    """ IN: Graph G, key x
        OUT: does a growth step in this point
             if space available, connects to new point w """
             
    if G[x].status_growing == False:
        raise ValueError ("Growth-error: Inactive point growing")
        
    new_cells = rates.cells_of_child(G,x)
    new_size = rates.radial_size(new_cells)
    
    SEARCH_DEPTH = int( 1.0 + 3.0 * new_size)
    whitelist = []
    whitelist = get_backward_whitelist(
                G=G, key=x, steps_left = SEARCH_DEPTH,
                already_branched = False, whitelist=whitelist)
    
    #sense forward for new space    
    for i in range(rates.num_of_senses(G,x)):
        new_direction = rates.growth_direction(G,x)
        new_coordinate = structures.PLUS ( G[x].position , new_direction )
        growth_activity = G.is_free(new_coordinate, new_size, x, whitelist)
        #CHECK_site_is_free( G, new_coordinate, new_size, x )
        
        if growth_activity == True:
            break
    
    # creates a new vertex if space is available    
    if growth_activity == True:
        
        #specify new vertex, copy values from predecessor
        y = G.add_vertex(new_coordinate)
        w=G[y]
        G.add_edge(x,y)
        w.growth_direction = new_direction
        w.trait = G[x].trait
        w.fitness = G[x].fitness
        
        w.radial_cells = new_cells
        w.radial_size = new_size
            
        #update activities of new vertex
        w.status_growing = True
        w.status_branching = True
        w.status_mutating = True
        w.status_competing = False
        w.status_radial_growing = True
        mutation.UPDATE_rate_mutation(G,y)
        radial_growth.UPDATE_rates_radial_growth(G,y)
        UPDATE_rates_growth_and_branch(G,y)
        # no competition update (same genotype on x and y)
        
        #deactivate old vertex after growth
        G[x].status_growing = False
        G[x].status_branching = False
        UPDATE_rates_growth_and_branch(G,x)
        mutation.UPDATE_rate_mutation(G,x)
        
    elif growth_activity == False:
        #set growth-activity to 0
        G[x].status_growing = False
        UPDATE_rates_growth_and_branch(G,x)
        mutation.UPDATE_rate_mutation(G,x)
    
    return
        


def EVENT_branch(G,x):
    """ IN: Graph G, key x
        OUT: does a branch step in this point
            if space available, connects to new point w and v
            more debug options can be done via the container-file"""
        
    if G[x].status_branching == False:
        raise ValueError ("Growth-error: Inactive point branching")
    
    Flag = False
    new_cells1, new_cells2 = rates.cells_of_branchings(G,x)
    new_size1, new_size2 = rates.radial_size(new_cells1), rates.radial_size(new_cells2)
    
    SEARCH_DEPTH = int( 1.0 + 3.0 * max(new_size1, new_size2) )
    whitelist = []
    whitelist = get_backward_whitelist(
                G=G, key=x, steps_left = SEARCH_DEPTH,
                already_branched = False, whitelist=whitelist)
    
    for i in range(rates.num_of_branch_senses(G,x)):
        # tries to find good branch-angle for 2 growing tips
        new_directions = rates.branch_directions(G,x)
        
        first_direction  = new_directions[0]
        first_coordinate = structures.PLUS ( G[x].position, first_direction )
        growth_activity_1 = G.is_free(first_coordinate, new_size1, x, whitelist)
        
        second_direction  = new_directions[1]
        second_coordinate = structures.PLUS ( G[x].position, second_direction )
        growth_activity_2 = G.is_free(second_coordinate, new_size2, x, whitelist)
        
        if growth_activity_1 == True and growth_activity_2 == True:
            Flag = True
            break
    
    # not enough space to branch, try growth if possible
    if Flag == False:
        G[x].status_branching = False
        UPDATE_rates_growth_and_branch(G,x)
        mutation.UPDATE_rate_mutation(G,x)
        if G[x].status_growing:
            EVENT_growth(G,x)
        return
    
    else:
        
    #branches to the first point
        y = G.add_vertex(first_coordinate)
        w=G[y]
        G.add_edge(x,y)
        w.growth_direction = first_direction
        w.trait = G[x].trait
        w.fitness = G[x].fitness
        
        w.radial_cells = new_cells1
        w.radial_size = new_size1
       
        #update activities of new vertex
        w.status_growing = True
        w.status_branching = True
        w.status_mutating = True
        w.status_competing = False
        w.status_radial_growing = True
        mutation.UPDATE_rate_mutation(G,y)
        radial_growth.UPDATE_rates_radial_growth(G,y)
        UPDATE_rates_growth_and_branch(G,y)
    
    # branches to the second point
        z = G.add_vertex(second_coordinate)
        v=G[z]
        G.add_edge(x,z)
        v.growth_direction = second_direction
        v.trait = G[x].trait
        v.fitness = G[x].fitness
        
        v.radial_cells = new_cells2
        v.radial_size = new_size2
        
        #update activities of new vertex
        v.status_growing = True
        v.status_branching = True
        v.status_mutating = True
        v.status_competing = False
        v.status_radial_growing = True
        mutation.UPDATE_rate_mutation(G,z)
        radial_growth.UPDATE_rates_radial_growth(G,z)
        UPDATE_rates_growth_and_branch(G,z)

    #deactivate growth/branch of old vertex
        G[x].status_growing = False
        G[x].status_branching = False
        UPDATE_rates_growth_and_branch(G,x)
        mutation.UPDATE_rate_mutation(G,x)

        #there is no competition from x to y or z
        #competition.UPDATE_rates_competition(G,x)
        return
    
#%%
def get_backward_whitelist(G, key, steps_left, already_branched, whitelist):
    """ Search the nearby graph for whitelisted vertices
        those will not "block" a growth/branch-event """
    if steps_left == 0:
        return whitelist
        
    whitelist.append(key)
    v = G.V[key]
        
        # if we reach the first branching point and have enough steps left,
        # some of the branching-structure is added to the whitelist
        
    edge1 = v.out_edge_1
    edge2 = v.out_edge_2
        
    if not already_branched:
        if edge1 and edge2:
            already_branched = True
            for out_edge in [edge1, edge2]:
                if out_edge not in whitelist:
                    whitelist = get_forward_whitelist (
                        G, out_edge, steps_left -1, whitelist)
                        
    if v.in_edge:                
        whitelist = get_backward_whitelist(
            G, v.in_edge, steps_left -1, already_branched, whitelist)
                
    return whitelist
        
def get_forward_whitelist(G, key, steps_left, whitelist):
    if steps_left == 0:
        return whitelist
        
    v = G.V[key]
    whitelist.append(key)
    if v.out_edge_1:
        whitelist = get_forward_whitelist(
            G, v.out_edge_1, steps_left -1, whitelist)
    if v.out_edge_2:
        whitelist = get_forward_whitelist(
            G, v.out_edge_2, steps_left -1, whitelist)            
        
    return whitelist