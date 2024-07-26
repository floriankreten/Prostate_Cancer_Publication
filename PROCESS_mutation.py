#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 14:39:03 2018

@author: floriankreten
"""

""" Everything you need for mutating and keeping track of all mutations
    Called from INPUT_rates when an event need to be executed
"""

""" the "ignore"-parameter is used when initializing the graph
    in this case, some update-steps are skipped
"""

import INPUT_rates
import PROCESS_competition
import PROCESS_growth
#import PROCESS_radial_growth
import numpy

from PROCESS_event_manager import UPDATE_EVENT_manager_mutation

def UPDATE_rate_mutation(G,x,*ignore):
    """ in: Graph G, key x
        updates the mutation-rate of the vertex according to rates-file"""
        
    new_rate = INPUT_rates.mutation(G,x)
    old_rate = G[x].rate_mutation
    
    G[x].rate_mutation = new_rate
    if new_rate != 0:
        G[x].status_mutating = True
    else:
        G[x].status_mutating = False
    
    # update graph-rate
    difference =  new_rate - old_rate
    
    G.graph_rate += difference
    if G.graph_rate < 0:
        if not ignore:
            raise ValueError ("TOTAL RATE < 0, ABORT Cancer")

    # needed for skipping some false updates when setting up the graph
    if ignore:
        return
    
    G = UPDATE_EVENT_manager_mutation(G,x,difference)
    
    return
    
        
        
def EVENT_mutation(G,x,*ignore):
    if G[x].rate_mutation == 0:
        if not ignore:
            print(G[x])
            raise ValueError( "Mutation should not be happening, rate = 0")
    
    # change trait (and fitness eventually) to the new one
    # do a bernoulli experiment to determine if weak or strong mutation
    p = numpy.random.random()
    if p <= INPUT_rates.p_rare_mut:
        # do rare mutation
        G.gen_tree.mutate_rare(G,x)
    else:
        # do "normal" mutation
        G.gen_tree.mutate_normal(G,x)
    
    #update competition and growth regarding the new type
    PROCESS_competition.UPDATE_rates_competition_environment(G,x,*ignore)
    PROCESS_growth.UPDATE_rates_growth_and_branch(G,x,*ignore)
    
    # PROCESS_radial_growth.UPDATE_rates_radial_growth(G,x)

    # mutation rate updated at last, might depend on growth = True/False
    UPDATE_rate_mutation(G,x,*ignore)
    
    return
    
#%%
    
class gen_tree:
    """ stores genotypical tree of the graph (is part of a graph G)
        the mutate-operations do everything you need
        only setting the first initial values has to be done explicitely
        the strong_mutations_log additionally stores all strong mutations
    """
     
    def __init__(self):
        # the dictionary stores all types
        # every type points store the name of the first free type-name
        self.book = dict()

    def __setitem__ ( self, key, newtype):
        self.book[key] = newtype
    
    def __getitem__( self ,  y  ):
        return self.book[y]
        
    def mutate_normal(self,G,x):
        # does a mutation at G[x]
        # at the moment being called, v=G[x] still has the old type
        # gives G[x] a new trait and updates the genealogical tree
        old_trait = G[x].trait
        new_trait = self.book[old_trait]   #mutation
        
        G[x].fitness = INPUT_rates.increase_fitness_normal(G,x)  # updates fitness
        
        G[x].trait = new_trait      # updates trait
        G.fitness_book[ new_trait ] = G[x].fitness
        self.book[ new_trait ] = new_trait + (1,)  #insert new trait
        
        #update the next free mutation-child: "next free name-slot"
        change_key = new_trait[-1] + 1
        self.book[ old_trait ] = new_trait[0:-1] + (change_key,)
        #############################################
        
    def mutate_rare(self,G,x):
        # does a rare mutation at G[x]
        # at the moment being called, v=G[x] still has the old type
        # gives G[x] a new trait and updates the genealogical tree
        old_trait = G[x].trait
        new_trait = self.book[old_trait]   #mutation
        
        G[x].fitness = INPUT_rates.increase_fitness_rare(G,x)  # updates fitness
        G[x].trait = new_trait      # updates trait
        G.fitness_book[ new_trait ] = G[x].fitness
        self.book[ new_trait ] = new_trait + (1,)  #insert new trait
        
        #update the old types next free mutation-child
        change_key = new_trait[-1] + 1
        self.book[ old_trait ] = new_trait[0:-1] + (change_key,)
        
        # add the strong mutation to the history-log
        G.strong_mutations_log[G[x].trait]= ( int (G.time)  )
        #############################################