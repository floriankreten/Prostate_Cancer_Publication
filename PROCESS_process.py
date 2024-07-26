#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 10:05:08 2018

@author: floriankreten
"""

""" Central module for running a  simulation """


""" Gillespie-Algorithm acting on a graph and its vertices
    the execution of the different events is handled in the respective modules
        competition, growth, mutation
    The event-manager of the graph is used for samplint the next event:
        a first random variable determines what type of event occurs,
        then an event is chosen from a weighted list of vertices by rejection-sampling
"""

import numpy.random
import gc

import INPUT_rates
import PROCESS_structures

import PROCESS_growth
import PROCESS_mutation
import PROCESS_competition
import PROCESS_radial_growth

import pandas
import EVAL_gen_tree_graph
import EVAL_over_time_plot_plotly as anim

import PROCESS_growth_mechanism

from PROCESS_cell_type_container import cell_type_container

#%%

def simulate_graph(border, spatial_size,
                   key="size", parameters=INPUT_rates.parameters):
    """ Input:  border = target size of graph
                key = "size" per default (for defining arbitrary stop-criteria)
                spatial_size = maximal size of a 3d-coordinates (needs an estimate)
        Out: Graph G
    """
    
    G = init_simulation_graph(spatial_size, parameters)
    G = continue_simulation(G, border, key)
        
    return G

#%%

def growth_step(G):
    """ Chooses the vertex where the growth step happens,
        directs to growth module, executes an event """
    E = G.EVENT_manager
    X_1 = E.num_growing
    X_2 = E.max_rate_growth
    
    x = numpy.random.randint(0,X_1)
    p = numpy.random.random() * X_2
    
    while E.step_list_growth[x].rate_growth < p:
        x = numpy.random.randint(0,X_1)
        p = numpy.random.uniform(0, X_2)
    
    key = E.step_list_growth[x].key
    
    PROCESS_growth.EVENT_growth(G,key)
    
def radial_growth_step(G):
    """ Chooses the vertex where the radial growth step happens,
        directs to radial growth module, executes an event """
    E = G.EVENT_manager
    X_1 = E.num_radial_growing
    X_2 = E.max_rate_radial_growth
    
    x = numpy.random.randint(0,X_1)
    p = numpy.random.random() * X_2
    
    while E.step_list_radial_growth[x].rate_radial_growth < p:
        x = numpy.random.randint(0,X_1)
        p = numpy.random.uniform(0, X_2)
    
    key = E.step_list_radial_growth[x].key
    PROCESS_radial_growth.EVENT_radial_growth(G,key)
    
def mutation_step(G):
    """ Chooses the vertex where the mutation happens,
        directs to mutation module, executes an event """
    E = G.EVENT_manager
    X_1 = E.num_mutating
    X_2 = E.max_rate_mutation
    if X_1 == 0:
        print("NO MUTATION SHOULD HAPPEN, NUM_MUT == 0")
        print("Rate of Mutation:", E.total_rate_mutation)
        print("Total Rate:", E.graph_rate)
        raise ValueError("mut, look upwards for detailed output")
    
    x = numpy.random.randint(0,X_1)
    p = numpy.random.random() * X_2
    
    while E.step_list_mutation[x].rate_mutation < p:
        x = numpy.random.randint(0,X_1)
        p = numpy.random.uniform(0, X_2)
    
    key = E.step_list_mutation[x].key
    
    PROCESS_mutation.EVENT_mutation(G,key) 
    
    
def branch_step(G):
    """ Chooses the vertex where the branch step happens,
        directs to growth module, executes an event """
    E = G.EVENT_manager
    X_1 = E.num_branching
    X_2 = E.max_rate_branch
    
    x = numpy.random.randint(0,X_1)
    p = numpy.random.random() * X_2
    
    while E.step_list_branch[x].rate_branch < p:
        x = numpy.random.randint(0,X_1)
        p = numpy.random.uniform(0, X_2)
    
    key = E.step_list_branch[x].key
    
    PROCESS_growth.EVENT_branch(G,key)
    
    
def competition_step(G):
    """ Chooses the vertex where the competition step happens,
        directs to competition module, executes an event """
    E = G.EVENT_manager
    X_1 = E.num_competing
    X_2 = E.max_rate_competition
    
    try:
        x = numpy.random.randint(0,X_1)
        p = numpy.random.random() * X_2
    except ValueError:
        print("Num competing", X_1)
        L = [ v for v in G.V if v.status_competing ]
        print(L)
        raise ValueError
    
    while E.step_list_competition[x].rate_competition < p:
        x = numpy.random.randint(0,X_1)
        p = numpy.random.uniform(0, X_2)
    
    key = E.step_list_competition[x].key
    
    PROCESS_competition.EVENT_competition_at_v(G,key,p) 
    
#%%

def random_increment_step(G):
    """ Samples an exponential random time-increment
        Samples what kind of event happens
        Updates the graph G (time+event)
    """
    
    E = G.EVENT_manager
    
    # First: choose time-step
    increm = numpy.random.exponential(1/E.graph_rate)
    G.time += increm
    
    # Second: choose what happens: mutation, competition, growth, branch
    
    p = numpy.random.random() * E.graph_rate
    # look at cumulative weights to decide what to do
    acc = E.total_rate_mutation
    if p <= acc:
        mutation_step(G)
        return
    
    acc += E.total_rate_growth
    if p <= acc:
        growth_step(G)
        return
    
    acc += E.total_rate_branch
    if p <= acc:
        branch_step(G)
        return
    
    acc += E.total_rate_radial_growth
    if p <= acc:
        radial_growth_step(G)
        return
    
    competition_step(G)
    return
     
    
#%%
                        
def init_simulation_graph(size, parameters=INPUT_rates.parameters):
    """ Initializes a graph G for a simulation
        Input:  - parameters and types of rates
                  per default defined in the INPUT_rates module
                  Recommended to change the functions only there!
                - estimate for the maximal coordinates of a vertex (infinity-norm)
    """
        
    # initialize rates
    INPUT_rates._set_rates()
    Length = INPUT_rates.length_of_container_boxes()
    G = PROCESS_structures.graph(Length, size) #delta, size for container-initialisation
    
    #there is a hidden vertex at (0,0,0), so the vertex-array starts at 1
    
    key_1=G.add_vertex( (0,0,0) )
    p=G[key_1]
    p.status_growing = True
    p.status_branching = True
    p.status_radial_growing = True
    p.growth_direction = PROCESS_growth_mechanism.normalize( (1,1,1) )
    p.in_edge=key_1  #easiest way for consistent implementation and no errors
    p.fitness = 1.0
    p.trait = (0,)   #starting vertex has type 0
    p.radial_cells = 1
    p.radial_size = INPUT_rates.radial_size(1)
    
    G.parameters = parameters # just a store-function for later purposes
    
    #initialize growth-rates and mutation-rates
    #update mutation, competition, growth in this order
    #the "ignore"-parameter tells the event manager that we are
    # initializing the graph, so it does not update after each insertion
    for key in ( [key_1 ]):
        G.fitness_book[ G[key].trait ] = G[key].fitness
        PROCESS_mutation.UPDATE_rate_mutation(G,key,"ignore")
        PROCESS_growth.UPDATE_rates_growth_and_branch(G,key,"ignore")
        PROCESS_competition.UPDATE_rates_competition(G,key,"ignore")
        PROCESS_radial_growth.UPDATE_rates_radial_growth(G, key, "ignore")
    
    G.EVENT_manager =  PROCESS_structures.EVENT_manager(G)
    
    #set initial values of the gen-tree
    G.gen_tree[ (0,) ] = (0,1)
    G.time = 0.0
    
    return G


def continue_simulation(G, border, key="size", gen_tree_crop = False, terminal_output = True):
    """
    continues a simulation up to a given border = Graph-size
    """
    
    # small pivot for stopping the simulations
    pivot = INPUT_rates.growth_r * 0.1   
    ################################################
    #Simulation
    
    print("Continuing simulation")
    reducestep = G.count + 500000
            
    if key == "size":
        while G.count < border:
            if G.EVENT_manager.total_rate_growth <= pivot:
                if G.EVENT_manager.total_rate_branch <= pivot:
                    print("Graph stopped growing at time", round(G.time,2) )
                    print("growth:", G.EVENT_manager.total_rate_growth,
                          "branch:", G.EVENT_manager.total_rate_branch,
                          "mutation:", G.EVENT_manager.total_rate_mutation,
                          "competition:", G.EVENT_manager.total_rate_competition)
                    break
            
            else:# random step
                random_increment_step(G)
            
            # reduce size of gen-tree after some steps
            if G.count >= reducestep and gen_tree_crop:
                G = reduce_gen_tree(G)
                reducestep += 500000
            
            
    else:
        raise AssertionError("No valid stopping condition")
        
    if terminal_output:
        diam = max( v.position[0] for v in G.V[1:]  )  
        print("Simulation done.\n Resulting size =", G.count,
                  "| Simulated time=", round(G.time, 2),
                  "\n", " Gen-Tree-size", len(G.gen_tree.book), "\n",
                      " radius =", round(diam, 2) )
        G.parameters["time"] = int(G.time)
        
    return G

#%%

def over_time_simulation(frame_points, spatial_size, key="time"):
    """ Simulation that pastes data into a panda-grid for plots over time
    
            creastes a grid of data points that can be used for
            showing the evolution over time of the process
            see EVAL_over_time_plot, the initial setup is managed there
    """
    G = init_simulation_graph(spatial_size)
    grid = pandas.DataFrame()
    
    # stopping criterion: growth stopped (early fluctiations, very rare)
    pivot = INPUT_rates.parameters["growth_per_fitness"] * 0.5
    j = 0
    
    print("Starting simulation")
    
    if key == "time":
        G.time = 0
        times = frame_points
        timestamp = times[j]
        
        if G.time >= timestamp:
            grid = anim.add_data_to_grid(G, timestamp, grid)
            j +=1
            timestamp = times[j]
        
        while G.time < times[-1]:
            if G.EVENT_manager.total_rate_growth <= pivot:
                print("Warning: Growth stopped after time", G.time)
                break
            
            random_increment_step(G)
            
            #update the grid after a "stepsize" has passed
            while G.time >= timestamp:
                print("time:", timestamp)
                grid = anim.add_data_to_grid(G, timestamp, grid)
                j +=1
                try:
                    timestamp = times[j]
                except IndexError:
                    print("Simulation done, evaluated", j, "framepoints")
                    return G, grid
                
    elif key == "size":
        
        sizes = frame_points
        current_size = sizes[j]
        
        if G.count >= current_size:
            grid = anim.add_data_to_grid(G, current_size, grid)
            j +=1
            current_size = sizes[j]
        
        # run simulation, stop at defined sizes
        while G.count < sizes[-1]:
            if G.EVENT_manager.total_rate_growth <= pivot:
                print("Warning: Growth stopped after time", G.time)
                break
            
            random_increment_step(G)
            
            #update the grid after the next goalpost has been reached
            while G.count >= current_size:
                print("size:", G.count)
                grid = anim.add_data_to_grid(G, current_size, grid)
                j +=1
                try:
                    current_size = sizes[j]
                except IndexError:
                    print("Simulation done, evaluated", j, "framepoints")
                    return G, grid
    
    else:
        print("key not implemented yet")
        return


#%%
# The following is useful for large simulations

def reduce_gen_tree(G):
    """ reduces the amount of entries in dictionaries
        gen_tree and fitness_book
        only entries after most recent common ancestor remain
    """
    try:
        tree = EVAL_gen_tree_graph.create_gen_graph(G, find_size=1, find_type="absolute")
    except AssertionError:
        print("Gen-tree unchanged")
        return G
    
            
    del_list = []
    ANC = tree.ancestor_generation
    #self.generations[INDEX - ANC] = set () OF NODES; NOT NAMES

    for name in G.fitness_book.keys():
        L = len(name)
        if L < ANC -1 :
            del_list.append(name)
        else:
            if name not in tree.V:
                del_list.append(name)

    for name in del_list:
        del(G.gen_tree.book[name])
        del(G.fitness_book[name])
        
    # needs to set up new dictionary
    #   (because the old one is deleted in a lazy way)
    new_tree_book = {}
    for key in G.gen_tree.book.keys():
        new_tree_book[key] = G.gen_tree.book[key]
    a = G.gen_tree.book
    G.gen_tree.book = new_tree_book
    del(a)
    
    new_fitness_book = {}
    for key in G.fitness_book.keys():
        new_fitness_book[key] = G.fitness_book[key]
    a = G.fitness_book
    G.fitness_book = new_fitness_book
    del(a)
    gc.collect()
    print("Gen-tree cut")
    
    return G


def resize_container(G, new_radius, container_length):
    """ Updates the cell-type container if needed (i.e. for large sim) """
    G.container = cell_type_container (container_length, new_radius)
    for v in G.V[1:]:
        G.container.insert(G,v.key)
    return G


#%%
""" some routines for creating graphs for testing etc """

def produce_100k_graph():
    border = 100000
    spatial_size = 200
    key = "size"
    parameters = INPUT_rates.parameters
    G = simulate_graph(border, spatial_size, key, parameters)
    return G

def produce_1k_graph():
    border = 1000
    spatial_size = 80
    key = "size"
    parameters = INPUT_rates.parameters
    G = simulate_graph(border, spatial_size, key, parameters)
    return G

def produce_10k_graph():
    border = 10000
    spatial_size = 100
    key = "size"
    parameters = INPUT_rates.parameters
    G = simulate_graph(border, spatial_size, key, parameters)
    return G