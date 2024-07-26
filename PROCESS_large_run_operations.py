"competition_factor"#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:32:35 2019

@author: floriankreten
"""

""" Does a simulation and evaluates/stores at certain step-points """


import PROCESS_process
import PROCESS_mutation
import INPUT_rates
import EVAL_single_graph_evaluation
import numpy.random

import os



def large_run(filename, eval_steps, rad_size = 800):
    """ Does a simulation up to desired size of graph
        Stores and evaluates at each eval_steps[i]
        Standard radius rad_size= 800 is enough for up to 20 Mio vertices """
    if len(eval_steps) == 0:
        raise Exception("Empty eval steps")
    if eval_steps[-1] >= 10000000:
        count = 1000000
        add = "Mio"
    else:
        count = 1000
        add = "k"
    
    step = eval_steps[0]
    G = PROCESS_process.simulate_graph(border = step, spatial_size = rad_size)
    file = filename + "/" + str( int(step/count)) + add
    
    EVAL_single_graph_evaluation.evaluation_routine(G, file)
    
    for step in eval_steps[1:]:
        G = PROCESS_process.continue_simulation(G, border = step)
        file = filename + "/" + str( int(step/count)) + add
        EVAL_single_graph_evaluation.evaluation_routine(G, file)
    
    #do final eval with biopsy
    directory = os.getcwd()
    dirName = str ( directory + "/output/" + file )
    print("Name before final eval", dirName)
    EVAL_single_graph_evaluation.final_routine(G,dirName)
    
    return G

def place_by_hand(filename, shooting_size, end_size=20000000, rad_size = 800, position = "random"):
    """ Simulates and places a strong mutation at a specified size-point
        Evaluates and stores the result (but no raw data/graph.db is stored)
        Input:  shooting_size is the size of the graph when the mutant is placed
                position = "random" (anywhere), "active" (random active vertex),
                                or "start" (root of the tree / first vertex)
    """
    
    if INPUT_rates.parameters["type_of_rate"] != "Place_strong_ones":
        raise Exception("Wrong rate-setup for placing strong ones")
    
    # simulate up to desired size
    G = PROCESS_process.simulate_graph(border = shooting_size, spatial_size = rad_size)
    
    # place a mutation at a certain vertex
    if position == "random":
        # random vertex
        height = G.count
        x = numpy.random.randint(0,height)
        filename += "_rand"
    elif position == "active":
        # random active vertex
        poss = [v for v in G.V[1:] if v.status_branching == True]
        height = len(poss)
        x = numpy.random.randint(0,height)
        x = poss[x].key
        filename += "_act"
    elif position == "start":
        # starting vertex
        filename += "_start"
        x = 1
    else:
        raise Exception("position definition error")    
    PROCESS_mutation.EVENT_mutation(G, x, "ignore")
    print("Set mutation")
        
    # continue up to endsize and evaluate result
    G = PROCESS_process.continue_simulation(G, border = end_size)
    EVAL_single_graph_evaluation.small_mutation_routine(G, filename)
    
    return G
        