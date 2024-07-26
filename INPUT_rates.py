#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 10:22:54 2018

@author: floriankreten
"""


""" This file contains all behavior-defining functions """

import PROCESS_growth_mechanism
# specify the parameters for the simulation
# this dictionary is saved in the graph and in the anc_tree (output-files)

# Memo about radial_growth: module is deactivated at the moment
# uncomment the corresponding "UPDATE"-lines in mutation,competition
#               and the input/output lines in the "output"-file



##### INPUT
parameters = dict()

# Specify before usage
parameters["type_of_rate"] = "Only_One_Type_Of_Drivers"

parameters["fitness_threshold"] = 10.0

# HERE ONLY ONE TYPE OF MUTATIONS
parameters["mutation_rate"] = 0.0001
parameters["fitness_increment"] = 0.25

# as a subset of all drivers, define the probability and impact of a strong driver
parameters["strong_driver_probability"] = 0.0
parameters["strong_fitness_increment"] = 0.0

parameters["rel_branching_freq"] = 0.2   # as proportion of all growth-events

parameters["num_of_growth_trials"] = 5

parameters["vertex_radius"] = 2.0

# internal parameter, do not change
parameters["container_box_len"] = 5

parameters["cells_per_node"] = 15

parameters["growth_per_fitness"] = 0.25   # should be in (0,1]




def _set_rates():
    """ Defines or updates rates, should be called before a simulation.
    
        We recommend to run a new session when changing parameters,
        though it is possible to call this function for re-defining rates.
    """

    global mu, p_rare_mut, weak_f_factor, str_f_factor, fit_MAX, \
            trials, branch_trials,\
            growth_r, branch_r, rad_r, totalg_r, cb_len, vertex_radius
            
    global length_of_container_boxes, radial_size, radial_cell_rate,\
            cells_of_child, cells_of_branchings
            
            
    global      growth, branch, mutation,\
                increase_fitness_normal, increase_fitness_rare,\
                competition, growth_direction, num_of_senses,\
                num_of_branch_senses, branch_directions
    
    if parameters["type_of_rate"] == "Only_One_Type_Of_Drivers":
                            
                    #MUTATION AND COMPETITION
        ##############################################

        # effective mutation rate
        mu_eff = parameters["mutation_rate"] * parameters["cells_per_node"]
        weak_f_factor = parameters["fitness_increment"] 
        
        # probability of a "rare" strong driver given a mutation event
        p_rare_mut = 0.0
        str_f_factor = 0.0
        
        
        #UNDERLYING MECHANICS OF BRANCHING PROCESS
        ###############################################

        """ Radial growth events are turned off, module not tested """
        rad_r = 0.0
        
        # growth, basic unit per fitness 
        growth_r = parameters["growth_per_fitness"]
        branch_r = parameters["rel_branching_freq"] * growth_r
        growth_r = growth_r - branch_r
        
        vertex_radius = parameters["vertex_radius"]
        cb_len = parameters["container_box_len"]
        
                     
        trials  = parameters["num_of_growth_trials"]
        # How often a branch can sense forward for free space
        # gets inactive if not is successful, else grows
        # has quadratically many tries compared to a normal growth event                
        branch_trials = int ( trials ** 2 )
        fit_MAX = parameters["fitness_threshold"]

        
        ###############################################
        
        
        def length_of_container_boxes():
            """ TODO: Define as function of max_poss_vertex_size /min_size """
            return cb_len
        
        def radial_size(num_cells):
            """ Radial "radius"=size of a vertex with num_cells radial cells """
            return vertex_radius
        
        def radial_cell_rate(G,x):
            """ Speed of radial growth of a vertex """
            return rad_r
            
        def cells_of_child(G,x):
            """ At a growth event: radial cells of child """
            return G[x].radial_cells
        
        def cells_of_branchings(G,x):
            """ At a branching event: radial cells of branches """
            a,b = int(G[x].radial_cells / 2), int(G[x].radial_cells / 2)
            a,b = max(a,1), max(b,1)
            return a,b
        
        def growth(G,x):
            """ in: Graph G, key x
                returns growth-rate of a growing vertex"""
            return G[x].fitness * growth_r
            
        def branch(G,x):
            """ in: Graph G, key x
                returns branch-rate of a growing vertex """
            return G[x].fitness * branch_r
        
        def mutation(G,x):
            """ returns new mutation-rate
                only used when mutating
                when called, G[x] already stores the new type"""
            return mu_eff * G[x].fitness
    
        def increase_fitness_rare(G,x):
            """ Called at the time of driver mutation,
                sets the fitness of the new genotype """
            raise ValueError("Rare Driver Mutation Called")
            return G[x].fitness
            
        def increase_fitness_normal(G,x):
            """ Called at the time of mutation,
                sets the fitness of the new genotype """
            f = G[x].fitness
            f = f * (1.0 + weak_f_factor * ( 1.0 - f/fit_MAX) )
            return f
        
        def competition(G,x,y):
            """ in: Graph G, keys x and y
                returns rate of genotypical change x->y """
        
            v = G[x]
            if v.trait == G[y].trait:
                return 0.0
            elif v.fitness <= G[y].fitness:
                return 0.0
            else:
                k = 1.0 # every vertex has an in-edge (ignoring the root)
                if v.out_edge_1:
                    k+=1.0
                if v.out_edge_2:
                    k+=1.0
                return (v.fitness-G[y].fitness) / k
        
        # possible new growth direction
        def growth_direction(G,x):
            """     input: Graph and key
                    Adds noise to the growth direction, returns vector """
            noise_strength = 60 #maximal strength of noise in degrees <= 90°
            new_direction = PROCESS_growth_mechanism.new_growth_direction (G[x].growth_direction,
                                                        noise_strength)
            return new_direction
        
        # number of forward-sensings given a growth-event
        def num_of_senses(G,x):
            return trials
        
        # number of forward-senses given a branch-event
        def num_of_branch_senses(G,x):
            return branch_trials
    
        # 2 possible branches in different direction
        def branch_directions(G,x):
            """ generates 2 random branch angles in opposing directions """
            alpha = 60 # determine the angle of bifurcation <= 90°
            a,b = PROCESS_growth_mechanism.branch_directions(G[x].growth_direction,alpha)
            return ( a,b )
        
    else:
        raise ValueError("Rates Definition failed")
