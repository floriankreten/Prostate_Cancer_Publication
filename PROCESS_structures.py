#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 09:39:54 2018

@author: floriankreten
"""

""" All Information is stored in one large graph, typically denoted as G
    
    G consists of   a list of vertices V
                    an event manager for sampling of random events
                    a dictionary of all mutations
                    the parameters used for a simulation
                    a container that allows to quickly check neighboring
                        (in space) vertices of a vertex G[x]
                    a "hidden vertex" indexed by 0, that should never be
                        accessed and turned to be useful
                        WARNING: This can cause trouble when this vertex is
                                 called explicitely
"""

import math


import PROCESS_mutation
import PROCESS_radial_growth
import PROCESS_competition
import PROCESS_growth

import INPUT_rates
from PROCESS_event_manager import EVENT_manager
from PROCESS_cell_type_container import cell_type_container


class vertex:
    """ stores all information of vertex:
        key, position, trait, fitness, adjacent edges, growth-direction
        RATES:       growth
                     branch
                     mutation
                     competition in_edge
                     competition out_edge 1
                     competition out_edge 2
                     total
    """
    
    __slots__ = ["position", "trait",
                              "fitness", "key",
                              
                              "in_edge", "out_edge_1",
                              "out_edge_2",
                              
                              "radial_size",
                              "radial_cells",
                              
                              "status_growing",
                              "status_branching",
                              "status_competing",
                              "status_mutating",
                              "status_radial_growing",
                              
                              "index_growth",
                              "index_branch",
                              "index_competition",
                              "index_mutation",
                              "index_radial_growth",
                              
                              "rate_mutation",
                              "rate_growth",
                              "rate_branch",
                              "rate_radial_growth",
                              "growth_direction",
                              "rate_competition",
                              "rate_competition_in_edge",
                              "rate_competition_out_edge_1",
                              "rate_competition_out_edge_2",
                              
                              "cluster_group",
                              "in_biopsy",
                              ]
    
    def __init__ (self):
        
        self.position=False     # coordinates in 3d
        self.trait=False        # trait in N*
        self.fitness=False      # fitness in R
        
        self.key= None          # key in vertex-list
        
        self.in_edge=False      # key of in-edge-vertex (parent)
        self.out_edge_1=False   # key of first out-edge-vertex (child)
        self.out_edge_2=False   # key of second out-edge-vertex (child)
        
        self.radial_cells = False
        self.radial_size = False
        
        self.status_growing= False    # indicator of growth
        self.status_branching = False # indicator of branching
        self.status_competing= False  # indicator of competition
        self.status_mutating = False  # indicator for possible mutations
        self.status_radial_growing = False
        
        # indices for the lists used for random sampling
        self.index_growth = None
        self.index_branch = None
        self.index_competition = None
        self.index_mutation = None
        self.index_radial_growth = None
        
        # rates of events happening at the vertex
        self.rate_mutation=False  #mutation-rate
        
        self.rate_growth=False      # growth rate (depends on type+fitness)
        self.rate_branch=False      # branch rate (depends on type+fitness)
        self.growth_direction=False # normed direction-vector (if growing)
        
        self.rate_competition_in_edge=False
        self.rate_competition_out_edge_1=False
        self.rate_competition_out_edge_2=False
        self.rate_competition = False
        self.rate_radial_growth = False
        
        self.cluster_group = None
        self.in_biopsy = None
        
    
    def reset(self):
            self.trait = (0,)
            self.fitness = 1.0
            
            self.status_growing= False    
            self.status_branching = False 
            self.status_competing= False  
            self.status_mutating = False  
            #self.status_radial_growing = False
            
            self.index_growth = None
            self.index_branch = None
            self.index_competition = None
            self.index_mutation = None
            #self.index_radial_growth = None
            
            self.rate_mutation=False
            self.rate_growth=False      
            self.rate_branch=False      
            self.growth_direction=False 
            
            self.rate_competition_in_edge=False
            self.rate_competition_out_edge_1=False
            self.rate_competition_out_edge_2=False
            self.rate_competition = False
            #self.rate_radial_growth = False
    

    def UPDATE_rates_complete(self, G, *ignore):
        """ Update rates around the vertex """
        if ignore:
            PROCESS_growth.UPDATE_rates_growth_and_branch(G,self.key,"ignore")
            PROCESS_competition.UPDATE_rates_competition(G,self.key,"ignore")
            PROCESS_mutation.UPDATE_rate_mutation(G,self.key,"ignore")
            PROCESS_radial_growth.UPDATE_rates_radial_growth(G,self.key,"ignore")
        else:
            PROCESS_growth.UPDATE_rates_growth_and_branch(G,self.key)
            PROCESS_competition.UPDATE_rates_competition(G,self.key)
            PROCESS_mutation.UPDATE_rate_mutation(G,self.key)
            PROCESS_radial_growth.UPDATE_rates_radial_growth(G,self.key)
            
        
    def __str__(self):
        """ Print function for vertices, rather long and detailed """
        
        output = "\n__" + str(self.key) + "__"
        output += "\n position: "
        output += str(self.position)
        output += "\n trait: "
        output += str(self.trait)
        output += "\n fitness: "
        output += str(self.fitness)
        
        output += "\n radial cells/radius: "
        output += str(self.radial_cells) + " | " + str(self.radial_size)
        output += "\n radial growth | rate: "
        output += str(self.status_radial_growing) + " | " + str(self.rate_radial_growth)
        
        output += "\n in_edge: "
        output += str(self.in_edge)
        output += "\n out_edges: "
        output += str(self.out_edge_1) + " | " + str( self.out_edge_2)
        output += "\n status_growing: "
        output += str(self.status_growing)
        output += "\n status_branching: "
        output += str(self.status_branching)
        output += "\n rate_growth: "
        output += str(self.rate_growth)
        output += "\n rate branch: "
        output += str(self.rate_branch)
        output += "\n growth_direction: "
        output += str(self.growth_direction)
        
        output += "\n status_competition: "
        output += str(self.status_competing)
        
        output += "\n competition-rates: "
        output += str(self.rate_competition_in_edge) + " | " \
                    + str(self.rate_competition_out_edge_1) + " | " \
                    + str(self.rate_competition_out_edge_2)
                    
        output += "\n mutation-rate: "
        output += str (self.rate_mutation)
        
        if self.in_biopsy != None:
            output += "\n in biopsy: "
            output += str(self.in_biopsy)
        
        if self.cluster_group != None:
            output += "\n cluster group: "
            if self.cluster_group != -42:
                output += str(self.cluster_group)
            else:
                output += "not classified"
        
        return output



class graph:
    """ Stores vertices in a list, edges are attributes of vertex-class
        .count gives entry of last vertex = num of vertices (first index = 1)
        needs list-keys of the vertices for all functions
        spatial structure saved in the cell-type-container"""
    
    # for better index-handling, there is a hidden vertex at (0,0,0), see
    #                                       process.init_simulation_graph
    # this is the only exception with key = 0 = False
    # hence indices start at 1 and True/False statements can be used
    
    def __init__ (self, container_length, size):
        INPUT_rates._set_rates()
        self.V = []
        self.count = -1
        self.graph_rate = 0.0
        
        self.EVENT_manager = None
        # needs to be initialized after setting up the initial condition
        
        self.container = cell_type_container (container_length, size)
        self.gen_tree = PROCESS_mutation.gen_tree()
        self.fitness_book = dict () #stores the fitness of every type
        self.add_hidden_vertex()    #hidden vertex, never use this one
        self.time = 0.0
        self.parameters = None  # stores the parameters used for simulation
        self.strong_mutations_log = dict()
    
    def __setitem__ (self, y, value):
        self.V[y] = value
    
    def __getitem__(  self ,  y  ):
        return self.V[y]
    
    def add_hidden_vertex(self):
        v = vertex()
        self.V.append(v)
        self.count += 1
        v.key = self.count
        v.position = ("Dont access hidden vertex with index 0 = False")
        return self.count
    
    def add_vertex(self,x):
        v = vertex()
        self.V.append(v)
        self.count += 1
        v.key = self.count
        v.position = x
        self.container.insert(self,v.key)
        return self.count #this is the key of the new vertex
    
    def add_edge(self,KEY_X,KEY_Y):
        """ adds an edge from x -> y"""
        if KEY_X == KEY_Y:
            raise ValueError ("Edge-Error: Start-vertex = End-vertex!")
        p=self.V[KEY_X]
        q=self.V[KEY_Y]
        if p.out_edge_1 == False:
            p.out_edge_1 = KEY_Y
        elif p.out_edge_2 == False:
            if p.out_edge_1 == KEY_Y:
                raise ValueError ("Doubled Edge")
            p.out_edge_2 = KEY_Y
        else:
            print("x Out_edges_mistake")
            
        if q.in_edge == False:
            q.in_edge = KEY_X
        else:
            print("y in_egde mistake")
            #if-statements for debugging
            
    def __str__(self):
        for p in self.V[1:]: #the hidden vertex is not shown
            print(p)
        print("\nnumber of vertices:", self.count, \
              "| Graph_rate:", self.graph_rate)
        return "\nGRAPH_FINISHED____XXXX____XXXX____"
    
    def is_free(self,x,new_size,father_key, whitelist):
        return self.container.is_free(self,x,new_size,father_key, whitelist)
    
    def appr_radius(self):
        R = max( v.position[0] for v in self.V[1:]  )
        return R
    
    
class mini_biopsy_graph(graph):
    """ A wrapper for easy evaluation of a biopsy
        After declaring it as minigraph, the existing routines can be used """
    def __init__(self, vertex_set, parent_graph):
        INPUT_rates._set_rates()
        self.V = []
        self.count = -1
        self.graph_rate = 0.0
        
        self.add_hidden_vertex()    #hidden vertex, never use this one
        self.time = 0.0
        self.parameters = None  # stores the parameters used for simulation
        self.V += list(vertex_set)
        self.count = len(vertex_set)
        self.fitness_book = parent_graph.fitness_book
        self.strong_mutations_log = parent_graph.strong_mutations_log
        




############################################################################
########## Operations which are useful but not defined for normal tuples   
def PLUS(a,b):
    """Difference of two arrays"""
    k = len(a)
    if k != len(b):
        raise ValueError("Länge", a, "!= Länge", b)
    else:
        c=[ 0 for i in range(k) ]
        for i in range(k):
            c[i] = a[i] + b[i]
        return tuple(c)

 
def MINUS(a,b):
    """Difference of two arrays"""
    k = len(a)
    if k != len(b):
        raise ValueError("Länge", a, "!= Länge", b)
    else:
        c=[ 0 for i in range(k) ]
        for i in range(k):
            c[i] = a[i]-b[i]
        return tuple(c)
        
def ABS (a,*b):
    """Absolute distance of two arrays"""
    if b:
        c = MINUS (a,b[0])
    else:
        c = a
    ABS = sum ( s**2 for s in c)
    return math.sqrt(ABS)

def ABS2 (a,*b):
    """Absolute distance^2 of two arrays"""
    if b:
        c = MINUS (a,b[0])
    else:
        c = a
    ABS = sum ( s**2 for s in c)
    return ABS