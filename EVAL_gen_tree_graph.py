#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:11:15 2019

@author: floriankreten
"""

""" Preferably use this module via EVAL_single_graph_evaluation """

""" a class to work with genealogical trees """

import EVAL_gen_tree_prep as eval_gr
from collections import defaultdict

import plotly.offline as py
import plotly.graph_objs as go

import random
import EVAL_data_output
from ordered_set import OrderedSet

import copy
from pathlib import Path


br_cutoff = 5



#%%

class node:
    
    def __init__ (self, name):
        self.name = name                # a tuple of integers
        self.fitness = None             # the fitness
        self.mean_fitness = None        # the mean_fitness of the subtree
        self.max_fitness = None         # the max_fitness of the subtree
        self.parent = None              # the parent_node (except root)
        self.children = None            # the child-nodes
        self.weight = None              # the amount of living cells
        self.branch_weight = None       # the amount of cells in the subtree
        self.generation = None          # the generation (as num of mutations)
        self.number = None              # each node has a number, internal only
        self.fittest_type = None        # one fittest subtype
        self.group = None               # the vertex assigned to the clustering
        self.colour = None              # a fixed colour for multiple plots
        self.spatial_center = None      # the coordinate of the center (only relevant for group nodes)
        self.is_found = None            # numbers of vertices with that type
                                        # ordered via biopsies

    def __str__(self):
        """ Print_function for nodes """
        
        output = "__" + str(self.name) + "__"
        output += "\n children: "
        output += str( [ child.name for child in self.children])
        output += "\n weight: "
        output +=  str(self.weight)
        output += "\n branch-weight: "
        output += str(self.branch_weight)

        
        output += "\n number_key: "
        output += str( self.number )
        
        if self.group:
            output += "\n group: "
            output += str(self.group.name)

        output += "\n fitness: "
        output += str(self.fitness)
        
        output += "\n fittest type: "
        output += str(self.fittest_type)
        
        output += "\n max_fitness: "
        output += str(self.max_fitness)
        
        output += "\n mean_fitness: "
        output += str(self.mean_fitness)
        
        output  += "\n is found: "
        output += str(self.is_found)
        
        if self.spatial_center:
            output += "\n spatial center: "
            output += str(self.spatial_center)
        
        return output
        
class ancestral_tree:
    # should be built using create_gen_graph(G)
    # Stores only the information concerning all living types
    #       -> types before the common ancestor are ignored
    
    def __init__(self):
        
        self.height = None      # height (in generations)
        self.generations = None # list of generations
        self.V = None           # dictionary of all nodes
        self.ancestor_generation = None     # generations of the ancestor
        self.type_of_number = None # internal only
        self.type_amount = None     # number of types
        self.distances = None       # not done yet
        self.total_weight = None    # size of tumour
        self.root = None            # root-node
        self.group_points = None    # nodes from clustering
        self.unclassified_maximal_fitness = None
        self.parameters = None      # for evaluation: stores the parameters
                                    #   which were used for the simulation
                                    
        self.strong_mutations_log = dict()
        
        return
    
    def __setitem__ (self, y, value):
        self.V[y] = value
    
    def __getitem__(  self ,  y  ):
        return self.V[y]
    
    def pack_graph_infos(self, G, tumour, ANC, in_generations, \
                 branch_weight, centroids, mean_fitness, max_fitness, fittest_type):
        
        self.height = len(in_generations) - ANC
        self.generations = []
        for i in range (self.height):
            self.generations.append( None )
        self.V = {}
        self.ancestor_generation = ANC
        self.type_of_number = []
        self.type_amount = None
        self.distances = None
        self.total_weight = None
        self.root = None
        self.group_points = None
        self.unclassified_maximal_fitness = None
        self.parameters = G.parameters
        
        # process data into usable graph-structure
        for INDEX in range (  ANC, len(in_generations) ) :
            
            in_gen = in_generations[ INDEX ]
            self.generations[INDEX - ANC] = set ()
            tree_gen = self.generations[INDEX - ANC]

            
            for TYP in in_gen:
                # add vertex to graph and generations, add information
                self.V[TYP] = node( TYP )
                v = self.V[TYP]
                v.name = TYP
                tree_gen.add(v)
                
                v.fitness = G.fitness_book[TYP]
                v.mean_fitness = mean_fitness[TYP]
                v.max_fitness = max_fitness[TYP]
                v.fittest_type = fittest_type[TYP]
                v.spatial_center = centroids[TYP]
                v.parent = tuple(TYP[:-1])
                v.children = in_gen[TYP]
                if TYP in tumour:
                    v.weight = tumour[TYP]
                else:
                    v.weight = False
                v.branch_weight = branch_weight[TYP]
                v.generation = INDEX - ANC
                
                if INDEX != len(v.name):
                    raise ValueError ("Input_Error: Generations were unordered")
        
        
        #get rid of eventually unnecessary parts of the gen-tree (which was cut)
        height_reduction = 0
        for k in range (self.height):
            
            if len(self.generations) >= 2:
            
                if len( self.generations[0] ) == 0:
                    self.generations.pop(0)
                
                elif len( self.generations[1] ) == 1 and len (self.generations[0]) == 1:
                    if len(self.generations[0] ) != 1:
                        raise ValueError("Creation of anc-graph failed, \n \
                                         something wrong around the ancestor")
                    
                    # when v is not the MRCA and has no living cells:
                    v = list ( self.generations[0] ) [0]
                    if v.weight == 0:
                        v = v.name
                        del ( self.V[v] )
                        self.generations.pop(0)
                        self.ancestor_generation += 1
                        height_reduction += 1
                 
                elif len( self.generations [0] ) > 1:
                    raise ValueError ( "Cutting tree failed at the root" )
            try:
                if len( self.generations[-1] ) == 0:
                    del(self.generations[-1])
                    height_reduction += 1
            except IndexError:
                print("WARNING: Gen-tree contains no node, \
                      pooling resulted in only 1 living type")
        
        # reset some parameters after cutting and do justice to the root
        self.height -= height_reduction
        self.root = list ( self.generations[0] ) [0]
        self.root.parent = self.root.name
        self.total_weight = float ( sum( size for size in tumour.values() ) )
        
        # assign to each gen-tree-vertex a unique number
        number = 0
        for v in self.V.values():
            self.type_of_number.append(v.name)
            v.number = number
            number += 1
        self.type_amount = number
        
        # finish by setting pointer of children and parents to actual nodes
        for v in self.V.values():
            par = v.parent
            v.parent = self.V[par]
            for i in range (len (v.children)):
                child = v.children[i]
                v.children[i] = self.V[child]
                
    def reduce_tree(self, subset_dict: dict):
        """ eliminates all nodes that are not part of the subset-dict
            weight-rebalancing has been done in get_hierarchy_sets
            in this step, only vertices are removed """

        subset = subset_dict
            
        # update V
        self.V = subset
        intersect = set( subset_dict.values() )
        
        # update generations
        for i in range(len(self.generations)):
            gen = self.generations[i]
            gen.intersection_update( intersect )
            self.generations[i] = gen
        
        i = max( [i for i in range(len(self.generations)) if len(self.generations[i]) > 0 ], default = 0 )
        self.generations = self.generations[:i+1]
            
        # update children
        for v in self.V.values():
            #print("subset", intersect)
            children = set( v.children )
            #print("before", children)
            children.intersection_update(intersect)
            
            v.children = list( children )
            #print("after", v.children)
            #raise ValueError()
        
        
        self.type_of_number = [] # internal only
        self.type_amount = None     # number of types
        
        self.group_points = None    # nodes from clustering
        self.unclassified_maximal_fitness = None
        
        # assign to each gen-tree-vertex a unique number
        number = 0
        for v in self.V.values():
            self.type_of_number.append(v.name)
            v.number = number
            number += 1
        self.type_amount = number
        
        # update height
        self.height = len(self.generations) 
        
           
    def branch_based_cluster_SIMPLE(self, cutoff = 0.05):
        """ Essentially removes all mutations with lof CCF
            All other mutations are put into genotype-clusters (group-points for further usage) """
        
        # delete rare subtypes = crop the tree
        vertices_to_keep = {k:v for k, v in self.V.items() if float(v.branch_weight)/self.total_weight >= cutoff   }
        self.reduce_tree(vertices_to_keep)
        
        # redefine the weights of the nodes after the cropping
        for k,v in self.V.items():
            children_weight = sum ( w.branch_weight for w in v.children )
            v.weight = v.branch_weight - children_weight
          
        group_points = set()
        # assign values to all vertices with positive weight
        for v in self.V.values():
            if v.weight > 0:
                v.group = v
                group_points.add(v)
            else:
                v.group = None
        
        self.group_points = OrderedSet(group_points)
                     
        
    def branch_based_cluster(self):
        """ Clusters the vertices of the tree
            Bottom-up search: points beneath the last branching
                are defined as clusters
            WARNING: POINTS ABOVE THE LEAVES ARE NOT CLUSTERED """
        group_points = set()
        for v in self.V.values():
            v.group = None
        
        # recursive looking for the nearest branching-point
        def look_for_group(self, v):
            if len(v.parent.children) == 1:
                if v != self.root:
                    marker = look_for_group(self,v.parent)
                if v == self.root:
                    marker = v
            if len(v.parent.children) > 1:
                marker = v
            # take the root as special case
            elif len(v.parent.children) == 0:
                if self.root == v:
                    marker = v
                else:
                    raise ValueError("Clustering: Empty parent -> children array")
            v.group = marker
            group_points.add(marker)
            return marker
        
        for i in range (self.height -1, -1, -1):
            # bottum-up search
            gen = self.generations[i]
            for v in gen:
                # define the group if the vertex is not grouped yet
                if v.group == None:
                    for child in v.children:
                        if child.group != None:
                            v.group = False
                if v.group == None:
                    look_for_group(self, v)
                    
                    
        self.group_points = OrderedSet(group_points)
        self.unclassified_maximal_fitness = \
                max ( [ v.fitness for v in self.V.values()
                    if v.weight > 0 and v.group == False  ], default = 0  )
    
    
    def set_colours(self):
        """ defines colours for the clonal groups """
        group_list = list(self.group_points)
        colours = [ 1.0 / len(group_list) * i for i in range(len(group_list)) ]
        random.shuffle(colours)
        for i in range (len(group_list)):
            group_list[i].colour = colours[i]
    
    
    ###########################
    def cluster_print(self, filename="test", title="", fullpathname=False, plotname=False, *additional_markers):
        """ Input: G, colours, group_points=list(cluster_points), filename """
        # take a look at type_clustering_print.cluster_plot
        # plots a gen-tree
        
        
        group_points = self.group_points
        group_list = list(group_points)
        # set coordinates
        subtree_size = {}
        for i in range ( self.height -1 , -1 , -1) :
            gen = self.generations[i]
            for v in gen:
                subtree_size[v] = max (
                        sum ( subtree_size[child] for child in v.children )
                        ,  1 )
                
        #assign (x,y) values to nodes
        coordinates = {}
        root_position = 0   # common ancestor is placed in the center
        root = self.root
        
        # creates colours if necessary
        if len(group_list) > 0:
            if group_list[0].colour == None:
                self.set_colours()
                
        # create new colours only if necessary
        if all ( [v.colour == None for v in self.group_points]   ):
            self.set_colours()
        colours = list ( v.colour for v in self.group_points  )
        
            
        #Finds position of the nodes, recursive node-placement
        #Definces coordinates-dictionary
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
                position += (subtree_size[child] * 0.5 * newfactor)
                place_subtree(self, child, position, factor)
                position += (subtree_size[child] * 0.5 * newfactor)
            
        place_subtree(self, root, root_position, 1)
        
        # find edge-positions for the tree 
        # add green markers to indicate living traits
        x_edges = []
        y_edges = []
        
        living_x = []
        living_y = []
        percent_living = []
        
        
        for v in self.V.values():
            start = coordinates[v]
            for child in v.children:
                end = coordinates[child]

                x_edges.append(start[0])
                x_edges.append(end[0])
                
                y_edges.append(start[1])
                y_edges.append(end[1])
                
                y_edges.append(None)
                x_edges.append(None)
            # extra marker for living traits with size larger 0.5%
            percent_val = round ( v.weight / self.total_weight * 100, 2)
            if percent_val > 0 :
                percent_val = round(percent_val, 2)
                percent_living.append ( str(percent_val) +"%" )
                    
                living_x.append ( start[0] )
                living_y.append ( start[1] )
                
        x_clustering = []
        y_clustering = []
        text_clustering = []
        strength_of_clustering = str ( int ( sum ( v.branch_weight for v in group_list ) \
                    / self.total_weight * 100) ) + "%"
        name_clustering = strength_of_clustering + \
                    " are clustered"
        
        for v in group_list:
            rel_weight = round ( v.branch_weight /self.total_weight *100,2)
            x_clustering.append ( coordinates[v][0] )
            y_clustering.append ( coordinates[v][1] )
            
            M = v.max_fitness
            MT = v.fittest_type
            
            
            text_clustering.append ( "Size:  " + str ( rel_weight ) +"%" +\
                                    " | Type:  " + str( v.name ) +
                                    " | Mean fitness:  " + str ( round (v.mean_fitness,2))
                                    + " <br /> " +
                            "Maximal fitness:  " + str( round (M, 2)) +\
                            " | Fittest type:  " + str (MT) 
                            )
            
        ##########################
        ##########################
        # if additional clustering results are given, add them to the figure
        # name should be changed if needed. default: types not found in biopsy
        if len(additional_markers) > 0:
            not_found_types = additional_markers[0]
            if len(additional_markers) >1 :
                print( len(additional_markers))
                raise ValueError("cluster_print failed, additional markers error")

            x_extra = []
            y_extra = []
                
            num = 0
            for v in list(not_found_types):
                num += v.branch_weight
                x_extra.append ( coordinates [v][0] )
                y_extra.append ( coordinates [v][1] )
                
            extra_element = go.Scatter(x= x_extra,
                              y= y_extra,
                              hoverinfo = 'none',
                              mode = 'markers',
                              marker = dict(size = 30, opacity = 0.9, color="black"),
                              name="Not found in biopsy"
                              )
        ##########################        
        ##########################

        # black edges between coalescants
        edges           =go.Scatter(
                            x=x_edges,
                            y=y_edges,
                            mode='lines',
                            name='edges can also indicate <br />' +
                                        'connections via extinct types',
                            line=dict( color='black', width = 3),
                            hoverinfo='none'
                            )
        
        # living traits have small green dots
        living_types    = go.Scatter (x= living_x,
                              y= living_y,
                              hoverinfo = 'text',
                              text = percent_living,
                              opacity = 0.9,
                              mode = 'markers',
                              name='Living genotypes',
                              marker = dict(
                                size = 12,
                                color = 'green' ) )
        # cluster nodes with colours as in the 3d-plot
        cluster_plot = go.Scatter(x= x_clustering,
                          y= y_clustering,
                          hoverinfo = 'text',
                          text = text_clustering,
                          mode = 'markers',
                          marker = dict(size = 25, opacity = 0.99,
                                                    colorscale = 'Jet',
                                                    cmin = 0,
                                                    cmax = 1,
                                                    color=colours),
                          name=name_clustering
                          )
 
        if len(colours) != len(x_clustering):
            raise ValueError( str("length of colours does not match, show:\n" +\
                  len(colours) + ", " + len(group_list))  )
        
        # put things together to a plotly figure and print it
        figure = {'data': [], 'layout': {} }
        noaxis = dict           (
                            autorange=True,
                            showgrid=False,
                            zeroline=False,
                            showline=False,
                            ticks='',
                            showticklabels=False )

        layout = go.Layout( xaxis = noaxis, yaxis = noaxis, title=title,
                           hovermode = 'closest')
        
        figure['data']=[edges, living_types, cluster_plot]
        if len(additional_markers) > 0:
            figure['data'].append(extra_element)
            
        figure['layout']= layout

        # Creates directory for storage
        if not fullpathname:
            dirName = EVAL_data_output.create_folder(filename)
        else:
            dirName = filename
            Path(dirName).mkdir(parents=True, exist_ok=True)
            
        if not plotname:
            name = str ( dirName + "/" + "anc_tree_clustered.html" )
        else:
            name = str ( dirName + "/" + plotname + ".html" )

        py.plot(figure, filename = name, auto_open=False ) 
        return
    
    def cluster_and_biopsy_print(self, filename, fullpathname = False, plotname = False):
        # special case using cluster_print for evaluating the biopsy
        not_found_clones = [ v for  v  in self.group_points if not v.is_found ]
        perc_found = int (100 - len(not_found_clones) / float(len(self.group_points))*100 )
        perc_found = str(perc_found) + "%"
        AN = "ANC: " + str(len(self.root.name)) + " | f " + str(self.root.fitness)
        MA = "MAX: " + str( max( len(key) for key in self.V.keys() ) ) +\
                " | f " + str( round( max( [v.fitness for v in self.V.values() ] ) , 2)  )
        
        title = "Clonal hierarchy and Clustering<br />" + \
                "Biopsy found " + perc_found + " of clustered types<br />" + \
                AN + "<br />" + MA
        self.cluster_print(filename, title, fullpathname, plotname, not_found_clones)


#%%

def create_gen_graph(G, find_size=1, find_type="absolute"):
    """ Input:  Graph G
        Returns a phylogenetic tree object
    
    if find_type == "relative":
                find_size gives in % the minimal branch_weight of the resulting graph
    if find_size == "absolute":
                find_size gives the minimal absolute branch_weight """
    
    if G.count <= find_size and find_type == "absolute":
        raise AssertionError("Cluster size ", find_size, " too large for given graph ", G.count)
    
    # get ancestral graph and additional information
    tumour, ANC, generations, branch_weight, centroids = \
            eval_gr.preparations(G)
    
    # preparations: delete small types
    # let pack_graph_infos do the rest  
    
    # sets minimal size-value for a "clone"
    total_weight = float ( len(G.V) - 1 )
    if find_type == "absolute":
        critical_size = find_size
    elif find_type == "relative":
        critical_size = int(find_size * 0.01 * total_weight) +1
    else:
        raise ValueError (" Cluster-type not specified: \
                              find_type must be absolute or relative")
    
    # assign mean-fitness to branches (=mutations)
    # and deletes mutations with too small num of cells
    mean_fitness = defaultdict(int)
    max_fitness = defaultdict(int)
    fittest_type = defaultdict(int)
    
    # assign mean-fitness-values
    for INDEX in range ( len(generations) -1, ANC-1, -1) :
        gen = generations[ INDEX ]
        for TYP in gen:
            if len( gen[TYP] ) >= 1 :
                mean_fitness[TYP] = sum ( 
                        mean_fitness[child] * (branch_weight[child] )
                        for child in gen[TYP] ) / float(branch_weight[TYP]) 
                
                #bugfix, max-function seems does not seem to work :S
                index = 0
                temp = max_fitness[gen[TYP][index]]
                for i in range (1,len(gen[TYP])):
                    poss_child = gen[TYP][i]
                    if max_fitness[poss_child] > temp:
                        index = i
                        temp = max_fitness[gen[TYP][i]]
                        
                max_fitness[TYP] = temp
                fittest_type[TYP] = fittest_type[gen[TYP][index]]
            else:
                max_fitness[TYP] = G.fitness_book[TYP]
                fittest_type[TYP] = TYP
            
            if TYP in tumour:
                mean_fitness[TYP] += \
                    G.fitness_book[TYP] * (tumour[TYP] / float(branch_weight[TYP]))
     

    # delete all types which are too small, below the critical size,
    #    adjust the tumour-weights 
    for INDEX in range ( len(generations) -1, ANC-1, -1) :
        gen = generations[INDEX]
        del_list = [ key for key in gen
                    if branch_weight[key] < critical_size  ]
        for TYP in del_list:
            del gen[TYP]
            #now throw out the entry in the parent-node and adjust its weight
            ancestor =  tuple(TYP[:-1]) 
            generations[ INDEX - 1 ][ancestor].remove(TYP)
            tumour[ancestor] += branch_weight[TYP]
            tumour.pop(TYP, None)
    
    # check the weights     
    a = generations[ANC]
    for TYP in a:
        root = TYP
    
    if total_weight != branch_weight[root] :
        print(total_weight, branch_weight[root])
        raise ValueError("Preprocess gen-graph Error, too much deletion")
     
    Tree = ancestral_tree()
    Tree.pack_graph_infos(G, tumour, ANC, generations, branch_weight, centroids, mean_fitness, max_fitness, fittest_type)
    
    Tree.strong_mutations_log = copy.deepcopy(G.strong_mutations_log)
    
    return Tree
