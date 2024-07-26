#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:11:15 2019

@author: floriankreten
"""

""" Preferably use this module via EVAL_single_graph_evaluation """

""" evaluates a graph by using a gen_tree_graph and creates a
    3d map of the different clones
    classify nodes produces an input for a subset-plot (see graph_plot_plotly) """

import EVAL_gen_tree_graph
import EVAL_data_output
import EVAL_tree_metrics


#%%
def classify_nodes ( G, find_size, find_type):
    """ Creates a tree with min_size find_size
        Clusters this tree
        Classifies the vertices of the graph """
    tree = EVAL_gen_tree_graph.create_gen_graph(G, find_size, find_type)
    tree.branch_based_cluster_SIMPLE()
    tree.set_colours()
    
    # G, group_list, tree = normal_classification_routine(G,tree)
    G, group_list, tree = SIMPLE_classification_routine(G,tree)
                
    return G, group_list, tree


def classify_nodes_from_directory(G, dirname, find_size=0.05, find_type="relative",
                                  cluster_type = "simple_branched"):
    """ Reads in a stored tree from a given directory
        Transforms this into a rough tree
        Clusters are defined as living traits
        Classifies vertices of the graph """
        
    tree = EVAL_data_output.get_tree_sql(dirname, fullpathname=True)
    
    # includes a simple branch-based clustering:
    tree = EVAL_tree_metrics._crop_and_set_weights(tree, cutoff = find_size)
    
    print("Tree loaded and cropped, simple_branched genotypes defined")
        
    
    tree.set_colours()
    G, group_list, tree = SIMPLE_classification_routine(G,tree)
        
                
    return G, group_list, tree
        
    
#%%
def SIMPLE_classification_routine(G,tree):
    """ given a cropped tree, a classification is done
        either: 1) group is assigned if the vertex has positive weight in the branch
        2) or group is assigned by going up the tree until such a group is found
        3) no group is assigned, which is a super rare bug i could not fix yet (?)
    """
    print("Classifying nodes of the graph for printing")
    
    # ordered in "ascending" order (from ancestor to leaves by length)
    cluster_nodes = tree.group_points
    group_list = []
    for v in tree.group_points:
        group_list.append(v.name)
    group_list = tuple(group_list)
    
    group_numbers = {}
    i=0
    for v in cluster_nodes:
        group_numbers[v.name] = i
        i+=1
    

    for v in G.V[1:]:
        v.cluster_group = 0
        real_genotype = v.trait
        
        for i in range ( len(real_genotype)  ):
            comparison = real_genotype[:i]
            if comparison in group_numbers.keys():
                v.cluster_group = group_numbers[comparison]
                #break
                
    return G, group_list, tree    