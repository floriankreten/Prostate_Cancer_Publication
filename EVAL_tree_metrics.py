#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:26:45 2024

@author: florian
"""

import copy
import math
from sortedcontainers import SortedSet

""" Metrics for evaluation of weighted phylogenetic trees """


def _crop_and_set_weights(untreated_tree, cutoff = 0.05):
    """ Remove rare subtypes with rel_freq < cutoff 
        Remove vertices with weight 0 that dont define a branching
        Internal preprocessing, use the cropped tree only for metric-calculation """
    
    tree = copy.deepcopy(untreated_tree)
    
    # delete rare subtypes
    vertices_to_keep = {k:v for k, v in tree.V.items() if float(v.branch_weight)/tree.total_weight >= cutoff   }
    tree.reduce_tree(vertices_to_keep)
    
    
    # redefine the weights of the nodes after the cropping
    for k,v in tree.V.items():
        children_weight = sum ( w.branch_weight for w in v.children )
        v.weight = v.branch_weight - children_weight
        
    group_points = set()
    # assign values to all vertices with positive weight
    for v in tree.V.values():
        if v.weight > 0:
            v.group = v
            group_points.add(v)
        else:
            v.group = None
            
    def __sort_order(item):
        return len(item.name)
    
    new_order = SortedSet( group_points, key = __sort_order )
    
    tree.group_points = new_order
    for v in tree.group_points:
        print(v.name)
                     
    
    # checksum for the new weights 
    assert  sum ( v.weight for k,v in tree.V.items() )  == tree.total_weight, "Tree cropping failed"
    
    return tree
    
#%%

def _Diversity_Index(tree):
    """ Input: cropped tree, see ""crop_and_set_weights"
        Output: Diversity Index , aka inverse Simpson Index D>=1 """
    return ( 1 / sum(  (float(v.weight)/tree.total_weight)**2 for k,v in tree.V.items() ) )


#%%

def _balance_index(tree):
    """ Input: cropped tree, see ""crop_and_set_weights"
        Output: Balance Index J_1c, see Noble et al. 
        For non-conservative Index J_1, comment out #1 """
    
    internal_nodes = [v for k,v in tree.V.items() if len(v.children) >= 1]
    
    # weighted sum over branchings
    WS = 0.0
    # normalization constant
    C = 0.0
    
    for v in internal_nodes:
        Total_Child_Size = float( v.branch_weight - v.weight )
        C = C + Total_Child_Size
        
        d = len(v.children)
        if d >= 2: # spotted a branching, only those give score >0
            W_v = 0.0
            for w in v.children:
                P = w.branch_weight / Total_Child_Size
                r = - P * math.log(P, d)
                W_v = W_v + r
                
            # Weighted branching score for vertex v
            v_score = W_v * Total_Child_Size
            
            # J_1c modification (Can also be commented out, see Noble)
            v_score = v_score * Total_Child_Size / v.branch_weight #1
            WS = WS + v_score
            
            
    # normalization constant
    C = float( sum ( v.branch_weight - v.weight for v in internal_nodes   ) )
    
    return WS/C

#%%

def calculace_DiversityD_and_BalanceJ(tree, cutoff = 0.05):
    cropped_tree = _crop_and_set_weights(tree, cutoff)
    D = _Diversity_Index(cropped_tree)
    J = _balance_index(cropped_tree)
    return D,J

