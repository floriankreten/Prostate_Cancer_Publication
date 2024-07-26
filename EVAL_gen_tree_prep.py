#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 13:57:53 2019

@author: floriankreten
"""

""" Some technical functions used when creating ancestral trees
    Typically, you can ignore them
"""

from collections import defaultdict


def tube_sample(G, diam, point):
    """ tubula biopsy-sample with mid-point in x1=point, x2=0
        the sample is a trait-dictionary with weights"""
        
    radius = (diam/2)**2
    tube_sample = defaultdict(int)
    
    for v in G.V[1:]:
        pos = v.position
        if (pos[0]-point)**2 + pos[1]**2 <= radius:
            tube_sample[v.trait] += 1
    
    return tube_sample

def preparations(G):
    """ Input: Graph G
        Does some prepararional work for the output
        Returns: Tumour (dictionary of weights of living types)
        Ancestral Graph
        branch-weight (dictionary of weights up to MRCA) """
    
    #creates dictionary of living traits: standard-size = 0
    tumour = defaultdict(int)
    
    # While scanning, also define centroids
    # Starts with 3d-coordinates [0,0,0]
    living_centroids = defaultdict( lambda: [0,0,0] )
    
    total_weight = G.count
    for v in G.V[1:]:
        trait = v.trait
        tumour[trait] += 1
        living_centroids[trait] = PLUS_arr(living_centroids[trait], v.position)
    
    # generates trait-space, alphabetically sorted by generation
    tumour_list = list(tumour.keys())
    tumour_list.sort ( reverse = True)
    tumour_list.sort (key= len, reverse=True )
    
    Ancestor_generation, generations, branch_weight, centroids = \
                find_ancestral_lines(tumour_list, tumour, living_centroids)

    
    #check for the MCA, look at split-points and weight of the subtrees
    try:           
        while len( generations[Ancestor_generation+1] ) == 1:
            g = generations[Ancestor_generation+1]
            for TYP in g:
                poss_anc = TYP
            if branch_weight[poss_anc] == total_weight:
                Ancestor_generation+=1
            else:
                break
    except IndexError:
        pass
    
    
    return tumour, Ancestor_generation, generations, branch_weight, centroids


def find_ancestral_lines(tumour_list, tumour, living_centroids, *ignore_weights):
    """ IN: Sorted tumour as in eval_clonal_structure
        OUT: Generation of common ancestor and sorted generations-graph"""
    
    # by graph search: find ancestors with living predecessors
    # sort types by generations,
    # value shows number of living predecessor-lines
    K = max ( len(TYP) for TYP in tumour_list  )
    generations = [ {} for i in range(K+1)]
    working_list = [ list(typ) for typ in tumour_list  ]
    
    #sort list into generations, returns L as oldest living generation
    if len(working_list) > 0:
        # now sort into the generations
        # in decending order
        for typ in working_list:
            LEN = len(typ)
            TYP = tuple(typ)
            if TYP not in generations[LEN]:
                generations[LEN][TYP] = []  # pointing to all known children
                                            # ancestor is known by own name
    else:
        raise ValueError("Error for ancestral lines: Empty Type-List")
    
    L = min(len(TYP) for TYP in tumour_list)
    
    # Trace back all living populations to their ancestors
    # the entry of every type points to the list of children
    for Index in range (K, L-1, -1):
        
        for TYP in generations[Index]:
            ancestor =  tuple(TYP[:-1])
            if ancestor in generations[Index -1]:
                generations[Index -1][ancestor].append(TYP)
            else:
                generations[Index -1][ancestor]  = [TYP]
    
    #now also trace back the "dead"/old part of the tree to find the common ancestor
    ANCESTOR = 0        
    for Index in range ( L-1,  0, -1  ):
        for TYP in generations[Index]:
            ancestor =  tuple(TYP[:-1])
            if ancestor in generations[Index -1]:
                generations[Index -1][ancestor].append(TYP)
            else:
                    generations[Index -1][ancestor]  = [TYP]
        if len(generations[Index]) == 1:
            ANCESTOR = Index
            break
        elif len(generations[Index]) ==0:
            raise ValueError ("Ancestor finding mistake")
            break
    
    # iteratively create centroids and branch_weights (bottom-up)
    centroids = {}    
    branch_weight = {}
    if not ignore_weights:
        for INDEX in range ( len(generations) -1, ANCESTOR-1, -1) :
            gen = generations[INDEX]
            for TYP in gen:
                if TYP in tumour:
                    branch_weight[TYP] =  tumour[TYP]
                    centroids[TYP]     =  living_centroids[TYP]
                else:
                    branch_weight[TYP] = 0
                    centroids[TYP]     = [0,0,0]
                #get weight of total subtree
                for child in gen[TYP]:
                    branch_weight[TYP] += branch_weight[child]
                    centroids[TYP] = PLUS_arr(centroids[TYP], centroids[child])
     
    # divide centroids by branch_weight for coordinate-center and finalize         
    for TYP in branch_weight:
        W = float(branch_weight[TYP])
        centroids[TYP] = tuple( [x / W for x in centroids[TYP] ])
               
    return ANCESTOR, generations, branch_weight, centroids
            

#%%

def is_found(generations, ANC, sample):
    """ In: generation, ancestor-gen, sample
        Out: Set of found mutations
        returns True if a mutation is found in the sampling process
        returns True iff a child with this mutation was in the sample
        (e.g. useful for biopsies)"""
        
    found = set()
    for INDEX in range ( len(generations) -1, ANC-1, -1) :
        gen = generations[ INDEX ]
        for TYP in gen:
            if TYP in sample:
                found.add(TYP)
            else:
                for child in gen[TYP]:
                    if child in found:
                        found.add(TYP)
                        break
    
    return found

def find_generations(G):
    """ finds the generations of the living traits of a graph G """
    #creates dictionary of living traits
    tumour = defaultdict(int)
    for v in G.V[1:]:
        trait = v.trait
        tumour[trait] += 1
    
    # generates trait-space, needs some sorting for creating generation-graph
    tumour_list = list(tumour.keys())
    tumour_list.sort ( reverse = True)
    tumour_list.sort (key= len, reverse=True )
    
    Ancestor_generation, generations, branch_weight = \
                find_ancestral_lines(tumour_list, tumour, "ignore_weights")
    try:
        while len( generations[Ancestor_generation+1] ) ==1:   #get rid of weird bug
            Ancestor_generation+=1
    except IndexError:
        pass
            
    return Ancestor_generation, generations

def PLUS_arr(a,b):
    """Sum of two arrays"""
    k = len(a)
    if k != len(b):
        raise ValueError("Länge", a, "!= Länge", b)
    else:
        c=[]
        for i in range(0, k):
            c.append(a[i]+b[i])
        return c