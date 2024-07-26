#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 15:33:47 2019

@author: floriankreten
"""

""" Preferably use this module via EVAL_single_graph_evaluation """

""" Two different types of biopsies are defined: circular and "punching"
    A biopsy is a list of sets, where each set is one "spatial component"
    of the biopsy
"""

import random
import math
import os
import gc

import pandas as pd
import numpy as np

import EVAL_graph_plot_plotly
import EVAL_gen_tree_graph
import EVAL_type_clustering_print

import EVAL_data_output


def take_4_quad_on_different_levels_biopsy(G, target_size=500,
            heights=[0.35, 0.65], num_samples = [2,2], ratio_diam_h = 15):
    
    # check if input is ok (needed since the format is standardized)
    for num in num_samples:
        if num != 2:
            raise ValueError("Num_samples for level_biopsy must be 4 (currently)")
    for h in heights:
        if not 0 <= h <= 1:
            raise ValueError("Heights must be normed to be in [0,1]")  
    if len(heights) != len(num_samples):
        raise ValueError("Dimension of num_samples does not match heights:\n", \
                         num_samples, heights)
    if len(heights) != 2:
        raise ValueError("Len heights must be equal to 1 atm")
        
    # FIRST STEP: understand height of the graph
    
    zmax = max ( v.position[2] for v in G.V[1:] )
    zmin = min ( v.position[2] for v in G.V[1:] )
    
    #diam = min (xmax, -xmin, ymax, -ymin, zmax, -zmin)
    z_radius = (zmax - zmin) / 2
    
    total_samples = []
    
    for i in range(len(heights)):
        # SECOND STEP: estimate local density
        height = heights[i]
        no_samples = num_samples[i]
        total_target_size = no_samples * target_size
        print("target_size: ", target_size, " num biopsies: ", no_samples)
        print("total target-size: ", total_target_size)
        
        z_level = zmin + 2* z_radius * height
        
        thickness = 4 * G.parameters["vertex_radius"]
            
        level_set = set ( v for v in G.V[1:]
                              if dist_height(v, z_level) <= 0.5 * thickness )
            
        #print("Size of first total level_set:", len(level_set))
            
        xmax = max ( v.position[0] for v in level_set )
        xmin = min ( v.position[0] for v in level_set )
                    
        ymax = max ( v.position[1] for v in level_set )
        ymin = min ( v.position[1] for v in level_set )
            
        xy_center = ( (xmax+xmin) * 0.5, (ymax+ymin) * 0.5 )
        
        xy_radius = min ( xmax-xmin, ymax-ymin  ) / 2
        xy_allowed_sup_from_middle = xy_radius / math.sqrt(2)
        xy_diam = 2 * xy_allowed_sup_from_middle
        
        
        # THIRD STEP: define size of biopsies, check if this fits with the graph
        max_quadrant = set (  v for v in level_set if
                              sup_dist_in_plain(v,xy_center) <= xy_allowed_sup_from_middle )
        
        print("Size of first max quadrant: ", len(max_quadrant))
        
        local_density = len(max_quadrant) / ( (2*xy_allowed_sup_from_middle)**2 * thickness )
        
        # size of a biopsy (length_x*length_y*height_z)
        target_height = ( target_size / 
                             (local_density * (ratio_diam_h**2)  ) ) ** (1.0/3.0)
        target_length = target_height * ratio_diam_h
        
        if target_length*0.9 >= xy_diam:
            raise ValueError("Graph not big enough for target size", target_size)
            
        est_size = local_density * (target_length**2) * target_height
        #print("Estimated size of one quadrant: ", est_size)
        
        # FOR BIOPSIES WITH HORIZONTAL DISTANCE, DELETE THE FOLLOWING LINE
        xy_allowed_sup_from_middle = target_length
        
        # FOURTH STEP: take two biopsies on the opposite site of the graph
        lower_left = ( xy_center[0] - xy_allowed_sup_from_middle, xy_center[1] - xy_allowed_sup_from_middle)
        upper_right = ( xy_center[0] + xy_allowed_sup_from_middle, xy_center[1] + xy_allowed_sup_from_middle)
        
        # separate the two boxes height if needed
        if target_length >= 0.5 * xy_diam:
            height_up = z_level + 0.5 * target_height
            height_down = z_level - 0.5 * target_height
        else:
            height_up = z_level
            height_down = z_level
        
        level_set_up = set ( v for v in G.V[1:]
            if  (dist_height(v, height_up) <= 0.5 * target_height)
            and (sup_dist_in_plain(v,xy_center) <= xy_allowed_sup_from_middle) )
        
        biopsy_up = set ( v for v in level_set_up 
            if sup_dist_in_plain(v, lower_left) < target_length)
        
        level_set_down = set ( v for v in G.V[1:]
            if  (dist_height(v, height_down) <= 0.5 * target_height)
            and (sup_dist_in_plain(v,xy_center) <= xy_allowed_sup_from_middle) )
        
        biopsy_down = set (v for v in level_set_down
            if sup_dist_in_plain(v, upper_right) < target_length )
        
        total_samples.append(biopsy_up)
        total_samples.append(biopsy_down)
        
        
    # LAST STEP: Update the vertices of the graph with biopsy-info
    G = update_graph_with_biopsy_info(G, total_samples)
        
    return G, total_samples


def take_4_quad_on_same_level_biopsies(G, target_size=500, heights= [0.5],
                             num_samples = [4], ratio_diam_h = 15):
    
    # check if input is ok (needed since the format is standardized)
    for num in num_samples:
        if num != 4:
            raise ValueError("Num_samples for level_biopsy must be 4 (currently)")
    for h in heights:
        if not 0 <= h <= 1:
            raise ValueError("Heights must be normed to be in [0,1]")  
    if len(heights) != len(num_samples):
        raise ValueError("Dimension of num_samples does not match heights:\n", \
                         num_samples, heights)
    if len(heights) != 1:
        raise ValueError("Len heights must be equal to 1 atm")
        
    # FIRST STEP: understand height of the graph
    
    zmax = max ( v.position[2] for v in G.V[1:] )
    zmin = min ( v.position[2] for v in G.V[1:] )
    
    #diam = min (xmax, -xmin, ymax, -ymin, zmax, -zmin)
    z_radius = (zmax - zmin) / 2

    height = heights[0]
    no_samples = num_samples[0]
    total_target_size = no_samples * target_size
    #print("target_size: ", target_size, " num biopsies: ", no_samples)
    #print("total target-size: ", total_target_size)
    
    z_level = zmin + 2* z_radius * height
    
    
    # SECOND STEP: estimate local density
    thickness = 4 * G.parameters["vertex_radius"]
        
    level_set = set ( v for v in G.V[1:]
                          if dist_height(v, z_level) <= 0.5 * thickness )
        
    #print("Size of first total level_set:", len(level_set))
        
    xmax = max ( v.position[0] for v in level_set )
    xmin = min ( v.position[0] for v in level_set )
                
    ymax = max ( v.position[1] for v in level_set )
    ymin = min ( v.position[1] for v in level_set )
        
    xy_center = ( (xmax+xmin) * 0.5, (ymax+ymin) * 0.5 )
    
    xy_radius = min ( xmax-xmin, ymax-ymin  ) / 2
    xy_allowed_sup_from_middle = xy_radius / math.sqrt(2)
    xy_diam = 2 * xy_allowed_sup_from_middle
        

    # THIRD STEP: define size of biopsies, check if this fits with the graph
    max_quadrant = set (  v for v in level_set if
                          sup_dist_in_plain(v,xy_center) <= xy_allowed_sup_from_middle )
    
    #print("Size of first max quadrant: ", len(max_quadrant))
    
    local_density = len(max_quadrant) / ( (2*xy_allowed_sup_from_middle)**2 * thickness )
    
    # density * height * (2 *sup_dist)^2 = total_target_size
    # sup_dist = radio_diam_h * height
    # solve for height and check 
    # if there is enough space within the current graph
    
    target_height = ( total_target_size / 
                         ( 4 * local_density * (ratio_diam_h**2)  ) ) ** (1.0/3.0)
    target_length = target_height * ratio_diam_h
    target_total_diam_needed = 2 * target_length
    
    #print("target_total_diam: ", target_total_diam_needed )
    #print("local density: ", local_density)
    #print("current allowed diam: ", xy_diam)
    
    if target_total_diam_needed * 0.9 >= xy_diam:
        raise ValueError("Graph not big enough for 4 quadratic biopsies")
    
    # shrink size if needed
    if target_total_diam_needed >= xy_diam:
        target_length = xy_diam * 0.49
        
    else:
        print("Graph: geometry OKOK")
        
    est_size = local_density * (target_length**2) * target_height
    #print("Estimated size of one quadrant: ", est_size)
        
    # FOURTH STEP: get 4 quadratic biopsies
    total_samples = [ set(), set(), set(), set() ]
    
    # if you want biopsies that have a specified horizontal distance,
    # redefine xy_allowed_sup accordingly (maximal distance: no change!)
    xy_allowed_sup_from_middle = target_length
    
    four_corners = [
        ( xy_center[0] - xy_allowed_sup_from_middle, xy_center[1] - xy_allowed_sup_from_middle),
        ( xy_center[0] - xy_allowed_sup_from_middle, xy_center[1] + xy_allowed_sup_from_middle),       
        ( xy_center[0] + xy_allowed_sup_from_middle, xy_center[1] + xy_allowed_sup_from_middle),
        ( xy_center[0] + xy_allowed_sup_from_middle, xy_center[1] - xy_allowed_sup_from_middle)
        ]
    
    new_level_set = set ( v for v in G.V[1:] if
        (dist_height(v, z_level) <= 0.5 * target_height) )
    
    new_max_quadrant = set( v for v in new_level_set if
            sup_dist_in_plain(v, xy_center) <= xy_allowed_sup_from_middle )
    
    #print("Size of new max quadrant: ", len(new_max_quadrant))
    
    # get four quadratic biopsies within new_max_quadrant
    for v in new_max_quadrant:
        for i in range(4):
            if sup_dist_in_plain(v, four_corners[i]) < target_length:
                total_samples[i].add(v)
                break
    
    # FINAL STEP: Update the vertices of the graph with biopsy-info
    G = update_graph_with_biopsy_info(G, total_samples)
        
    return G, total_samples


def take_five_quad_circ_biopsies(G, target_size=500, heights=[0.5],
                             num_samples = [5], ratio_diam_h = 15):
    

    for num in num_samples:
        if num != 5:
            raise ValueError("Num_samples for level_biopsy must be 5 (currently)")
    for h in heights:
        if not 0 <= h <= 1:
            raise ValueError("Heights must be normed to be in [0,1]")
        
    if len(heights) != len(num_samples):
        raise ValueError("Dimension of num_samples does not match heights:\n", \
                         num_samples, heights)
    else:
        K = len(heights)

    # FIRST STEP: get size-parameters, understand graph (approximately a ball)
    
    zmax = max ( v.position[2] for v in G.V[1:] )
    zmin = min ( v.position[2] for v in G.V[1:] )
    
    #diam = min (xmax, -xmin, ymax, -ymin, zmax, -zmin)
    z_radius = (zmax - zmin) / 2
    total_samples = []
    
    
    # iterate over several "heights"
    for k in range(K):
        
        height = heights[k]
        num_sampl = num_samples[k]
        
        # SECOND STEP: get slice of the level-set at height z_level
        #              to understand the local density of the graph
        
        z_level = zmin + 2 * z_radius * height
        thickness = 4 * G.parameters["vertex_radius"]
        
        level_set = set ( v for v in G.V[1:]
                          if dist_height(v, z_level) <= 0.5 * thickness )
        
        #print("Size of first total level_set:", len(level_set))
        
        xmax = max ( v.position[0] for v in level_set )
        xmin = min ( v.position[0] for v in level_set )
                
        ymax = max ( v.position[1] for v in level_set )
        ymin = min ( v.position[1] for v in level_set )
        
        xy_center = ( (xmax+xmin) * 0.5, (ymax+ymin) * 0.5 )
        xy_radius = min ( xmax-xmin, ymax-ymin  ) * 0.45 # small dampening
        
        
        # estimate the density of the full triangle (with diam 2 max_diam)
        # use this to generate the final size parameters (if possible)
        
        
        full_circle = set (   v for v in level_set if
                          plain_dist2(v,xy_center) <= xy_radius**2 )
        
        #print("Size of first circle: " + str(len(full_circle )))
        density = len(full_circle) / ( math.pi * xy_radius**2 * thickness )
        
        #print("Density: ", density)
        
        
        # diameter = height * ratio_diam_h    is the "diameter" of one of the 5 biopsies
        # say diameter ~ radius of the outer circle (see output picture)
        # target_size * num_sampl = density * height * ( pi * r^2)
        
        height = ( (num_sampl * target_size) /  (ratio_diam_h**2 * density * math.pi ) ) **(1.0/3.0)
        radius = height * ratio_diam_h
        
        
        if radius > xy_radius:
            print("Needed radius" + str(radius) )
            print("Radius " + str(xy_radius)) 
            output = "Needed size " + str(target_size) + \
                    " can not be reached with ratio diam/h = " + str(ratio_diam_h)
            raise ValueError(output)
            
                
        adjusted_circ = set ( v for v in G.V[1:] if
                          plain_dist2(v,xy_center) <= radius**2 and
                          dist_height(v, z_level) <= 0.5 * height )
        
        print("Size of full level set: " + str(len(adjusted_circ)) )
        

        
        # THIRD STEP: take quadratic biopsies (if level_set is large enough)
        # Divide plane into 4 quadrants
        def quadrant(v):
            a = v.position[0] - xy_center[0] 
            b = v.position[1] - xy_center[1]
            if a >= 0:
                if b >= 0:
                    q = 0
                else:
                    q = 1
            else:
                if b >= 0:
                    q = 2
                else:
                    q = 3
            return q
        
        
        center_rad2 = target_size / (density * height * math.pi)
        center_set = set (v for v in adjusted_circ if
                           plain_dist2(v, xy_center) <= center_rad2  )
        
        # four quadratic sets
        samples = [ set(), set(), set(), set()]
        for v in adjusted_circ:
            if v not in center_set:
                q = quadrant(v)
                samples[q].add(v)
                
        samples.append(center_set)
                   
        for samp in samples:
            total_samples.append(samp)     
            print("Length of indiv. sample: " + str(len(samp)) )
            
        print("Density times height", density*height)
             
        
    # LAST STEP: Update the vertices of the graph with biopsy-info
    G = update_graph_with_biopsy_info(G, total_samples)
        
    return G, total_samples



#%%


def take_five_quadratic_biopsies_with_center_cut(G, target_size=500, heights=[0.5],
                             num_samples = [5], ratio_diam_h = 15):
    

    for num in num_samples:
        if num != 5:
            raise ValueError("Num_samples for level_biopsy must be 5 (currently)")
    for h in heights:
        if not 0 <= h <= 1:
            raise ValueError("Heights must be normed to be in [0,1]")
        
    if len(heights) != len(num_samples):
        raise ValueError("Dimension of num_samples does not match heights:\n", \
                         num_samples, heights)
    else:
        K = len(heights)

    # FIRST STEP: get size-parameters, understand graph (approximately a ball)
    
    zmax = max ( v.position[2] for v in G.V[1:] )
    zmin = min ( v.position[2] for v in G.V[1:] )
    
    #diam = min (xmax, -xmin, ymax, -ymin, zmax, -zmin)
    z_radius = (zmax - zmin) / 2
    total_samples = []
    
    
    # iterate over several "heights"
    for k in range(K):
        
        height = heights[k]
        num_sampl = num_samples[k]
        
        # SECOND STEP: get slice of the level-set at height z_level
        #              to understand the local density of the graph
        
        z_level = zmin + 2 * z_radius * height
        thickness = 4 * G.parameters["vertex_radius"]
        
        level_set = set ( v for v in G.V[1:]
                          if dist_height(v, z_level) <= 0.5 * thickness )
        
        print("Size of total level_set:", len(level_set))
        
        xmax = max ( v.position[0] for v in level_set )
        xmin = min ( v.position[0] for v in level_set )
                
        ymax = max ( v.position[1] for v in level_set )
        ymin = min ( v.position[1] for v in level_set )
        
        xy_center = ( (xmax+xmin) * 0.5, (ymax+ymin) * 0.5 )
        xy_radius = min ( xmax-xmin, ymax-ymin  ) * 0.48 # small dampening
        
    
        
        # estimate the density of the full level_set (with diam 2 max_diam)
        # use this to generate the final size parameters (if possible)
        max_diam = xy_radius / math.sqrt(2)
        full_rectangle = set (   v for v in level_set if
                          sup_dist_in_plain(v,xy_center) <= max_diam )
        print("Size of first rectangle: " + str(len(full_rectangle)))
        density = len(full_rectangle) / ( (2 * max_diam)**2 * thickness )
        
        height = ( (target_size * num_sampl) / ( density *  4 *ratio_diam_h**2  ) ) **(1.0/3.0)
        diam = height * ratio_diam_h
        rad = 2**(0.5) * diam
        

        if rad > xy_radius:
            raise ValueError("Needed size " + str(target_size) + \
                    " can not be reached with ratio diam/h = " + str(ratio_diam_h) )
                
        
        adjusted_box = set ( v for v in G.V[1:] if
                          sup_dist_in_plain(v,xy_center) <= diam and
                          dist_height(v, z_level) <= height/2 )
        
        print("Size of second, adjusted rectangle: " + str(len(adjusted_box)) )
        
        # THIRD STEP: take quadratic biopsies (if level_set is large enough)
        # Divide plane into 4 quadrants
        def quadrant(v):
            a = v.position[0] - xy_center[0] 
            b = v.position[1] - xy_center[1]
            if a >= 0:
                if b >= 0:
                    q = 0
                else:
                    q = 1
            else:
                if b >= 0:
                    q = 2
                else:
                    q = 3
            return q
        
        
        center_len = ( target_size / (density * height ) ) **(0.5)
        center_set = set (v for v in adjusted_box if
                           sup_dist_in_plain(v, xy_center) <= center_len * 0.5 )
        
        # four quadratic sets
        samples = [ set(), set(), set(), set()]
        for v in adjusted_box:
            if v not in center_set:
                q = quadrant(v)
                samples[q].add(v)
                
        samples.append(center_set)
                   
        for samp in samples:
            total_samples.append(samp)     
            print("Length of indiv. sample: " +str(len(samp)) )
            
        
    # LAST STEP: Update the vertices of the graph with biopsy-info
    G = update_graph_with_biopsy_info(G, total_samples)
        
    return G, total_samples

#%%

def take_chips_biopsy(G, target_size=500, heights=[0.5],
                             num_samples = [1], ratio_rad_h = 15):
        
    if len(heights) != len(num_samples):
        raise ValueError("Dimension of num_samples does not match heights:\n", \
                         num_samples, heights)
    else:
        K = len(heights)


    # FIRST STEP: get size-parameters, understand graph (approximately a ball)
    
    zmax = max ( v.position[2] for v in G.V[1:] )
    zmin = min ( v.position[2] for v in G.V[1:] )
    
    #diam = min (xmax, -xmin, ymax, -ymin, zmax, -zmin)
    z_radius = (zmax - zmin) / 2
    total_samples = []
    
    
    # iterate over several "heights"
    for k in range(K):
        
        height = heights[k]
        num_sampl = num_samples[k]
        
        z_level = zmin + 2 * z_radius * height
        thickness = 4 * G.parameters["vertex_radius"]
        
        # SECOND STEP: get slice of the level-set at height z_level
        level_set = set (   v for v in G.V[1:]
                          if dist_height(v, z_level) <= 0.5 * thickness )
        
        print("Size of total level_set:", len(level_set))
        
        xmax = max ( v.position[0] for v in level_set )
        xmin = min ( v.position[0] for v in level_set )
                
        ymax = max ( v.position[1] for v in level_set )
        ymin = min ( v.position[1] for v in level_set )
        
        xy_center = ( (xmax+xmin) * 0.5, (ymax+ymin) * 0.5 )
        xy_radius = min ( xmax-xmin, ymax-ymin  ) * 0.5
        
        def get_circ_subset(level_set, center, squared_radius):
            subset = set (   v for v in level_set if
                          plain_dist2( v, center) <= squared_radius
                              )
            return subset
        
        # estimate the needed radius (uniform for these samples)
        R_max = xy_radius * 0.98 # cut away the outmost 2%
        center_circle = get_circ_subset(level_set, xy_center, R_max**2)
        
        density = len(center_circle) / ( math.pi * R_max**2 * thickness)
        # Size of the biopsy depends on a single parameter, since the ratio
        # rad / height is given as an input. Uses the previously estimated density
        height =  ( target_size / ( density * math.pi * ratio_rad_h **2 ) ) **(1.0/3.0)
        r_pref = height * ratio_rad_h
        
        level_set = set (   v for v in G.V[1:]
                          if dist_height(v, z_level) <= 0.5 * height )
        
        gc.collect()
        
        print("Size of new level set", len(level_set))
        print("Density estimate", density)
        print("Height, r_pref", round(height,2), round(r_pref,2))
        

        # more of a test case than an actual sample
        if num_sampl == 1:
            sample = get_circ_subset(level_set, xy_center, r_pref**2  )
            total_samples.append(sample)   
            
        elif num_sampl >= 2:
            # check if the preferred radius is valid for the target-size
            # STOP IF THIS WAS NOT POSSIBLE-> Samples = []

            alpha = 2 * math.pi / num_sampl
            R_poss = R_max - r_pref # maximal possible (also best) spread
            if r_pref + 2 >= R_poss * math.sin(alpha/2):
                raise ValueError("Targeted biopsy-size can not be reached for \
                                 chosen height", height, "num_samples",num_sampl)

            beta = ( 2 * math.pi / len(heights) ) * k # rotate over heights
            for i in range(num_sampl):
                angle = alpha * i + beta
                X_rot = 1 * math.cos(angle) * R_poss
                Y_rot = 1 * math.sin(angle) * R_poss
                xc = xy_center[0] + X_rot
                yc = xy_center[1] + Y_rot
                sample = get_circ_subset(level_set, [xc,yc], r_pref**2  )
                total_samples.append(sample)
        gc.collect()
    
    
    
    # LAST STEP: Update the vertices of the graph with biopsy-info
    G = update_graph_with_biopsy_info(G, total_samples)
        
    return G, total_samples


#%%

def take_level_circle_biospy(G, target_size, heights=[0.5],
                             num_samples = [1], thickness = False):
    """ takes "num_samples" in [1,2,3,4] round biopsies
        at given "height", real-valued in [0,1], relative to graph height
        automatically tries to match the "target_size",
            that gives the number of vertices per biopsy"""
    
    for num in num_samples:
        if num not in [1,2,3,4]:
            raise ValueError("Num_samples for level_biopsy should be in [1,2,3,4]")
        
    if len(heights) != len(num_samples):
        raise ValueError("Dimension of num_samples does not match heights:\n", \
                         num_samples, heights)
    else:
        K = len(heights)


    # FIRST STEP: get size-parameters, understand graph (approximately a ball)
    xmax = max ( v.position[0] for v in G.V[1:] )
    xmin = min ( v.position[0] for v in G.V[1:] )
            
    ymax = max ( v.position[1] for v in G.V[1:] )
    ymin = min ( v.position[1] for v in G.V[1:] )
    
    zmax = max ( v.position[2] for v in G.V[1:] )
    zmin = min ( v.position[2] for v in G.V[1:] )
    
    #diam = min (xmax, -xmin, ymax, -ymin, zmax, -zmin)
    z_radius = (zmax - zmin) / 2
    total_samples = []
    
    # iterate over several "heights"
    for k in range(K):
        inner_scaling = 1
        height = heights[k]
        num_sampl = num_samples[k]
        
        z_level = zmin + 2 * z_radius * height
        z_from_zero = abs ( 2 * (0.5 - height) )
        x_rel_radius = math.sqrt(1 - z_from_zero**2) * \
                            inner_scaling * (xmax-xmin) * 0.5
        y_rel_radius = math.sqrt(1 - z_from_zero**2) * \
                            inner_scaling * (ymax-ymin) * 0.5                  
                            
        xy_rel_radius = min(x_rel_radius, y_rel_radius)
        
        xy_center = ( (xmax+xmin) * 0.5, (ymax+ymin) * 0.5 )
        
        if not thickness:
            thickness = 2 * G.parameters["vertex_radius"]
        
        # SECOND STEP: get slice of the level-set at height z_level
        def get_level_set(thickness, squared_radius):
            level_set = set (   v for v in G.V[1:]
                          if dist_height(v, z_level) <= thickness
                          and plain_dist2( v, xy_center) <= squared_radius
                              )
            return level_set
        
        # first try full level set
        level_set = get_level_set(thickness, float("inf"))
        
        # adjust size if necessary (and possible)
        adj_val = (target_size * num_sampl) / len(level_set)
        tries = 0
    
        while adj_val < 0.95 or adj_val > 1.05:
            if tries <6:
                inner_scaling = math.sqrt(adj_val) * inner_scaling
                curr_radius_sq = (xy_rel_radius * inner_scaling)**2
                level_set = get_level_set(thickness, curr_radius_sq)
                adj_val = (target_size * num_sampl) / len(level_set)
                tries = tries +1
            else:
                break
        if adj_val > 1.25 or adj_val < 0.8:
            print("Could not match required target size,",
                  "try to lower the height or increase thickness manually")
            print("Ratio target/current size:", round(adj_val,2))
        
        # THIRD STEP: divide level-set into samples
        if num_sampl == 1:
            samples = [level_set]
            
        elif num_sampl >= 2:
            
            # slice the x-y-plane into N (angular) pizza-pieces
            def vertex_angle_sort(v):
                Y = v.position[1]
                X = v.position[0]
                sliceno = np.int32((np.pi + np.arctan2(Y, X)) * (num_sampl / (2*np.pi)))
                return sliceno
            
            samples = []
            for i in range(num_sampl):
                samples.append( set() )
            for v in level_set:
                i = vertex_angle_sort(v)
                samples[i].add(v)
            for samp in samples:
                total_samples.append(samp)
        
        
    # LAST STEP: Update the vertices of the graph with biopsy-info
    G = update_graph_with_biopsy_info(G, total_samples)
        
    return G, total_samples


#%%
        

def take_circular_biopsy(G, diameter = 90.0, amount = 2 ):
    """ Input: Graph G
        takes a graph/tumour and does a biopsy
        returns n samples as lists of vertices"""

    if amount == 0:
        raise ValueError("Cannot take zero biopsies")    

    # first step: understand geometry of the graph
    # assumption: ~ ellipsoid
    
    xmax = max ( v.position[0] for v in G.V[1:] )
    xmin = min ( v.position[0] for v in G.V[1:] )
            
    ymax = max ( v.position[1] for v in G.V[1:] )
    ymin = min ( v.position[1] for v in G.V[1:] )
    
    zmax = max ( v.position[2] for v in G.V[1:] )
    zmin = min ( v.position[2] for v in G.V[1:] )
    
    rad_x = (xmax - xmin) / 2.0
    rad_y = (ymax - ymin) / 2.0
    rad_z = (zmax - zmin) / 2.0

    # approximated center of the tumour
    zero_of_graph = np.array ( [ (xmax+xmin) / 2.0, (ymax+ymin) / 2.0, (zmax+zmin) / 2.0 ] )
    
    # second step: choose n deterministic (x,y,z) points
    # as centres of the circular biopsies
    
    r_min = min(rad_x, rad_y, rad_z)
    if diameter >= r_min / 1.5:
        diameter = r_min / 1.5
        print("Graph too small, diameter of biopsy changed to {}".format(round(diameter,0)))
    radius_of_biopsy = diameter / 2.0
    
    angles = np.linspace(start = 0.0, stop = 2 * np.pi , num = amount, endpoint = False)
    centres = np.array ( [ [np.cos(alpha) * (rad_x - radius_of_biopsy),
                                np.sin(alpha) * (rad_y - radius_of_biopsy),
                                0.0] \
                             + zero_of_graph
                       for alpha in angles ] )
    if amount >= 4:
        print("Warning: There might be work to do with the geometry of the biopsies")

        
    
    # third step: take n samples
    samples = []
    
    for i in range(amount):
        check = radius_of_biopsy ** 2.0
        centre = centres[i]
        
        sample = set ( v for v in G.V[1:]
                            if dist2( v, centre ) <= check )
        samples.append(sample)
    
    # LAST STEP: Update the vertices of the graph with biopsy-info
    G = update_graph_with_biopsy_info(G, samples)
        
    return G, samples

def show_biopsy(G, samples, file_name, show = False, fullpathname = False):
    """ Use single_graph_eval.sketch_biopsy for full functionality"""
    
    """ Input: Graph G and samples/biopsies
        Plots the position of the biopsies within the graph """
    print("   Plotting samples:")
    print(" Size of biopsies:", [len(s) for s in samples]  )
    megaset = set().union(*samples)
    EVAL_graph_plot_plotly.subset_plot(G, megaset, file_name, show,
                                       title = "biopsy_sketch",
                                       fullpathname=fullpathname,
                                       biopsy_colour_code = True)


def compare_biopsy_and_clustering(samples, tree):
    """ Input:  Samples / biopsies and a genealogical tree corresponding
                to the same graph G (not needed as input)
        Evaluates which genotypes are found in a biopsy
        returns a set of all clusters which were not found """
    
    group_list = tuple(v.name for v in tree.group_points)
    
    is_found = {}
    for group in tree.group_points:
        is_found[group.name] = [0 for i in range(len(samples))]
    
    uncl = "unclassified"
    for j in range(len(samples)):
        sample = samples[j]
        for v in sample:
            index = v.cluster_group
            if index != uncl:
                is_found[group_list[index]] [j] += 1
    
    for g in tree.group_points:
        g.is_found = is_found[g.name]
        
    for g in tree.group_points:
        if g.is_found == [0, 0, 0, 0, 0]:
            g.is_found = False
                        
    return tree



def compare_biopsy_and_all_mutations(samples, tree):
    """ Input:  Samples / biopsies and a genealogical tree corresponding
                to the same graph G (not needed as input)
        Evaluates which mutations are found in a biopsy
        
        returns how many cells with a mutation
        are found within the biopsies"""
        
    for v in tree.V.values():
        v.is_found = [ 0 for k in range(len(samples)) ]
       
    # scan biopsies for found traits
    for k in range(len(samples)):
        sample = samples[k]
        for node in sample:
            try:
                genotype = tree.V[node.trait]
                genotype.is_found[k] += 1
            except:
                pass
            
    # bottom-up: add child-values for each node to the weight of the node
    for i in range (tree.height -1, -1, -1):
        gen = tree.generations[i]
        for v in gen:
            for child in v.children:
                for k in range(len(samples)):
                    v.is_found[k] += child.is_found[k]
            
    # check that everything is fine: only if going deep into the evaluation,
    # that is: no cutoff for the evaluated mutations. 
    # root-node has is_found values equal
    # to the total sizes of the biopsies
    # assert( tree.root.is_found == [ len(s) for s in samples  ]  )
                        
    return tree


def do_and_print_biopsy_clustering(G, file_name, find_size, find_type, spatial_plot,
                 hierarchical_plot, biopsy_plot, save_biopsy=False):
    """ Input:  Graph G, file_name,
        find_size is a number, and find_type = "relative" or "absolute"
            refers to the cutoff in the ancestral graph
        takes a biopsy, does a clonal clustering, compares them,
        spatial_plot and hierarchical_plot can be "on" or "off" """
    
    #if type(G) == small_graph and spatial_plot == "on":
    #    raise ValueError("spatial plot not possible with saved graphs")
    
    # take biopsy    
    G, samples = take_circular_biopsy(G=G)
    
    # show 3D-structure of biopsy
    if biopsy_plot == "on":
        show_biopsy(G, samples, file_name)
    
    # creates a tree, clusters and plots the 3D-structure of the clusters
    if spatial_plot == "on":
        tree = EVAL_type_clustering_print.cluster_plot(G,\
                            file_name, find_size, find_type, spatial_plot,
                            "off")
    elif spatial_plot == "off":
        tree =  EVAL_gen_tree_graph.create_gen_graph(G, find_size, find_type)
        tree.branch_based_cluster()
        tree.set_colours()
    else:
        raise ValueError("Do and print biopsy: Specify what you want please!")
    
    tree = compare_biopsy_and_clustering(G, samples, tree)
    
    # clonal hierarchy-plot, shows goodness of Biopsy
    if hierarchical_plot == "on":
        tree.cluster_and_biopsy_print(file_name)
        
    if save_biopsy:
        print("ERROR: do_and_print_biopsy_clustering not connected to biopsy storage")
        
    return G, tree

def take_punching_biopsy(G):
    """ takes a graph and does a biopsy
    
        size of the biopsy is adjusted to the size of the graph
        in the case of the investigated graphs of size |V| = 20mln,
        this procedure results in the claimed 1:8 relation
        
        returns 5 samples as lists of vertices"""
    
    # FIRST STEP: understand size of the graph and adapt parameters

    # assumption: ~ roughly a ball
    diam_0 = 50     # diameter of each punch
    length_0 = 400  # length of each punch
    
    xmax = max ( v.position[0] for v in G.V[1:] )
    xmin = min ( v.position[0] for v in G.V[1:] )
            
    ymax = max ( v.position[1] for v in G.V[1:] )
    ymin = min ( v.position[1] for v in G.V[1:] )
    
    zmax = max ( v.position[2] for v in G.V[1:] )
    zmin = min ( v.position[2] for v in G.V[1:] )
    
    R = min (xmax, -xmin, ymax, -ymin, zmax, -zmin)
    r = R * 0.5
    possible_length = math.sqrt( R **2 - r**2) * 1.8
    
    # define length and radius of sample
    if possible_length < length_0:
        length_sample = possible_length
        print("WARNING: Biopsy used smaller length: "\
                      + str( int(length_sample) ) )
    else:
        length_sample = length_0
        
    possible_rad = r * 0.4
    if possible_rad < diam_0 * 0.5:
        rad_sample = possible_rad
        print("WARNING: Biopsy used smaller diameter: "\
                      + str( int(rad_sample*2) ) )
    else:
        rad_sample = diam_0 * 0.5
    
    # define weighted center of the tumour
    zero = ( (xmax+xmin)*0.5, (ymax+ymin)*0.5, (zmax+zmin)*0.5 )
    
    # SECOND STEP: choose five deterministic (x,y) entry points

    # and random z-coordinates "heights"
    
    points = [ ( zero[0], zero[1]   ),
               ( 0.7  * r + zero[0], -0.7 * r + zero[1] ),
               ( -0.7 * r + zero[0], -0.7 * r + zero[1]),
               ( -0.7 * r + zero[0],  0.7 * r + zero[1]),
               (  0.7 * r + zero[0],  0.7 * r + zero[1])
            ]
    # correspond to mid-point and 4 tubulas around the midpoint (anti-clockwise)
    
    # maximal squared distance in x-y-direction
    check = rad_sample**2
    
    heights = []
    for point in points:
        tempmin = min ( v.position[2] for v in G.V[1:]
                if plain_dist2( v, point ) <= check )
        tempmax = max ( v.position[2] for v in G.V[1:]
                if plain_dist2( v, point ) <= check )
        
        # length of random interval for choosing entry point in z
        rand = tempmax - tempmin - length_sample
        
        if rand > 0:
            rand = random.uniform(0, rand)
            height = rand + tempmin
        else:
            height = tempmin
            
        heights.append(height)
    
    
    # THIRD STEP: take five samples
    samples = []
    
    for i in range(5):
        point = points[i]
        height_low = heights[i]
        height_up = height_low + length_sample
        
        sample = set (   v for v in G.V[1:]
                            if plain_dist2( v, point) <= check
                            and height_low <= v.position[2] <= height_up
                      )
        samples.append(sample)
    
    # updates the vertices of the graph with biopsy-info
    for i in range(len(samples)):
        sam = samples[i]
        for v in sam:
            v.in_biopsy = i + 1
            
    for v in G.V[1:]:
        if v.in_biopsy == None:
                v.in_biopsy = False
    
    return G, samples

  

#%%
def absoluteFilePaths(directory):
    """ Walks through a directory """
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            yield os.path.abspath(os.path.join(dirpath)), os.path.abspath(os.path.join(dirpath, f))

#%%
def get_all_mut(directory, MioSize = 20):
    """ Walks through a directory """
    for direct, file in absoluteFilePaths(directory):
        if direct.endswith( str(MioSize) + "Mio" ) and file.endswith("all_mutations.xlsx"):
            all_mut = pd.read_excel(file)
            break
    return all_mut

#%%
def get_all_mut_biopsies(directory, MioSize = 20):
    """ Loads biopsies from the correct subfolders """
    for direct, file in absoluteFilePaths(directory):
        if direct.endswith( str(MioSize) + "Mio" ):
            D = direct
            break
    # go to right subfolder of biopsies
    for direct, file in absoluteFilePaths(D):
        if direct.endswith("biopsies_punched"):
            D = direct
            break
    biopsies = []
    # get biopsies-all-mutations
    for direct, file in absoluteFilePaths(D):
        if direct[-8:-1] == "biopsy_" and file.endswith("all_mutations.xlsx"):
            table = pd.read_excel(file)
            biopsies.append(table)
    return biopsies
    
#%%
def get_found_in_biopsy_values(trait: str, biopsies):
    """ Walks over all biopsies and finds if a trait was found """
    # convert str to tuple
    trait = eval(trait)
    
    # default: trait not found
    value = [0 for i in range(len(biopsies)) ]
    T = len(trait)
    
    # go through biopsies
    for k in range(len(biopsies)):
        b = biopsies[k]
        gen_form = b["Genetic Formula"]
        for i in range( len(gen_form) ):
            # convert to tuple and compare
            formula = gen_form[i]
            # if None: do nothing
            # else: if type or subtype is found, stop (ok for top-down)
            if type (formula) == str:
                formula = eval(gen_form[i])
                if T <= len(formula):
                    if formula[:T] == trait:
                        value[k] = int(b["Cells"][i])
                        break
    
    # if not found: return "Not found "
    if value == [0 for i in range(len(biopsies)) ]:
        value = "Not found"
    return value

#%%
def find_all_mutations_in_biopsies(graphdir, cutoff_perc = 0.01, MioSize = 20):
    # load all_mutations file from biopsies
    biopsies = get_all_mut_biopsies(graphdir, MioSize)
    
    if len(biopsies) != 5:
        raise ValueError("Could not load all biopsies, found only " +str( len(biopsies)) )
    
    # load information about all_mutations in tumor from excel-sheet
    M = get_all_mut(graphdir, MioSize)
    # Get indices of != None rows
    K = len(M["Index"])
    valid_rows = [ i for i in range(K) if not np.isnan(M["Index"][i])]
    big_rows = [ i for i in valid_rows if M["Percent"][i] >= cutoff_perc ]
    
    # check if the mutation is found in the the biopsies
    found_in_biopsy = [ None for i in range(K)]
    
    for i in big_rows:
        trait = M["Genetic Formula"][i]
        entry = get_found_in_biopsy_values(trait, biopsies)
        found_in_biopsy[i] = entry
    
    # add new information to dataframe
    M["Found in Biopsy"] = found_in_biopsy
    
    # replace the old all_mut with the updated version
    filename = graphdir + "/" + str(MioSize) + "Mio/all_mutations.xlsx"
    M.to_excel(filename, index=False)
    print("DONE WITH", filename, "\n")   
    
#%%
def important_mutations_in_biopsies(graphdir, MioSize = 20,
                                    cutoff_perc = 1.0,
                                    mode = "all_mutations",
                                    biopsy_sens_perc = 5.0):
    """ IN: full directory and graph-size
            mode = "all_mutations" or "clustering"
            biopsy sensitivity in percent
        OUT: biopsy_clonal table in same directory/Mio:
                information about found mutations / genotypes """
    
    # load info about biopsies
    biopsies = get_all_mut_biopsies(graphdir, MioSize)
    if len(biopsies) != 5:
        raise ValueError("Could not load all biopsies, found only " +str( len(biopsies)) )
    biopsy_sizes = [ b["Cells"][0] for b in biopsies ]
    
    # load information about all_mutations in tumor from excel-sheet
    ALL_mut = get_all_mut(graphdir, MioSize)
    
    # load gen-tree for spatial information
    tree = EVAL_data_output.get_tree_sql(graphdir + "/" + str(MioSize) + "Mio",
                                         fullpathname=True)
    
    # Define importand rows (deepcopy of sub-frame)
    if mode == "all_mutations":
        # take all mut above threshold
        BIG_mut = ALL_mut[ ALL_mut["Percent"]>= cutoff_perc ]
        BIG_mut = BIG_mut[ BIG_mut["Genetic Formula"] != "(0,)" ].copy()
        # use branch-based clustering as defined in gen_tree_graph

        tree.branch_based_cluster()
        groups = set( v.name for v in tree.group_points)
        
        def key_in_groups(v):
            if type(v) == str:
                key = eval(v)
                ret = key in groups
                return ret
            else:
                return False
        
        new_keys = set ( v for v in ALL_mut["Genetic Formula"] if key_in_groups(v) )
        criterion = ALL_mut["Genetic Formula"].map(lambda x: x in new_keys)
        BIG_mut = ALL_mut[ criterion ].copy()
    
    def binfound(x):
        if x >= biopsy_sens_perc:
            return "A"
        elif 0 < x:
            return "a"
        else:
            return 0
    
    
    # get info about found/not found in biopsy
    found_perc_info = []
    found_binary_info = []
    for found_str in BIG_mut["Found in Biopsy"]:
        if found_str == "Not found":
            found_perc_info.append("")
            found_binary_info.append("")
        else:
            # convert string to list
            found_vector = eval( found_str )
            if len(found_vector) == len(biopsies):
                # get percent of biopsy-values
                perc_found_vector = [round(100*found_vector[i] / biopsy_sizes[i],1) \
                                             for i in range(len(biopsies))]
                found_perc_info.append(perc_found_vector)
                found_bin = [binfound(x) for x in perc_found_vector]
                found_binary_info.append(found_bin)
            else:
                raise ValueError("Number of biopsies does not match")
      
    found_spatial_info = []
    for mut_str in BIG_mut["Genetic Formula"]:
        mut_name = eval(mut_str)
        v = tree[mut_name]
        vectr = v.spatial_center
        vectr = tuple( [ round(i,2) for i in vectr ] )
        found_spatial_info.append(vectr)

    BIG_mut["Biopsy percent size"] = found_perc_info
    BIG_mut["Biopsy threshold "+str(int(biopsy_sens_perc)) + "%"] = found_binary_info
    BIG_mut["Spatial center"] = found_spatial_info
    
    out = BIG_mut [ ["Percent", "Genetic Formula", "Found in Biopsy",\
                     "Biopsy percent size", "Biopsy threshold "+str(int(biopsy_sens_perc)) + "%",\
                         "Genetic Formula", "Mean fitness", "Maximal fitness", "Spatial center", "Fittest types"] ]

    # store new dataframe
    if mode == "all_mutations":
        filename = graphdir + "/" + str(MioSize) + "Mio/biopsy_clonal_all.xlsx"
    if mode == "clustering":
        filename = graphdir + "/" + str(MioSize) + "Mio/biopsy_clonal_cluster.xlsx"
        
    out.to_excel(filename, index=False)
    print("DONE WITH", filename, "\n")   
    return BIG_mut

def all_modes_important_mutations_in_biopsies(graphdir, MioSize):
    all_mut = important_mutations_in_biopsies(graphdir, cutoff_perc = 1.0,
                                    biopsy_sens_perc = 5.0, MioSize = MioSize,
                                    mode = "all_mutations")
    cluster = important_mutations_in_biopsies(graphdir, cutoff_perc = 1.0,
                                    biopsy_sens_perc = 5.0, MioSize = MioSize,
                                    mode = "clustering")
    return
    
#%%

def plain_dist2( v, point):
    """ asymmetrical input: vertex v and coordinate-point 
        distance in (x,y)-plane """
    return ( v.position[0] - point[0] )**2 + ( v.position[1] - point[1] )**2

def dist2 ( v, point):
    """ asymmetrical input: vertex v and coordinate-point """
    return  ( v.position[0] - point[0] )**2 +\
            ( v.position[1] - point[1] )**2 +\
            ( v.position[2] - point[2] )**2
            
def dist_height (v, z_level):
    return ( abs(v.position[2] - z_level) )

def sup_dist_in_plain(v, point):
    a = abs ( v.position[0] - point[0]  )
    b = abs ( v.position[1] - point[1]  )
    return max(a,b)

def update_graph_with_biopsy_info(G, total_samples):
    #Update the vertices of the graph with biopsy-info
    
    # delete old info (not needed, but good for safety)
    for v in G.V[1:]:
        v.in_biopsy = None
    
    # add new info
    for i in range(len(total_samples)):
        sam = total_samples[i]
        for v in sam:
            v.in_biopsy = i + 1
            
    for v in G.V[1:]:
        if v.in_biopsy == None:
                v.in_biopsy = False
  
    return G