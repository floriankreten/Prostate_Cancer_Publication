#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 14:16:51 2018

@author: floriankreten
"""

""" A container is used in the graph_class (see PROCESS_structures)

    Sorts the vertices of a graph into a 3d-grid
    This is needed to check if a new point in space is free or blocked 
    The passed "whitelist" is defined in the growth / radial_growth modules
        depending on the growth-event, certain neighboring vertices
        should be ignored as "blocking"- vertices"

"""

import PROCESS_structures
import math


class cell_type_container:
    """ needed to check if a point in space is blocked
        Input:
            size (of simulation) = maximal inf_norm of a coordinate of a vertex
            container_length = length of a box / block of the grid
        Warning: is explicitely written for dimension three
    """
    def __init__(self, container_length, size):
        if size <= container_length:
            raise ValueError ("Container input: size <= absorption")
        if container_length < 1.5:
            raise Warning ("Container input: Container Size too small")
        if container_length > 20:
            raise Warning("Container input: absorption too high")
        if size < 1:
            raise ValueError ("Container input: size invalid")
        size = int(size)
        self.size = size
        self.L = max ( container_length, 4 ) # size of one grid-point-container
        self.D = int ( (self.size/self.L + 1) ) #grid-points are from (-D,D)
        self.container = [ [  [ [] for x in range(2*self.D)]
                    for y in range(2*self.D) ] for z in range(2*self.D) ]
    
    def __setitem__(self, xyz, value):
        """ values in the container are shifted by D
            this allows negative entries up to -D
        """
        x,y,z = xyz
        for i in xyz:
            if abs (i) >= self.size:
                raise ValueError ("Container: Coordinate of new entry exceeds size")
        self.container[x+self.D][y+self.D][z+self.D] = value
    
    def __getitem__(self, xyz):
        x,y,z = xyz
        for i in xyz:
            if abs (i) >= self.size:
                raise ValueError ("Container: Call exceeds size")
        try:
            v = self.container[x+self.D][y+self.D][z+self.D]
        except IndexError:
            print("Call exceeds container size")
            raise IndexError
        
        return v
    
    def insert(self,G,key): ### inserts a vertex into a grid-point
        a = G[key].position
        x = math.floor ( a[0]/self.L )
        y = math.floor ( a[1]/self.L )
        z = math.floor ( a[2]/self.L )
        self[x,y,z].append (key)
        
    def __str__(self):
        output = "Container with grid-size " + str(self.L)
        for x in range(-self.D,self.D):
            for y in range(-self.D,self.D):
                for z in range(-self.D,self.D):
                    if self[x,y,z] != []:
                        output += "\n" + str( (self.L*x,self.L*y,self.L*z) ) + " ->  "
                        for KEY in self[x,y,z]:
                            output += " | " + str(KEY)
        return output
    
    def point_procedure(self, G, key, pos, new_size, whitelist):
        """ Checks if a given point "pos" touches a vertex with index "key"
            within the absorption radius delta
            
            points who lie in on same or close sibling branch (=whitelist)
            dont matter for collisions (up to a certain search-depth)
            T =  max ( 2, delta + 1   ), see growth-modules for definition """
            
        v = G.V[key]
        old_pos = v.position
        old_size = v.radial_size
        deltapow2 = (old_size+new_size) **2
        if PROCESS_structures.ABS2 (old_pos,pos) < deltapow2:
            if key not in whitelist:
                return False
        else:
            return True
        
    def check_box_is_free(self,G,pos,new_size,whitelist,x,y,z):
        for key in self[x,y,z]:
            var = self.point_procedure(G, key, pos, new_size, whitelist)
            if var == False:
                return False
    
    def is_free(self, G, pos, new_size, father_key, whitelist):
        """ checks, if point pos=(x,y,z) is blocked by an other vertex
            returns True if space is free
            if pos=False, the current node is trying to expand (pos=father_key.pos)
        """
        
        if new_size > 3:
            print(G[father_key])
            raise TypeError("Radial size shall not exceede 3 for geometrical reasons!")
        
        
        if pos == None:
            pos = G.V[father_key].position
            raise ValueError("YOU SHALL NOT PASS")
        
        #get to a corresponding grid point
        for i in pos:
            if abs(i) >= self.size:
                raise ValueError ("Container error: Graph exceeds coordinate size")
        x = math.floor ( pos[0]/self.L )
        y = math.floor ( pos[1]/self.L )
        z = math.floor ( pos[2]/self.L )
        
        x0 = pos[0] - x
        y0 = pos[1] - y
        z0 = pos[2] - z
        
        #checks the corredsponding grid-container
        var = self.check_box_is_free(G,pos,new_size,whitelist,x,y,z)
        if var == False:
            return False
        
        "Check neighboring sites"
        upper = self.L - z0 < new_size  and  z < self.D 
        lower = z0 < new_size           and self.D > -z
        right = self.L - y0 < new_size  and  y < self.D 
        left  = y0 < new_size           and self.D > -y
        behind= self.L - x0 < new_size  and  x < self.D 
        front = x0 < new_size           and self.D > -x
        
        neighb_tests = [ behind, right, upper,
                        
                        front, left, lower,
                        
            behind and right, behind and upper, right and upper,
            
            front and left, front and lower, left and lower,
                        
            front and right, front and upper, left and upper,
            
            right and lower, behind and lower, behind and left,
            
            behind and right and upper, front and left and lower,
            
            front and left and upper, front and right and lower, behind and left and lower,
            
            behind and right and lower, behind and left and upper, front and right and upper
            
            
            ]
        
        neighb_boxes  = neighboring_boxes(x,y,z)
            
        for i in range(len(neighb_tests)):
            if neighb_tests[i]:
                x,y,z = neighb_boxes[i]
                test = self.check_box_is_free(G,pos,new_size,whitelist,x,y,z)
                if test == False:
                    return False

        return True
    
    
def neighboring_boxes(x,y,z):
    """ The 9 neighboring boxes of a box with coordinates x,y,z """
    neighb_boxes        = [
                        [x+1,y,z],
                        [x,y+1,z],
                        [x,y,z+1],
                        
                        [x-1, y, z],
                        [x, y-1, z],
                        [x, y, z-1],
                        
                        [x+1, y+1, z],
                        [x+1, y, z+1],
                        [x, y+1, z+1],
                        
                        [x-1, y-1, z],
                        [x-1, y, z-1],
                        [x, y-1, z-1],
                        
                        [x-1, y+1, z],
                        [x-1, y, z+1],
                        [x, y-1, z+1],
                        
                        [x, y+1, z-1],
                        [x+1, y, z-1],
                        [x+1, y-1, z],
                        
                        [x+1, y+1, z+1],
                        [x-1, y-1, z-1],
                        
                        [x-1, y-1, z+1],
                        [x-1, y+1, z-1],
                        [x+1, y-1, z-1],
                        
                        [x+1, y+1, z-1],
                        [x+1, y-1, z+1],
                        [x-1, y+1, z+1]     ]
    return neighb_boxes
    