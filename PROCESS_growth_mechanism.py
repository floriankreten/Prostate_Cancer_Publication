#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 11:08:03 2018

@author: floriankreten
"""

""" Geometric calculations of a growth- or branch-step.
    new_growth_direction and branch_directions are called from 
        INPUT_rates at a growth or branch-event,
        they only return the new random dircetions
"""

import math
import numpy.random
import PROCESS_structures

def random_radius(noise_strength):
    """ IN: noise_strength in degrees
        OUT: radius according to random angle chosen
             uniformly from (-noise_strength, noise_strength)"""
            
    # convert angle in degrees to angle in radians
    noise_strength = math.radians(noise_strength)
    # choose a random angle within the range of the noise strength
    noise = numpy.random.uniform(-noise_strength, noise_strength)
    
    radius = math.sqrt ( ( 1 /   (( math.cos(noise))**2)   ) - 1 )
    return radius

def normalize (a):
    """ normalizes a vector in R^d """
    norm = sum ( x**2 for x in a)
    try:
        norm = 1.0 / math.sqrt(norm)
    except ZeroDivisionError:
        raise ZeroDivisionError ("Normalization failed, norm == 0")
    
    output = [ x * norm for x in a ]
    return tuple (output)

def first_orth_vector(x):
    """ IN: vector in R^3 which is not equal to zero
        OUT: normalized orthogonal vector"""
    if abs ( x[0] ) + abs ( x[1] ) > 0:
        u = (-x[1], x[0] , 0  )
    else:
        u = ( -x[2], 0, x[0] )
    return (normalize(u))

def second_orth_vector(a,b):
    """ in: vector in R^3
        returns a normalized vector orthogonal to the two others
        up to factor 100 faster than library methods"""
    # cross-product
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
    # normalize
    return normalize(c)

def new_growth_direction(x,noise_strength):
    """ input: noise_strength in degrees and given direction x
        output: new random-direction"""
    beta = numpy.random.uniform ( 0, 2*math.pi )
    
    #equations for the plane
    u = first_orth_vector(x)
    v = second_orth_vector(x,u)
    r = random_radius(noise_strength)
    
    #new growth direction
    r_vector = [ r * ( math.cos(beta) * u[i] + math.sin(beta) * v[i] )
                for i in (0,1,2) ]
    output = PROCESS_structures.PLUS ( x, r_vector)
    
    return normalize(output)


def branch_directions(x,angle):
    """ input: given direction x, fixed branching  angle in degrees
        returns: two new growth directions"""
    alpha = math.radians (angle)
    #equations for the plane
    beta = numpy.random.uniform ( 0, 2*math.pi )
    u = first_orth_vector(x)
    v = second_orth_vector(x,u)
    r = math.sqrt ( ( 1 /  (  (math.cos(alpha))**2) ) - 1 )
    
    #calculates first growth direction
    r_vector = [ r * ( math.cos(beta) * u[i] + math.sin(beta) * v[i] )
                for i in (0,1,2) ]
    out_1 = PROCESS_structures.PLUS ( x, r_vector)
    out_1 = normalize (out_1)
    
    #calculates second growth direction
    out_2 = PROCESS_structures.MINUS ( x, r_vector)
    out_2 = normalize (out_2)
    
    return (out_1,out_2)
