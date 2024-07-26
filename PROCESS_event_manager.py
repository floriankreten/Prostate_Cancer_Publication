#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 15:41:35 2019

@author: floriankreten
"""

""" An event-manager is part of a graph
    It keeps track of all possible events and corresponding rates """

""" Each event can be executed in O(1)
    This can be done since the different types of events are stored in lists.
    Their ordering does not matter, so deletions can be done in O(1),
        because deletions can be shifted to the end of the list.
    After any event, the rates need to be updated only in the direct
    neighborhood of a vertex G[x], making updates rather easy and fast.
    """

class EVENT_manager:
    """ manages the rates of the different random events
        events are pre-sorted: growth, branch, mutation, competition
        sampling of a Gillespie-step is done in three steps:
            symple type of event -> symple vertex
                -> execute event (e.g. choose an edge)
    """
    
    def __init__ (self, G):
        
        # the rates of G have to be up to date BEFORE this is called
        # initialize them first, using the optional "ignore"-parameter
        
        self.total_rate_mutation = sum ( v.rate_mutation for v in G.V[1:] )
        self.total_rate_growth = sum ( v.rate_growth for v in G.V[1:]  )
        self.total_rate_branch = sum ( v.rate_branch for v in G.V[1:] )
        self.total_rate_competition = sum (v.rate_competition for v in G.V[1:])
        self.total_rate_radial_growth = sum ( v.rate_radial_growth for v in G.V[1:] )
        
        self.graph_rate =     self.total_rate_branch \
                            + self.total_rate_growth \
                            + self.total_rate_competition \
                            + self.total_rate_mutation \
                            + self.total_rate_radial_growth
        
        self.step_list_mutation = []
        self.step_list_growth = []
        self.step_list_branch = []
        self.step_list_competition = []
        self.step_list_radial_growth = []
        
        self.num_mutating = 0
        self.num_growing = 0
        self.num_branching = 0
        self.num_competing = 0
        self.num_radial_growing = 0
        
        # creates lists of active vertices (with rates > 0)
        for v in G.V[1:]:
            
            if v.status_mutating:
                self.step_list_mutation.append(v)
                v.index_mutation = self.num_mutating
                self.num_mutating += 1
                
            if v.status_growing:
                self.step_list_growth.append(v)
                v.index_growth = self.num_growing
                self.num_growing +=1
                
                
            if v.status_branching:
                self.step_list_branch.append(v)
                v.index_branch = self.num_branching
                self.num_branching +=1
                
                
            if v.status_competing:
                self.step_list_competition.append(v)
                v.index_competition = self.num_competing
                self.num_competing += 1
                
            if v.status_radial_growing:
                self.step_list_radial_growth.append(v)
                v.index_radial_growth = self.num_radial_growing
                self.num_radial_growing += 1
        #

        self.max_rate_growth =      max ( v.rate_growth for v in G.V[1:]  )
        self.max_rate_branch =      max ( v.rate_branch for v in G.V[1:]  )
        self.max_rate_competition = max ( v.rate_competition for v in G.V[1:]  )
        self.max_rate_mutation =    max ( v.rate_mutation for v in G.V[1:]  )
        self.max_rate_radial_growth = max ( v.rate_radial_growth for v in G.V[1:] )
        
    def __str__(self):
            # Print_function for graph_manager
        output = "\n__EVENT MANAGER__"
        output += "\nRate of mutation:" + str(self.total_rate_mutation)
        output += "\nRate of growth:" + str(self.total_rate_growth)
        output += "\nRate branch:" + str(self.total_rate_branch)
        output += "\nRate competition:" + str(self.total_rate_competition)
        
        output += "\nTotal rate:" + str(self.graph_rate) + "\n"
        
        return output
        
        
def UPDATE_EVENT_manager_growth(G,x, diff_growth, diff_branch):
    """ updates the event manager after a growth-event at a vertex G[x] """
    
    # updates branch-rates first
    G = UPDATE_EVENT_manager_branch(G,x,diff_branch)
    
    # update EVENT-manager
    E = G.EVENT_manager
    E.graph_rate = G.graph_rate
    
    # 1) update growth-rates
    # dumb case: vertex just needs to be ignored
    if G[x].status_growing == False and G[x].index_growth == None:
        return G
    
    # first case: vertex gets deactivated after he was avtive before
    if G[x].status_growing == False and G[x].index_growth != None:
        # update list of growing vertices
        change_key = G[x].index_growth
        E.step_list_growth[change_key] =  E.step_list_growth[-1]
        E.step_list_growth[change_key].index_growth = change_key
        E.num_growing -= 1
        G[x].index_growth = None
        del (E.step_list_growth[-1])
        # update rates
        E.total_rate_growth += diff_growth
        E.total_rate_growth = max(E.total_rate_growth, 0.0)
        
    # second case: vertex is active and changes
    elif G[x].status_growing == True and G[x].index_growth != None:
        # update rates
        E.total_rate_growth += diff_growth
        E.total_rate_growth = max(E.total_rate_growth, 0.0)
        E.max_rate_growth = max ( [E.max_rate_growth, G[x].rate_growth] )
        
    # third case: vertex gets activated
    elif G[x].status_growing == True and G[x].index_growth == None:
        # update list
        E.step_list_growth.append(G[x])
        G[x].index_growth = E.num_growing
        E.num_growing += 1
        # update rates
        E.total_rate_growth += diff_growth
        E.max_rate_growth = max ( [E.max_rate_growth, G[x].rate_growth] )
        
    return G
    

def UPDATE_EVENT_manager_branch(G,x,diff_branch):
    """ updates the event manager after a branch-event at a vertex G[x] """
    # update EVENT-manager
    E = G.EVENT_manager
    E.graph_rate = G.graph_rate
    
    # 1) update branch-rates
    # dumb case: vertex just needs to be ignored
    if G[x].status_branching == False and G[x].index_branch == None:
        return G
    
    # first case: vertex gets deactivated after he was avtive before
    if G[x].status_branching == False and G[x].index_branch != None:
        # update list of growing vertices
        change_key = G[x].index_branch
        E.step_list_branch[change_key] =  E.step_list_branch[-1]
        E.step_list_branch[change_key].index_branch = change_key
        E.num_branching -= 1
        G[x].index_branch = None
        del (E.step_list_branch[-1])
        # update rates
        E.total_rate_branch += diff_branch
        E.total_rate_branch = max (  E.total_rate_branch, 0.0 )
        
    # second case: vertex is active and changes
    elif G[x].status_branching == True and G[x].index_branch != None:
        # update rates
        E.total_rate_branch += diff_branch
        E.total_rate_branch = max (  E.total_rate_branch, 0.0 )
        E.max_rate_branch = max ( [E.max_rate_branch, G[x].rate_branch] )
        
    # third case: vertex gets activated
    elif G[x].status_branching == True and G[x].index_branch == None:
        # update list
        E.step_list_branch.append(G[x])
        G[x].index_branch = E.num_branching
        E.num_branching += 1
        # update rates
        E.total_rate_branch += diff_branch
        E.max_rate_branch = max ( [E.max_rate_branch, G[x].rate_branch] )
        
    return G
    

def UPDATE_EVENT_manager_mutation(G,x,difference):
    """ updates the event manager after a mutation-event at a vertex G[x] """
    # update Event-manager
    E = G.EVENT_manager
    E.graph_rate = G.graph_rate
    
    # dumb case: vertex is ignored
    if G[x].status_mutating == False and G[x].index_mutation == None:
        return
    
    # first case: vertex is active and changes
    if G[x].status_mutating == True and G[x].index_mutation != None:
        # update rates
        E.total_rate_mutation += difference
        E.max_rate_mutation = max ( [E.max_rate_mutation, G[x].rate_mutation] )
        
    # second case: vertex gets activated
    elif G[x].status_mutating == True and G[x].index_mutation == None:
        # update list
        E.step_list_mutation.append(G[x])
        G[x].index_mutation = E.num_mutating
        E.num_mutating += 1
        # update rates
        E.total_rate_mutation += difference
        E.max_rate_mutation = max ( [E.max_rate_mutation, G[x].rate_mutation] )
    
    # third case: vertex gets deactivated
    elif G[x].status_mutating == False and G[x].index_mutation != None:
        # update list
        change_key = G[x].index_mutation
        E.step_list_mutation[change_key] =  E.step_list_mutation[-1]
        E.step_list_mutation[change_key].index_mutation = change_key
        E.num_mutating -= 1
        G[x].index_mutation = None
        del (E.step_list_mutation[-1])
        # update rates
        E.total_rate_mutation += difference
        
    return G
    

def UPDATE_EVENT_manager_competition(G,x,diff_total):
    """ updates the event manager after a competition-event at a vertex G[x] """
    # update Event-manager                    
    E = G.EVENT_manager
    E.graph_rate = G.graph_rate
    
    # dumb case: ignore the vertex, was not competing and still isnt
    if G[x].status_competing == False and G[x].index_competition == None:
        return
    
    #first case: vertex gets activated
    elif G[x].status_competing == True and G[x].index_competition == None:
        # update list
        E.step_list_competition.append(G[x])
        G[x].index_competition = E.num_competing
        E.num_competing += 1
        # update rates
        E.total_rate_competition += diff_total
        E.max_rate_competition = max (
                        [E.max_rate_competition, G[x].rate_competition] )
        
    # second case: vertex gets deactivated
    elif G[x].status_competing == False and G[x].index_competition != None:
        # update list
        change_key = G[x].index_competition
        E.step_list_competition[change_key] =  E.step_list_competition[-1]
        E.step_list_competition[change_key].index_competition = change_key
        E.num_competing -= 1
        G[x].index_competition = None
        # update rates
        E.total_rate_competition += diff_total
        E.total_rate_competition = max ( 0.0, E.total_rate_competition )
        del (E.step_list_competition[-1])
        
    # third case: vertex changes while being active (by mutation)
    elif G[x].status_competing == True and G[x].index_competition != None:
        # update rates
        E.total_rate_competition += diff_total
        E.max_rate_competition = max (
                        [E.max_rate_competition, G[x].rate_competition] )
        
    return G

def UPDATE_EVENT_manager_radial_growth(G,x,difference):
    """ updates the event manager after a radial-growth-event at a vertex G[x] """

    # update event manager
    E = G.EVENT_manager
    E.graph_rate = G.graph_rate
    
    # dumb case: vertex is ignored
    if G[x].status_radial_growing == False and G[x].index_radial_growth == None:
        return
    
    #first case: vertex gets activated
    elif G[x].status_radial_growing == True and G[x].index_radial_growth == None:
        # update list
        E.step_list_radial_growth.append(G[x])
        G[x].index_radial_growth = E.num_radial_growing
        E.num_radial_growing += 1
        # update rates
        E.total_rate_radial_growth += difference
        E.max_rate_radial_growth = max (
                        [E.max_rate_radial_growth, G[x].rate_radial_growth] )
        
    # second case: vertex gets deactivated
    elif G[x].status_radial_growing == False and G[x].index_radial_growth != None:
        # update list
        change_key = G[x].index_radial_growth
        E.step_list_radial_growth[change_key] =  E.step_list_radial_growth[-1]
        E.step_list_radial_growth[change_key].index_radial_growth = change_key
        E.num_radial_growing -= 1
        if E.num_radial_growing < 0:
            print(G[x])
            raise ValueError("Num Radial < 0 !")
        G[x].index_radial_growth = None
        # update rates
        E.total_rate_radial_growth += difference
        E.total_rate_radial_growth = max ( 0.0, E.total_rate_radial_growth )
        del ( E.step_list_radial_growth[-1] )
        
    # third case: vertex changes while being active (i.e. after growth-event)
    elif G[x].status_radial_growing == True and G[x].index_radial_growth != None:
        # update rates
        E.total_rate_competition += difference
        E.max_rate_radial_growth = max (
                        [E.max_rate_radial_growth, G[x].rate_radial_growth] )
        
    return G    

    
    return