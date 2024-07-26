#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:45:40 2019

@author: floriankreten
"""

""" Two sorts of objects can be stored in databases: ancestral trees and graphs
    for reading in from a distant source (not /output), pass 
   "fullpathname"=True in **kwargs """

""" WARNING: Radial Size of a vertex v is not used at the moment,
    v.size always taken from INPUT_rates """

import INPUT_rates

import os
import json
import math
import gc
import sqlite3
from pandas import DataFrame
import pandas as pd


# for rebuilding the graph
import PROCESS_structures
import PROCESS_mutation
import PROCESS_process
import PROCESS_radial_growth

import PROCESS_growth
import PROCESS_competition
import EVAL_gen_tree_graph

from pathlib import Path
#Path("/my/directory").mkdir(parents=True, exist_ok=True)


#%%
# Preliminaries
##################################
def create_folder(file_name, *file_type):
    """ Creates directory for storage: output/file_name/file_type
        or output/file_name if no second argument is passed """
    directory = os.getcwd()
    dirName = str ( directory + "/output/" + file_name )
    
    Path(dirName).mkdir(parents=True, exist_ok=True)
    
    # create an additional subfolder if this argument is passed
    if file_type:
        file_type = file_type[0]
        dirName = str ( directory + "/output/" +  file_name + "/" + file_type)
        Path(dirName).mkdir(parents=True, exist_ok=True)
    
    return (dirName)
        

def absoluteFilePaths(directory):
    """ Walks through files in directory in generator-format """
    for dirpath,_,filenames in os.walk(directory):
       for f in filenames:
           yield os.path.abspath(os.path.join(dirpath, f))

#%%
# Store and extract functions for graphs
########################################
    
def save_graph_sql(G, filename, fullpathname = False):
    """ 
        1) reduces the size of the gen_tree (for 2 dicts used in the graph)
        2) stores the vertex-information in SQL-database
        3) store the graph.gen_tree and graph.fitness_book in the database
        4) Optional i makes is for numbering the saves of the graph within
                the same folder
        WARNING: DELETES ALL OLD ENTRIES FIRST
    """
    # reduce the gen-tree up to the MRCA
    # G = PROCESS_process.reduce_gen_tree(G)
    
    # Creates directory for storage
    if not fullpathname:
        dirName = create_folder(filename)
    else:
        dirName = filename
        Path(dirName).mkdir(parents=True, exist_ok=True)
    
    # prepare the database
    name = str ( dirName + "/" + "graph.db" )
    
    try:
        os.remove(name)
    except:
        pass
    
    database = sqlite3.connect(name)
    cursor = database.cursor()
    
    # vertex-table
    cursor.execute( """ DROP TABLE if exists vertices """)
    create_vertex_table = """   CREATE TABLE vertices ( 
                                vertex_key INTEGER,
                                vertex_pos_x FLOAT,
                                vertex_pos_y FLOAT,
                                vertex_pos_z FLOAT,
                                
                                in_edge INTEGER,
                                out_edge_1 INTEGER,
                                out_edge_2 INTEGER,
                                
                                fitness FLOAT, 
                                trait TEXT,
                                
                                radial_size FLOAT,
                                radial_cells INTEGER,
                                status_radial_growing INTEGER,
                                
                                status_growing INTEGER,
                                status_branching INTEGER,
                                cluster_group INTEGER,
                                in_biopsy INTEGER
                                ) """
    cursor.execute(create_vertex_table)
    database.commit()
    
    # how do save a vertex from the graph in the database
    insert_vertex_comm = """ INSERT INTO vertices (
                    vertex_key,
                    vertex_pos_x, vertex_pos_y, vertex_pos_z,
                    in_edge, out_edge_1, out_edge_2, fitness,
                    trait, status_growing, status_branching,
                    cluster_group, in_biopsy, radial_size, radial_cells,
                    status_radial_growing)
                    VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) """
    
    def insert_node(v):
        # transfer Cluster: ["unclassified", None, Int] to [-42, -1, Int]
        if type(v.cluster_group)== int:
            cluster = int(v.cluster_group)
        elif type(v.cluster_group) == str:
            cluster = -42
        elif v.cluster_group == None:
            cluster = -1
        
        if v.in_biopsy != None:
            in_biop= int(v.in_biopsy)
        else:
            in_biop = -1
        
        cursor.execute(insert_vertex_comm, [
                    v.key,
                    v.position[0],
                    v.position[1],
                    v.position[2],
                    v.in_edge,
                    v.out_edge_1,
                    v.out_edge_2,
                    v.fitness,
                    json.dumps(v.trait),
                    int(v.status_growing),
                    int(v.status_branching),
                    cluster,
                    in_biop,
                    v.radial_size,
                    v.radial_cells,
                    v.status_radial_growing
                    ])
    
    # gen-tree-table
    cursor.execute( """ DROP TABLE if exists genotypes """)
    create_genotype_table = """ CREATE TABLE genotypes (
                                gen_key TEXT,
                                free_child TEXT,
                                fitness FLOAT
                                ) """
    
    cursor.execute(create_genotype_table)
    insert_gen_tree_comm = """ INSERT INTO genotypes (
                            gen_key,
                            free_child,
                            fitness )
                        VALUES(?,?,?) """
                        
    # parameters used
    cursor.execute ( """ DROP TABLE if exists information """)
    insert_parameters_comm = """ CREATE TABLE information (
                            parameters TEXT ) """
    cursor.execute(insert_parameters_comm)
    parameters = json.dumps(G.parameters)
    cursor.execute ( """ INSERT INTO information (
                            parameters )
                            VALUES(?) """, [parameters] )
    database.commit()
    
    # history of strong mutations
    cursor.execute ( """ DROP TABLE if exists strong_mutations_log """)
    insert_strong_mutations_log = """ CREATE TABLE strong_mutations_log (
                                        mutation_key TEXT,
                                        mutation_time FLOAT
                                    ) """
    cursor.execute (insert_strong_mutations_log)
    insert_log_comm = """ INSERT INTO strong_mutations_log (
                            mutation_key,
                            mutation_time )
                        VALUES(?,?) """
    for key in G.strong_mutations_log.keys():
        time_of_mut = G.strong_mutations_log[key]
        cursor.execute( insert_log_comm , [ json.dumps(key), float(time_of_mut) ] )
    database.commit()
    
    # how do save genotype information in the database
    def insert_genotype(gen_key):
        name = json.dumps(gen_key)
        if gen_key in G.gen_tree.book:
            free_child = json.dumps(G.gen_tree.book[gen_key])
        else:
            free_child = json.dumps(-1)
        fitness = G.fitness_book[gen_key]
        cursor.execute ( insert_gen_tree_comm,[
                    name, free_child, fitness
                    ])
        
    for v in G.V[1:]:
        insert_node(v)
    database.commit()
    for genotype in G.fitness_book:
        insert_genotype(genotype)
    database.commit()
    database.close()
    
    print("Saved graph at ", name)
    
    return G
    ###################################
    
    
def get_gen_tree_of_graph(filename, **kwargs):
    """ used by get_graph_sql only
        go to /output/filename/graph.db,
        takes specific information from there """

    directory = os.getcwd()
    dirName = str ( directory + "/output/" + str(filename) )
    
    name = str ( dirName + "/" + "graph.db" )
    
    if kwargs is not None:
        if "fullpathname" in kwargs.keys():
            name = filename + "/" + "graph.db"
    
    db = sqlite3.connect(name)
    cursor = db.cursor()
    fitness_book = dict()
    gen_tree_book = dict()
    
    # get values from database
    cursor.execute ( """ SELECT gen_key, free_child, fitness FROM genotypes """ )
    values = cursor.fetchall()
    
    for i in range( len(values) ):
        coll = values[i]
        
        gen_key = tuple ( json.loads(coll[0]) )
        child = tuple ( json.loads (coll[1]) )
        fitness = coll[2]
        fitness_book[gen_key] = fitness
        gen_tree_book[gen_key] = child
    
    db.close()
    del(values)
    gc.collect()
    
    return fitness_book, gen_tree_book

def get_graph_info(filename, fullpathname = False):
    """ used by get_graph_sql only
        Loads graph-infos from a database """
    # go to output/filename/graph.db
    
    if fullpathname:
        dirName = filename
    else:
       # go to output/filename
        directory = os.getcwd()
        dirName = str ( directory + "/output/" + filename )
    
    # find tree        
    for f in absoluteFilePaths(dirName):
        if f.endswith("specifications.xlsx"):
            file_name = f
            break

    import pandas as pd

    df = pd.read_excel(file_name)
    
    return df

def get_graph_sql(filename, container=False, rates=False, *container_size,**kwargs):
    """ Loads a graph from a database
        container/rates can be true/false for initializing those parts of the graph
        container-size is an optional estimate
        i=0 for normal usage 
        use kwarg "fullpathname" for different directory:
                -> leads to "fullpathname/graph.db """
                
    
    
    # prepare the dictionary of genotypes (insert comes later)
    fitness_book, gen_tree_book = get_gen_tree_of_graph(filename, **kwargs)

    # go to output/filename/graph.db
    directory = os.getcwd()
    dirName = str ( directory + "/output/" + str(filename) )
    
    name = str ( dirName + "/" + "graph.db" )

    if kwargs is not None:
        if "fullpathname" in kwargs.keys():
            name = filename + "/" + "graph.db"
    
    # for speed-up: needs more ram! use db_in otherwise
    db_in = sqlite3.connect(name)
    db = sqlite3.connect(':memory:')
    
    db_in.backup(db)
    cursor = db.cursor()
    
    cursor.execute ( """ SELECT gen_key, free_child, fitness FROM genotypes """ )
    values = cursor.fetchall()
    
    for i in range( len(values) ):
        coll = values[i]
        
        gen_key = tuple ( json.loads(coll[0]) )
        child = tuple ( json.loads (coll[1]) )
        fitness = coll[2]
        fitness_book[gen_key] = fitness
        gen_tree_book[gen_key] = child
    
    del(values)
    gc.collect()
    
    
    # get first information about the graph
    cursor.execute ( """ SELECT vertex_pos_x, vertex_pos_y, vertex_pos_z from vertices """)
    positions = cursor.fetchall()
    
    cursor.execute ( """SELECT vertex_key FROM vertices """)
    keys =  [ a[0] for a in cursor.fetchall() ]
    
    cursor.execute (""" SELECT parameters FROM information """)
    parameters = json.loads ( cursor.fetchall()[0][0] )
    
    # initialize graph
    try:
        delta = parameters["container_box_len"]
    except KeyError:
        delta = 10
    
    if container_size:
        container_size = container_size[0]
    else:
        container_size = 0
    size = max( container_size, \
            1.5 * max ( PROCESS_structures.ABS(pos, (0,0,0) ) for pos in positions) )
    size = min(size, 2000)
    size = max(size, 100)
    
    G = PROCESS_structures.graph(delta, size)
    G.parameters = parameters
    try:
        G.time = G.parameters["time"]
    except Exception:
        pass
    
    # place all vertices (with positions only)
    for i in range(len(positions)):
        pos = positions[i]
        G.add_vertex(pos)
        if G.count != keys[i]:
            raise ValueError()
    
    del(positions)
    del(keys)
    gc.collect()
    
    
    # insert the gen-tree
    G.fitness_book = fitness_book
    G.gen_tree = PROCESS_mutation.gen_tree()
    G.gen_tree.book = gen_tree_book
    
    # get all infos about the vertices and place them into the graph-structure
    # adjacent edges
    cursor.execute(""" SELECT in_edge, out_edge_1, out_edge_2 from vertices """)
    edges = cursor.fetchall()
    for i in range( G.count ):
        v = G[i+1]
        coll = edges[i]
        v.in_edge = coll[0]
        v.out_edge_1 = coll[1]
        v.out_edge_2 = coll[2]
        
    del(edges)
    gc.collect()
    
    # trait and fitness
    cursor.execute( """ SELECT trait, fitness from vertices """)
    values = cursor.fetchall()
    for i in range( G.count ):
        v = G[i+1]
        coll = values[i]
        v.trait = tuple( json.loads( coll[0] ) )
        v.fitness = coll[1]
    
    del(values)
    gc.collect()
    
    # growth/branching
    cursor.execute ( """ SELECT status_growing, status_branching
                        FROM vertices """)
    values = cursor.fetchall()
    for i in range( G.count ):
        v = G[i+1]
        coll = values[i]
        v.status_growing = bool(coll[0])
        v.status_branching = bool(coll[1])
        
    del(values)
    gc.collect()

    # radial size and activiy, deactivated atm
    for i in range( G.count ):
            v = G[i+1]
            v.radial_size = INPUT_rates.parameters["vertex_radius"]
    # cursor.execute ( """ SELECT radial_size, radial_cells, status_radial_growing
    #                         FROM vertices """)
    """
    values = cursor.fetchall()
    for i in range( G.count ):
        v = G[i+1]
        coll = values[i]
        v.radial_size = coll[0]
        v.radial_cells = coll[1]
        v.status_radial_growing = bool(coll[2])
            
    del(values)
    gc.collect()
    """
    
    # biopsy and clustering
    cursor.execute ( """ SELECT cluster_group, in_biopsy FROM vertices """ )
    values = cursor.fetchall()
    for i in range (G.count):
        
        v = G[i+1]
        coll = values[i]
        cluster_group = coll[0]
        # [-42, -1, Int] to ["unclassified", None, Int]
        if cluster_group >= 0:
            v.cluster_group = cluster_group
        elif v.cluster_group == -1:
            v.cluster_group = None
        else:
            v.cluster_group = "unclassified"
        
        in_biopsy = coll[1]
        if in_biopsy != -1:
            v.in_biopsy = in_biopsy
    
    del(values)
    gc.collect()
    
    # get strong_mut_log
    cursor.execute ( """ SELECT mutation_key, mutation_time
                        FROM strong_mutations_log """)
    values = cursor.fetchall()
    for str_mut in values:
        key  = json.loads( str_mut[0] )
        time = str_mut[1]
        G.strong_mutations_log[tuple(key)] = time
    
    
    del(values)
    gc.collect()
    
    
    if container or rates:
        G = reinitialize_graph_for_simulations(G, container, rates)
        print("Graph reinitialized")
        
    db_in.close()
    db.close()
    return G


def reinitialize_graph_for_simulations(G, container, rates):
    """ after extracting the graph, some informations for the simulations 
        are still missing: rates, EVENT_manager, container
        those can easily be obtained from the state of the graph
        this function does the necessary steps """
    
    
    # updates the container    
    for key in range( 1, G.count+1 ):
        G.container.insert(G, key)
    
    # initializes all rates and growth directions 
    if rates:
        for v in G.V[1:]:
            v.UPDATE_rates_complete(G, "ignore")
    
        # set growth directions of growing vertices    
        for v in G.V[1:]:
            if v.status_growing or v.status_branching:
                pos =  v.position
                parent =  G.V[v.in_edge]
                par_pos = parent.position
                v.growth_direction = PROCESS_structures.MINUS(pos, par_pos)
            
        
    # sets up the EVENT_manager
    G.EVENT_manager = PROCESS_structures.EVENT_manager(G)
    G.graph_rate = G.EVENT_manager.graph_rate
    return G



#####################
#####################

# Store and extract-functions for ancestral trees
#################################################

def save_tree_sql(tree, filename, fullpathname = False):
    """ Creates an sql-database and stores a genealogical tree """
    
    # Creates directory for storage
    if not fullpathname:
        dirName = create_folder(filename)
    else:
        dirName = filename
        Path(dirName).mkdir(parents=True, exist_ok=True)
    
    # prepare the database
    name = str ( dirName + "/" + "tree.db" )
    try:
        os.remove(name)
    except:
        pass
    db = sqlite3.connect(name)
    cursor = db.cursor()    
    
    # prepare database
    ##############################
    # create node-table
    cursor.execute( """ DROP TABLE if exists nodes """)
    create_node_table = """   CREATE TABLE if not exists nodes ( 
                                            name TEXT,
                                            fitness FLOAT,
                                            mean_fitness FLOAT,
                                            max_fitness FLOAT,
                                            parent INTEGER,
                                            children TEXT,
                                            weight INTEGER,
                                            branch_weight INTEGER,
                                            generation INTEGER,
                                            number INTEGER,
                                            fittest_type TEXT,
                                            group_node INTEGER,
                                            colour TEXT,
                                            spatial_center_x FLOAT,
                                            spatial_center_y FLOAT,
                                            spatial_center_z FLOAT,
                                            is_found INTEGER
                                            ) """
    
    cursor.execute( create_node_table )
    db.commit()
    
    # create group-node-table
    cursor.execute ( """ DROP TABLE if exists group_nodes """ )
    cursor.execute ( """ CREATE TABLE if not exists group_nodes (
                            name TEXT,
                            number INTEGER
                            ) """ )
    db.commit()
    
    # create information-table
    cursor.execute ( """ DROP TABLE if exists tree_information """ )
    create_info_table = """ CREATE TABLE tree_information 
                            (
                        height  INTEGER,
                        ancestor_generation INTEGER,
                        type_amount INTEGER,
                        total_weight INTEGER,
                        root INTEGER,
                        unclassified_maximal_fitness FLOAT,
                        parameters TEXT
                        ) """
    cursor.execute(create_info_table)
    db.commit()
    
    
    
    # fill the database with all sorts of things
    ##############################
    
    # store nodes
    insert_node_comm = """ INSERT INTO nodes (
                                            name,
                                            fitness,
                                            mean_fitness,
                                            max_fitness,
                                            parent,
                                            children,
                                            weight,
                                            branch_weight,
                                            generation,
                                            number,
                                            fittest_type,
                                            group_node,
                                            colour,
                                            spatial_center_x,
                                            spatial_center_y,
                                            spatial_center_z,
                                            is_found )
                        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) """
    def insert_node(n):
        # convert children to readable string
        children = json.dumps( [ child.number for child in n.children ] )
        
        # converts tht is_found-information to int -1,0,1 = None, False, True
        if n.is_found != None:
            is_found = int( n.is_found )
        else:
            is_found = -1
            
        # stores the number of the group_node if defined
        if n.group == None:
            group_node = -2
        elif n.group == False:
            group_node = -1
        else:
            group_node = n.group.number
            
        # stores the spatial center if defined
        if n.spatial_center:
            spatial_center = n.spatial_center
        else:
            spatial_center = [None, None, None]
        
        # insert into database
        cursor.execute ( insert_node_comm, [
                                            json.dumps(n.name),
                                            n.fitness,
                                            n.mean_fitness,
                                            n.max_fitness,
                                            n.parent.number,
                                            children,
                                            n.weight,
                                            n.branch_weight,
                                            n.generation,
                                            n.number,
                                            json.dumps(n.fittest_type),
                                            group_node,
                                            json.dumps(n.colour),
                                            spatial_center[0],
                                            spatial_center[1],
                                            spatial_center[2],
                                            is_found
                         ] )
    
    for v in tree.V.values():
        insert_node(v)
    db.commit()
    
    # store information
    insert_info_command = """ INSERT INTO tree_information (
                        height,
                        ancestor_generation,
                        type_amount,
                        total_weight,
                        root,
                        unclassified_maximal_fitness,
                        parameters )
                            VALUES (?,?,?,?,?,?,?) """
    cursor.execute ( insert_info_command, 
                            [
                            tree.height,
                            tree.ancestor_generation,
                            tree.type_amount,
                            tree.total_weight,
                            tree.root.number,
                            tree.unclassified_maximal_fitness,
                            json.dumps(tree.parameters)
                            ])
    
    # store the group-nodes (if a clustering was done before)
    if tree.group_points != None:
        for v in tree.group_points:
            cursor.execute ( """INSERT INTO group_nodes (
                                    name, number)
                                    VALUES(?,?) """ ,[
                                    json.dumps(v.name),
                                    v.number
                                    ])
    
    db.commit()

    # history of strong mutations
    cursor.execute ( """ DROP TABLE if exists strong_mutations_log """)
    insert_strong_mutations_log = """ CREATE TABLE strong_mutations_log (
                                        mutation_key TEXT,
                                        mutation_time FLOAT
                                    ) """
    cursor.execute (insert_strong_mutations_log)
    insert_log_comm = """ INSERT INTO strong_mutations_log (
                            mutation_key,
                            mutation_time )
                        VALUES(?,?) """
    for key in tree.strong_mutations_log.keys():
        time_of_mut = tree.strong_mutations_log[key]
        cursor.execute( insert_log_comm , [ json.dumps(key), float(time_of_mut) ] )
    db.commit()
    

    db.close()
    
    print("Saved tree at ", name)
    
    return


class OperationalError(Exception):
    pass


def get_tree_sql(filename, fullpathname = False):
    """ Loads a gen tree from either the output/filename
        or from the fullpath=filename if fullpathname == True """
    
    if fullpathname:
        dirName = filename
    else:
       # go to output/filename
        directory = os.getcwd()
        dirName = str ( directory + "/output/" + filename )
    
    # find tree        
    for f in absoluteFilePaths(dirName):
        if f.endswith("tree.db"):
            file_name = f
            break

    # get database
    db = sqlite3.connect(file_name)
    cursor = db.cursor()
    
    # create tree
    tree = EVAL_gen_tree_graph.ancestral_tree()

    cursor.execute ( """ SELECT 
                            height,
                            ancestor_generation,
                            type_amount,
                            total_weight,
                            root,
                            unclassified_maximal_fitness,
                            parameters FROM tree_information """)
    info = cursor.fetchall()[0]

    # unpack general information
    tree.height = info[0]
    tree.ancestor_generation = info[1]
    tree.type_amount = info[2]
    tree.total_weight = info[3]
    tree.root =  info[4] # later gets replaced by actual node
    tree.unclassified_maximal_fitness = info[5]
    tree.parameters = json.loads( info[6] )

    # prepare tree for nodes
    tree.V = dict()
    tree.generations = []
    for i in range ( tree.height ):
        tree.generations.append ( set () )
    
    tree.type_of_number = []
    for i in range ( tree.type_amount ):
        tree.type_of_number.append( None )
    
    # unpack nodes
    # get info
    cursor.execute ( """ SELECT
                                                name,
                                                fitness,
                                                mean_fitness,
                                                max_fitness,
                                                parent,
                                                children,
                                                weight,
                                                branch_weight,
                                                generation,
                                                number,
                                                fittest_type,
                                                group_node,
                                                colour,
                                                spatial_center_x,
                                                spatial_center_y,
                                                spatial_center_z,
                                                is_found
                        FROM nodes """ )
    nodes = cursor.fetchall()

    # insert into gen-tree
    for node in nodes:
        name = tuple (  json.loads(node[0]) )
        v = EVAL_gen_tree_graph.node(name)
        v.fitness = node[1]
        v.mean_fitness = node[2]
        v.max_fitness = node[3]
        v.parent = node[4]
        v.children = json.loads(node[5])
        v.weight =  node[6]
        v.branch_weight = node[7]
        v.generation = node[8]
        v.number = node[9]
        v.fittest_type = tuple ( json.loads(node[10]) )
        
        # if a clustering was made, the nodes get the number of the corresponding
        # clustering-node (later replaced by node itself)
        if node[11] == -2:
            v.group_node = None
        if node[11] == -1:
            v.group_node = False
        else:
            v.group_node = node[11]
            
        v.colour = json.loads( node[12] )
        
        # if info about the spatial center is stored, it gets unpacked
        if node[13] == None:
            v.spatial_center = None
        else:
            v.spatial_center = ( node[13], node[14], node[15] )
        
        # if a biopsy was made, the it is stored if the vertex was found or not
        if node[16] == -1:
            v.is_found = None
        elif node[16] == 0:
            v.is_found = False
        else:
            v.is_found = int(node[16])
            
        tree.V[v.name] = v
        tree.type_of_number[v.number] = v
        tree.generations[v.generation].add(v)
        
    # replace root-key by actual node
    tree.root = tree.type_of_number[tree.root]
    
    # insert info about the group-nodes
    cursor.execute ( """ SELECT number FROM group_nodes """)
    group_point_keys = cursor.fetchall()
    group_point_keys = [ key[0] for key in group_point_keys ]
    tree.group_points = set()
    for i in group_point_keys:
        if type(i) != int:
            print(i)
        tree.group_points.add( tree.type_of_number[i] )
    
    # second iteration through nodes: replace keys by nodes
    for v in tree.V.values():
        v.parent = tree.type_of_number[ v.parent ]
        v.children = [ tree.type_of_number[w] for w in v.children  ]
            
        # sets the group if the vertex is clustered
        if v.group:
            v.group = tree.type_of_number[v.group]
            
    # get history of strong mutations
    # get strong_mut_log
    cursor.execute ( """ SELECT mutation_key, mutation_time
                        FROM strong_mutations_log """)
    values = cursor.fetchall()
    for str_mut in values:
        key  = json.loads( str_mut[0] )
        time = str_mut[1]
        tree.strong_mutations_log[tuple(key)] = time
    del(values)
    
    db.close()
    return tree

#%%
def gen_clusters_to_excel(tree, filename, store_graph_info = True, fullpathname = False):
    """ IN: Evaluated (clustered, biopsy) gen_tree
        OUT: Writes the information about the clusters into an excel-sheet """
        
    if not fullpathname:
        dirname = create_folder(filename)
    else:
        dirname = filename
        Path(dirname).mkdir(parents=True, exist_ok=True)
    
    if tree.group_points == None:
        raise ValueError("Writing to excel failed, gen_tree not evaluated")
        
    
    # set of clustered genotypes
    genotypes = tree.group_points
    
    # extract different information to lists
    names = list ( v.name for v in genotypes )
    indices = [ i for i in range(1,len(genotypes)+1) ]
    sizes = list ( v.branch_weight for v in genotypes )
    sizes_percent = [ round ( 100.0 * size / tree.total_weight,2) for size in sizes ]
    max_fitnesses = list ( v.max_fitness for v in genotypes )
    fittest_types = list ( v.fittest_type for v in genotypes  )
    centroids = list( v.spatial_center for v in genotypes ) 
    mean_fitnesses = list ( round(v.mean_fitness,2) for v in genotypes )
    found_in_biopsy = list ( v.is_found for v in genotypes )
    
    
    
    root = tree.root    
    # one empty row, unclustered types, root=common ancestor
    names.append( None )
    names.append( None )
    names.append(root.name)
    
    indices.append( None)
    indices.append( "Non-classified")
    indices.append( "Common Ancestor")
    
    sizes.append ( None )
    uncl = tree.total_weight - sum ( v.branch_weight for v in genotypes)
    sizes.append ( uncl )
    sizes.append(root.branch_weight)
    
    sizes_percent.append ( None )
    sizes_percent.append ( round(100.0 * uncl / tree.total_weight, 2) )
    sizes_percent.append ( round(100.0 * root.branch_weight / tree.total_weight, 2) )
    
    max_fitnesses.append ( None )
    max_fitnesses.append ( tree.unclassified_maximal_fitness )
    max_fitnesses.append ( root.max_fitness )
    
    centroids.append (None)
    centroids.append (None)
    centroids.append (None)
    
    mean_fitnesses.append(None)
    mean_fitnesses.append(None)
    mean_fitnesses.append( round(root.mean_fitness,2) )
    
    fittest_types.append (None)
    fittest_types.append(None)
    fittest_types.append(None)
    
    for i in range( len(found_in_biopsy) ):
        if found_in_biopsy[i] == 0:
            found_in_biopsy[i] = "Not found"
    found_in_biopsy = found_in_biopsy + [None, None, None]
    
    # create panda dataframe for easy usage
    DF = DataFrame({
            "Index": indices,
            "Num of cells": sizes,
            "Percent": sizes_percent,
            "Genetic Formula": names,
            "Fittest types": fittest_types,
            "Maximal fitness": max_fitnesses,
            "Mean fitness": mean_fitnesses,
            "Spatial Center": centroids,
            "Found in Biopsy": found_in_biopsy
            
            })
    
    gen_name = dirname + "/" + "clustered_genotypes.xlsx"
    param_name = dirname + "/" + "specifications.xlsx"
    #tree.parameters["common_ancestor"] = tree.root.name
    
    DF.to_excel(gen_name,index=False)
    
    if store_graph_info:
        spec = list (tree.parameters.keys() )
        para = list (tree.parameters.values() )
        PM = DataFrame({
                "Specification": spec,
                "Parameter": para
                })
        PM.to_excel(param_name, index = False)
        
        # output for strong-mutation-log
        log = tree.strong_mutations_log
        names = list ( str(k) for k in log.keys() )
        times = list ( round(float(t),2) for t in log.values() )
        HL = DataFrame ( {"Name of mutation": names, "Time of mutation": times})
        log_name = dirname + "/" + "strong_mut_log.xlsx"
        HL.to_excel(log_name, index = False) 
    
    return
  
#%%

def all_mutations_to_excel(tree, filename, fullpathname = False, num_mut_info = True, mutations_in_biopsy_info = False):
    """ IN: gen_tree
        OUT: Writes information about all mutations to excel-sheet """
    
    if not fullpathname:
        dirname = create_folder(filename)
    else:
        dirname = filename
        Path(dirname).mkdir(parents=True, exist_ok=True)
    
    if tree.V == None:
        raise ValueError("Writing to excel failed, gen_tree has no mutations")
        
    
    # set of clustered genotypes
    genotypes = list ( tree.V.values() )
    genotypes.sort(key = lambda v: v.branch_weight)
    genotypes.reverse()
    
    # extract different information to lists
    names = list ( v.name for v in genotypes )
    indices = [ i for i in range(1,len(genotypes)+1) ]
    sizes = list ( v.branch_weight for v in genotypes )
    sizes_percent = [ round ( 100.0 * size / tree.total_weight,3) \
                     for size in sizes ]
    num_mutations = list ( len(v.name) for v in genotypes  )
    max_fitnesses = list ( round(v.max_fitness,2) for v in genotypes )
    fittest_types = list ( v.fittest_type for v in genotypes  )
    max_mutations = list (  len(v) for v in fittest_types )
    mean_fitnesses = list ( round(v.mean_fitness,2) for v in genotypes )
    
    if mutations_in_biopsy_info :
        found_in_biopsy = list ( tuple(v.is_found) for v in genotypes )
    else:
        found_in_biopsy = list ( 1.0 for v in genotypes )
    
    empty_column = [ None for v in genotypes ]
    
    all_columns = [ names, indices, sizes, sizes_percent, max_fitnesses,\
        fittest_types, mean_fitnesses, max_mutations,\
             num_mutations, empty_column, found_in_biopsy]
      
    for col in all_columns:
        col.insert(1, None)
    K,L = 1,0
    for size_threshold in [10.0, 1.0, 0.1, 0.01, 0.001]:
        checks = [i for i in range(K+1,len(genotypes)) \
            if sizes_percent[i] >= size_threshold ]
        if checks != []:
            L = checks[-1]
        else:
            L = K
        if L != K:
            L+=1
            for col in all_columns:
                col.insert(L, None)
            K = L
    
    # create panda dataframe for easy usage
    DF = DataFrame({
            "Index": indices,
            "Percent": sizes_percent,
            "Cells": sizes,
            "Mutations": num_mutations,
            "Mut of fittest": max_mutations,
            "": empty_column,
            "Mean fitness": mean_fitnesses,
            "Maximal fitness": max_fitnesses,
            "Genetic Formula": names,
            "Fittest types": fittest_types,
            "Found in Biopsy": found_in_biopsy
            })
    
    all_mut_name = dirname + "/" + "all_mutations.xlsx"
    
    DF.to_excel(all_mut_name,index=False)
    
    if num_mut_info:
    # also write additional hierarchy-information
        num_mutations_to_excel(tree, filename, fullpathname)
    
    return
    
#%%
    
def num_mutations_to_excel(tree, filename, fullpathname = False):
    """ IN: Evaluated gen_tree
        OUT: Writes the information about the clusters into an excel-sheet """
    
    if not fullpathname:
        dirname = create_folder(filename)
    else:
        dirname = filename
        Path(dirname).mkdir(parents=True, exist_ok=True)
    
    num_mut_name = dirname + "/" + "num_mut.txt"
    
    if tree.V == None:
        raise ValueError("Writing to excel failed, gen_tree has no mutations")
        
    f = open(num_mut_name, 'w') 
    S = 0
    for gen in tree.generations:
        for v in gen:
            num_mut = len(v.name)
            break
        S += sum ( v.weight for v in gen )
        print("At most", num_mut, "mutations:", S, file = f)

    print("\n", file = f)
    print("More then X mutations", file = f)
    S = sum ( v.weight for v in tree.V.values() )
    for gen in tree.generations:
        for v in gen:
            num_mut = len(v.name)
            break
        S -= sum ( v.weight for v in gen )
        print(num_mut, ":", S, file = f)
            
    f.close()
    
    return

#%%
""" INPUT/OUTPUT concerning biopsies """
def biopsy_single_store(sample, sample_name, graph_directory, name_prefix = False):
    """ Stores a biopsy as Excel Sheet (via panda DF) """
    
    keys = [ v.key for v in sample]
    # create panda dataframe for easy usage
    DF = DataFrame({
            "Vertex Keys": keys
            })
    
    if not name_prefix:
        biopsydir = graph_directory + "/biopsies"
    else:
        biopsydir = graph_directory + "/" + str(name_prefix)
    Path(biopsydir).mkdir(parents=True, exist_ok=True)
    
    filename = biopsydir + "/" + sample_name + ".xlsx"
    DF.to_excel(filename, index=False)
    
def biopsies_store(samples, graph_directory, name_prefix = False):
    """ Stores all biopsies given via samples:
        Stores the vertex-keys as xlsx """
    
    names = ["biopsy_" + str(i) +"/" + "biopsy_" + str(i) for i in range(1, len(samples)+1)]
    
    for i in range(len(samples)):
        s = samples[i]
        n = names[i]
        
        biodir = graph_directory + "/" + str(name_prefix) + "/biopsy_" + str(i+1) +"/"
        Path(biodir).mkdir(parents=True, exist_ok=True)
    
        biopsy_single_store(s, n, graph_directory, name_prefix=name_prefix)

        
def biopsies_get(G, graph_directory, name_prefix = False):
    """ Load biopsies of a given graph and a directory
        Graph as info is needed because the biopsy-files
        give only keys of vertices"""
    if not name_prefix:
        biopsydir = graph_directory + "/biopsies" + "/storage"
    else:
        biopsydir = graph_directory + "/" + str(name_prefix) + "/"
    samples = []
    for file in os.scandir(biopsydir):
        if file.path.endswith(".xlsx"):
            print(file.path)
            table = pd.read_excel(file)
            keys = table["Vertex Keys"]
            vertices = [G[i] for i in keys]
            samples.append(vertices)
    
    return samples
