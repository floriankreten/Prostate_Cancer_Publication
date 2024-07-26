#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 11:23:20 2019

@author: floriankreten
"""

""" A module where all important evaluation-schemes are bundled
    The data is typically searched and stored in ./output/ next to the main.py
    For loading data from other sources, use the "fullpathname" argument 
"""

import EVAL_biopsy
import EVAL_data_output
import EVAL_type_clustering_print
import EVAL_graph_plot_plotly
import EVAL_gen_tree_graph
import EVAL_gen_tree_output_simple
import EVAL_tree_metrics

import PROCESS_structures

import gc
import numpy as np
import pandas as pd
from ast import literal_eval
    


class table_of_evaluation:
    """ A class that brings together all possible methods for
        storage, graphical and excel-output
        needs a graph G as input
        can also load a pre-defined tree (generates it otherwise)
        For an impression, load a graph and call .complete_eval_and_storage() 
        For handling also biopsies, see child-> biopsy_evaluation_table below """
        
    def __init__(self, filename: str):
        self.graph = None
        self.gen_tree = None
        self.clusters_found = None
        self.biopsies = None
        self.all_biopsies = None
        
        self.graphnodes_cluster_list = None
        self.vertices_by_cluster_groups = None
        self.cluster_colours = None
        
        self._filename = filename
        
    def get_graph(self, graph):
        """ Loads a graph G into the table """
        self.graph = graph
        
    def get_graph_from_path(self):
        """ Load stored graph from a specific location """
        self.graph = EVAL_data_output.get_graph_sql(self._filename, fullpathname=True)
    
    def get_gen_tree(self, find_size=0.05, find_type="relative"):
        """ Creates a gen-graph, minimal size of mutation is 
            either relative or absolute """
        if self.graph == None:
            raise ValueError("Insert graph first")
        self.graph, self.graphnodes_cluster_list, self.gen_tree = \
            EVAL_type_clustering_print.classify_nodes(
                        self.graph, find_size, find_type)
        if self.biopsies:
             self._compare_biopsy_and_clonal_table()
             
    def _vertices_of_type(self):
        """ Lists of vertices, ordered by genotypes (defined in anc-tree) """
        if not self.gen_tree:
            assert False, "Define tree before trait-plot (for labeling vertices)"

        group_list = self.graphnodes_cluster_list
        lists = [ [] for i in group_list]

        new_list = [ v for v in self.graph.V if isinstance(v.cluster_group, int) ]
        
        for v in new_list:
            group = v.cluster_group
            lists[group].append(v)
        
        self.vertices_by_cluster_groups = lists
        
    def spatial_trait_and_anc_plot(self, title = "", show = False, 
                                       pre_subset = False,
                                       fullpathname = False,
                                       show_tree = True):
        
        """ Spatial plot of different genotypes
            Sketch only, 3D-plot engine can not handle millions of nodes """
        
        # order clustered vertices
        self._vertices_of_type()
        lists = self.vertices_by_cluster_groups
        
        # select for pre-subset if wanted
        if pre_subset:
            new_lists = []
            
            for L in lists:
                new_list = [v for v in L if v in pre_subset ]
                new_lists.append(new_list)
                
            lists = new_lists
                
            
        # thin out plot if necessary (engine gets problems for >50k nodes)
        subset_plot = []
        wanted_size = 50000.0
        current_size = sum ( len(l) for l in lists)
        quotient = wanted_size / current_size
        
        if quotient >= 1.0:
            for L in lists:
                subset_plot += L
        
        else:
            for i in range(len(lists)):
                gen_list = lists[i]
                S = len(gen_list)
                length_wanted = int(S * quotient) + 1
                length_wanted = max(length_wanted, 200)
                if length_wanted > S:
                    V = gen_list
                else:
                    V = np.random.choice(a = gen_list, size = length_wanted)
                    V = np.transpose(V)
                    V = list(V)
                    
                #print("num of vertices in reduced plot:", len(V))
                subset_plot += V
        
        EVAL_graph_plot_plotly.subset_plot(self.graph, subset_plot, self._filename,
                            show = show, color_key = self.cluster_colours,
                            title = title, fullpathname = fullpathname)
        print("For graph-plot, we use colors ", self.cluster_colours)
        
        
        if show_tree:
            print("Printing phylogeny according to spatial plot")
            EVAL_gen_tree_output_simple.print_tree_simple(self.gen_tree,
                                                          self._filename,
                                                          fullpathname =True,
                    minimal_size = 0.05, size_mode= "relative",
                    biopsy_extra = False,  use_colors=True)
    
    def get_circular_biopsies(self, num = 2):
        """ DEPRECATED
            Takes virtual biopsies
            Biopsies are opposing balls, either 2 or 3 """
        if not num in [2,3]:
            raise ValueError("Number of biopsies must be 2 or 3, gave", num)
        if self.graph == None:
            raise ValueError("Insert graph first")
        self.graph, self.biopsies = EVAL_biopsy.take_circular_biopsy(
            self.graph, diameter = 90, amount = num)
        if self.gen_tree:
            self._compare_biopsy_and_clonal_table()
            
    def evaluate_tree_metrics(self):
        if not self.gen_tree:
            raise ValueError("Tree needed for tree metrics")
        D, J = EVAL_tree_metrics.calculace_DiversityD_and_BalanceJ(self.gen_tree)
        return D, J

        
    def sketch_biopsy_clones(self, show = False, fullpathname = False):
        """ A sketch of the biopsy within the tumor """
        if not self.all_biopsies:
            self._throw_biopsies_together()
            
        self.spatial_trait_and_anc_plot(title = "biopsy_clonal", show = show, 
                                       pre_subset = self.all_biopsies,
                                       fullpathname = fullpathname)
        
    def sketch_biopsy_plot(self, show = False, fullpathname = False):
        """ A sketch of the biopsy within the tumor """
        if not self.all_biopsies:
            self._throw_biopsies_together()
        EVAL_biopsy.show_biopsy(
            self.graph, self.biopsies, self._filename, show, fullpathname=fullpathname)
        
    def ancestral_graph_plot(self, fullpathname=False, plotname = False):
        """ Prints the ancestral tree / graph
            Marks all genotypes that are not found in the biopsies """
        if not self.gen_tree:
            print("Calculating ancestral tree ...")
            self.get_gen_tree()
        
        name = self._filename
        self.gen_tree.cluster_and_biopsy_print(
            filename = name, fullpathname =fullpathname, plotname=plotname)
        
    def clonal_centers_plot(self, show = False):
        """ DEPRECATED
            Clonal centers are plotted as balls in R³ """
        if not self.gen_tree:
            print("Calculating ancestral tree ...")
            self.get_gen_tree()
        EVAL_type_clustering_print.clonal_centers_plot(
            self.gen_tree, self._filename, show)
        
    def save_genotypes(self, store_graph_info = True, fullpathname=False):
        """ Stores information of the clustered genotypes into an excel-file"""
        EVAL_data_output.gen_clusters_to_excel(self.gen_tree, self._filename, store_graph_info,
                                              fullpathname=fullpathname)
        
    def save_graph(self, fullpathname = False):
        """ Stores the graph
            WARNING: Can take up to 10gb for very large graphs """
        EVAL_data_output.save_graph_sql(self.graph, self._filename,
                                        fullpathname= fullpathname )
        
    def store_all_mutations(self, fullpathname= False, cutoff = 0.02):
        """ Stores information about all mutations in a nice ancestral tree.db """
        if not self.graph:
            raise ValueError("Define graph pls!")
        tree = EVAL_gen_tree_graph.create_gen_graph(self.graph, find_size=cutoff, find_type="relative")
        EVAL_data_output.all_mutations_to_excel(tree, self._filename, fullpathname=fullpathname)
        EVAL_data_output.save_tree_sql(tree, self._filename, fullpathname=fullpathname )
        all_mut_enhancer(self._filename + "/all_mutations.xlsx")
        del(tree)
        gc.collect()
        
        
        
    def complete_eval_and_storage(self):
        """ Easy to use rather complete evaluation of a graph
            Needs only self.graph and self.filename to be defined """
        self.get_gen_tree()
        #self.get_circular_biopsies()
        #self.complete_mut_pattern()
        #self.sketch_biopsy_plot
        self.ancestral_graph_plot()
        
        title = "size" + str(round(self.graph.count/1000000)) + "Mio |  radius " + str(int( self.graph.appr_radius() ))
        self.spatial_trait_and_anc_plot(title=title, fullpathname = True)
        
        self.save_graph()
        self.save_genotypes()
        self.store_all_mutations()
        
    def _compare_biopsy_and_clonal_table(self):
        self.gen_tree = EVAL_biopsy.compare_biopsy_and_clustering(
            self.biopsies, self.gen_tree)
        self.cluster_colours = list ( v.colour for v in self.gen_tree.group_points )
        
    def _throw_biopsies_together(self):
        if not self.biopsies:
            self.get_circular_biopsies()
        self.all_biopsies = set().union(*self.biopsies)



#%%
    
class biopsy_evaluation_table(table_of_evaluation):
    """ Advanced version of the table_of_evaluation that is used
        for making and evaluating the 5 punching biopsies
        again needs a graph as input
        use .precise_clonal_eval for an impression/complete call"""
        
    def get_punching_biopsies(self):
        """ Make punching biopsies of a given graph """
        if self.graph == None:
            raise ValueError("Insert graph first")
        self.graph, self.biopsies = EVAL_biopsy.take_punching_biopsy(
                                                        self.graph)
        if self.gen_tree:
            self._compare_biopsy_and_clonal_table()
            
    def get_level_circ_biopsies(self, target_size, heights=[0.5],
                                thickness = False, num_samples=[1]):
        if self.graph == None:
                    raise ValueError("Insert graph first")
        self.graph, self.biopsies = EVAL_biopsy.take_level_circle_biospy(
                        G=self.graph, target_size = target_size,
                        heights = heights, num_samples = num_samples,
                        thickness=thickness)
        
        if self.gen_tree:
            self._compare_biopsy_and_clonal_table()
            
    def take_five_quad_circ_biopsies(self, target_size = 500, heights = [0.5],
                               num_samples = [4], ratio_diam_h = 15):
        if self.graph == None:
                    raise ValueError("Insert graph first")
        self.graph, self.biopsies = EVAL_biopsy.take_five_quad_circ_biopsies(
                        G=self.graph, target_size=target_size,
                        heights=heights, num_samples = num_samples,
                        ratio_diam_h = ratio_diam_h)
        
    def take_five_quad_biopsies_with_center_cut(self, target_size = 500, heights = [0.5],
                               num_samples = [4], ratio_diam_h = 15):
        if self.graph == None:
                    raise ValueError("Insert graph first")
        self.graph, self.biopsies = EVAL_biopsy.take_five_quadratic_biopsies_with_center_cut(
                        G=self.graph, target_size=target_size,
                        heights=heights, num_samples = num_samples,
                        ratio_diam_h = ratio_diam_h)
            
    def get_chips_biopsy(self, target_size=500, heights=[0.5],
                             num_samples = [1], ratio_rad_h = 15):
        if self.graph == None:
                    raise ValueError("Insert graph first")
        self.graph, self.biopsies = EVAL_biopsy.take_chips_biopsy(
                        G=self.graph, target_size = target_size,
                        heights = heights, num_samples = num_samples,
                        ratio_rad_h = ratio_rad_h)
        
    def take_4_quad_on_same_level_biopsies(self, target_size=500, heights=[0.5],
                             num_samples = [1], ratio_diam_h = 15):
        if self.graph == None:
                    raise ValueError("Insert graph first")
        self.graph, self.biopsies = \
                    EVAL_biopsy.take_4_quad_on_same_level_biopsies(
                        G=self.graph, target_size = target_size,
                        heights = heights, num_samples = num_samples,
                        ratio_diam_h = ratio_diam_h)
                    
    def take_4_quad_on_different_levels_biopsy(self, target_size=500,
                heights=[0.35, 0.65], num_samples = [2,2], ratio_diam_h = 15):
        if self.graph == None:
                    raise ValueError("Insert graph first")
        self.graph, self.biopsies = \
                    EVAL_biopsy.take_4_quad_on_different_levels_biopsy(
                        G=self.graph, target_size = target_size,
                        heights = heights, num_samples = num_samples,
                        ratio_diam_h = ratio_diam_h)
        
    
    def get_gen_tree_from_directory(self,
                                    find_size=0.05,
                                    find_type="relative",
                                    cluster_type = "simple_branched"):
        
        """ Creates a new ancestral tree for a given graph
            Uses the stored tree from the given path to save computation time """
            
        if self.graph == None:
            raise ValueError("Insert graph first")
            
        self.graph, self.graphnodes_cluster_list, self.gen_tree = \
            EVAL_type_clustering_print.classify_nodes_from_directory(
                        G = self.graph,
                        dirname = self._filename,
                        find_size=find_size,
                        find_type=find_type,
                        cluster_type = cluster_type)
            
        self.cluster_colours = list ( v.colour for v in self.gen_tree.group_points )
        
        if self.biopsies:
             self._compare_biopsy_and_clonal_table()
             
    def compare_biopsy_and_all_mutations(self):
        self.gen_tree = EVAL_biopsy.compare_biopsy_and_all_mutations(
            self.biopsies, self.gen_tree)
            
    def store_all_mutations_of_biopsies(self, biopsy_name="default_biopsy"):
        """ Stores information about mutations in biopsies """
        if not self.graph:
            raise ValueError("Define graph pls!")
        
        for i in range(1, len(self.biopsies) +1 ):
            subset = self.biopsies[i-1]
            biopsydir = self._filename + "/" + biopsy_name + "/biopsy_" + str(i)
            minigraph = PROCESS_structures.mini_biopsy_graph(subset, self.graph)
            minitree = EVAL_gen_tree_graph.create_gen_graph(minigraph, find_size=1, find_type="absolute")
            EVAL_data_output.all_mutations_to_excel(minitree, biopsydir,
                                            fullpathname=True, num_mut_info = False)
            del(minitree, minigraph)
            gc.collect()
            
        if self.all_biopsies:
            subset = self.all_biopsies
            biopsydir = self._filename + "/" + biopsy_name + "/combined_biopsies"
            minigraph = PROCESS_structures.mini_biopsy_graph(subset, self.graph)
            minitree = EVAL_gen_tree_graph.create_gen_graph(minigraph, find_size=1, find_type="absolute")
            EVAL_data_output.all_mutations_to_excel(minitree, biopsydir,
                                            fullpathname=True, num_mut_info = False)
            del(minitree, minigraph)
            gc.collect()
            
        
            
    def precise_clonal_eval(self, load_graph = False, biopsy="level",target_size = 66000):
        """ WARNING: Needs a gen-tree stored at self._filename
            See table_of_evaluation above for further information
            Evaluates the ancestral graph
            Does a biopsy
            Compares biopsy and graph: how much info about genotypes
                do we get from biopsies """
        
        # values are chosen for 20Mln graphs and 70000 vertex nodes per biopsy
        
        if load_graph:
            self.get_graph_from_path()
        
        if biopsy == "level":
            print("Doing level-circ biopsy")
            self.get_level_circ_biopsies(target_size=target_size,
                                         heights=[0.3,0.7], num_samples = [3,2],
                                         thickness =  2 * self.graph.parameters["vertex_radius"])    
            
        elif biopsy == "quad_circ":
            print("Doing quad_circ level biopsy")
            self.take_five_quad_circ_biopsies(target_size=target_size,
                                         heights=[0.5], num_samples = [5],
                                         ratio_diam_h = 15)
            
        elif biopsy == "quad_center_cut":
            print("Doing quad with center cut biopsy")
            self.take_five_quad_biopsies_with_center_cut(target_size=target_size,
                                         heights=[0.5], num_samples = [5],
                                         ratio_diam_h = 15)
        
        elif biopsy == "4_quad_on_same_level":
            print("Doing 4 quad on same level biopsy")
            self.take_4_quad_on_same_level_biopsies(target_size=target_size,
                                         heights=[0.5], num_samples = [4],
                                         ratio_diam_h = 15)
            
        elif biopsy == "4_quad_on_different_levels":
            print("Doing 4 quad on different levels biopsy")
            self.take_4_quad_on_different_levels_biopsy(target_size=target_size,
                                         heights=[0.35, 0.65], num_samples = [2,2],
                                         ratio_diam_h = 15)
            
        
        else:
            self.get_punching_biopsies()
            
        self.get_gen_tree_from_directory()
        
        self.sketch_biopsy_plot(fullpathname = True)
        self.ancestral_graph_plot(fullpathname=True)
        
        # do not change the order
        self.sketch_biopsy_clones(fullpathname = True)
        self.spatial_trait_and_anc_plot(fullpathname = True) 
        
        self.store_punching_biopsies(fullpathname=True)
        self.save_genotypes(fullpathname=True)
        self.store_all_mutations_of_biopsies()
        
        # second evaluation: rough upper part of the tree
        self.get_gen_tree_from_directory(
                                    find_size=0.05,
                                    find_type="relative",
                                    cluster_type = "rough")
        self.ancestral_graph_plot(fullpathname=True, plotname = "tree_rough_biopsy")
        print("Clonal eval complete")
        
        
    def clonal_eval_2023(self, load_graph = False, biopsy="level",target_size = 66000):
        """ WARNING: Needs a gen-tree stored at self._filename
            See table_of_evaluation above for further information
            Evaluates the ancestral graph
            Does a biopsy, tries to match target_size for each biopsy-core
            Compares biopsy and graph: how much info about genotypes
                do we get from biopsies """
                        
        
        # values are chosen for 20Mln graphs and 70000 vertex nodes per biopsy
        
        if load_graph:
            self.get_graph_from_path()
        
        elif biopsy == "4_quad_on_same_level":
            print("Doing 4 quad on same level biopsy")
            self.take_4_quad_on_same_level_biopsies(target_size=target_size,
                                         heights=[0.5], num_samples = [4],
                                         ratio_diam_h = 15)
            
        elif biopsy == "4_quad_on_different_levels":
            print("Doing 4 quad on different levels biopsy")
            self.take_4_quad_on_different_levels_biopsy(target_size=target_size,
                                         heights=[0.35, 0.65], num_samples = [2,2],
                                         ratio_diam_h = 15)

        else:
            raise ValueError("Only two types of biopsies allowed: \
                             4_quad_on_same_level, 4_quad_on_different_levels")
                             
        self._throw_biopsies_together()
            
        # Update all filenames to biopsy-type, for storage only
        name_save = self._filename
        def update_name_to_biopsy():
            self._filename = self._filename + "/" + biopsy
        def revert_name():
            self._filename = name_save
        
        # load gen-tree from directory and insert biopsy-info
        self.gen_tree = EVAL_data_output.get_tree_sql(self._filename, fullpathname=True)
        

        print("Updating gen-tree with biopsy")
        self.compare_biopsy_and_all_mutations()
        
        
        # 3D plot of the biopsy within the tumor
        # All_mutations file with infos about mutations in biopsies
        update_name_to_biopsy()
        EVAL_data_output.all_mutations_to_excel(self.gen_tree, self._filename,
                fullpathname=True, mutations_in_biopsy_info = True,
                num_mut_info = False) # works like a charme
        self.sketch_biopsy_plot(fullpathname = True)
        EVAL_gen_tree_output_simple.print_tree_simple(self.gen_tree, self._filename, fullpathname = True,
                minimal_size = 0.01, size_mode= "relative", biopsy_extra = True) # works like a charme
        revert_name()
    
        
        # store raw biopsy info
        EVAL_data_output.biopsies_store(self.biopsies, self._filename, name_prefix = biopsy)
        
        # store excel-sheets with mutation-infos for each biopsy
        self.store_all_mutations_of_biopsies(biopsy)
        
        
        print("Clonal eval complete")
        
#%%    
        
def evaluation_routine(G, filename):
    """ Input: Graph G and filename
        Does an evaluation of the evolution, no biopsy here """
    E = biopsy_evaluation_table(filename)
    E.get_graph(G)
    E.complete_eval_and_storage()
    del(E)
    gc.collect()
    return


def Routine_2024(G,fullpathname,target_size = 67000):
    """ Input: Graph G and filename
        and gen-tree stored under filename
        Does a quite extensive evaluation and a certain type of biopsy """
    print("Doing final eval")
    E = biopsy_evaluation_table(fullpathname)
    print(E._filename)
    E.graph = G
    
    E.clonal_eval_2023(biopsy="4_quad_on_same_level",target_size= target_size)
    
    del(E)
    gc.collect()
    return

def small_mutation_eval_routine(G, filename):
    """ Input: Graph G and filename
        Stores information about mutations in excel-sheets
        and active boundary in a .txt """
    E = table_of_evaluation(filename)
    E.get_graph(G)
    E.store_all_mutations()
    del(E)
    gc.collect()
    return

#%%
def genotype_in_biopsy_routine(graphdir):
    """ Adds information about genotypes found in biopsy
        needs graphdir as input, therein valid biopsies and all_mut_file
    """
    EVAL_biopsy.find_all_mutations_in_biopsies(graphdir = graphdir)
    print("Found in Biopsy added, starting eval")
    for M in ["clustering", "all_mutations"]:
        EVAL_biopsy.important_mutations_in_biopsies(graphdir, MioSize = 20,
                                        cutoff_perc = 1.0,
                                        mode = M,
                                        biopsy_sens_perc = 5.0)
        
def all_mut_enhancer(filename):
    """ Input: all_mutations sheat at the given filename (complete path needed)
        adds the information "cells of type" (genotype) to each row """
    # Load old excel-sheet
    table = pd.read_excel(filename)
    
    # Length of common ancestor
    origin = table.at[0, "Genetic Formula"]
    print(origin)
    ANC = len( literal_eval(origin) )
    print(ANC)
    
    # convert formula to name of parent
    parents_name = [ np.nan for i in range(len(table)) ]
    for i in range(len(table)):
        name = table.at[i, "Genetic Formula"]
        if type(name) == float:
            continue
        else:
            # string to tuple
            nametup = list( literal_eval(name) )
            if len(nametup) > ANC:
                parent = str (tuple ( nametup[:-1]  ))
                parents_name[i] = parent
    print(parents_name[:10])
    
    # index of a certain genetic formula
    formula_position = {}
    for i in range(len(table)):
        name = table.at[i, "Genetic Formula"]
        if type(name) == float:
            continue
        else:
            formula_position[name] = i
    
    # index of the parent
    parents_position = [ np.nan for i in range(len(table))]
    for i in range(len(table)):
        parent = parents_name[i]
        if type(parent) == float:
            continue
        else:
            parent_number = formula_position[parent]
            parents_position[i] = parent_number

    # first: cells = branch_weight as before
    def conv(i):
        if type(table.at[i, "Genetic Formula"]) ==float:
            return np.nan
        else:
            return int(table.at[i, "Cells"])
    cells_of_formula = [ conv(i) for i in range(len(table)) ]
    
    # substract branch-weight of child from the parent-weight
    for i in range( len(table) - 1, -1, -1) :
        if type(parents_position[i]) == float:
            continue
        else:
            parent_id = parents_position[i]
            cells_of_formula[parent_id] -= int(table.at[i, "Cells"])
    S = [ i for i in cells_of_formula if type(i) == int]
    print("New total sum:", sum(S))
    
    # add new information to table and store
    table["Cells of formula"] = cells_of_formula
    table.to_excel(filename,index=False)
    return table

