#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 10:13:03 2018

@author: floriankreten
"""

""" Testet with Python 3.7.13 and plotly 5.9.0 """


""" Simplified parameter setup for visualization:
        Small simulation with high rate of driver mutations
"""

# First step: define behavior via the rates-file
# All parameters and rates are defined/can be chosen in the INPUT_rates
# You should start a new instance/reset the kernel after changing the parameters
import INPUT_rates

# fixate the simulation-framework after all parameters are set:
INPUT_rates._set_rates()

#%% Small colourful example without further usage, good for module-testing :)
import EVAL_over_time_plot_plotly
EVAL_over_time_plot_plotly.sim_and_plot_with_size_frames()


#%% Simulation with 3D plots at defined sizes
# sizes should not exceed 100.000 vertices for this visualization
# bigger plots use a "lazy" scheme and show in fact only subsets
import PROCESS_process

# need initial estimate of the spatial size (as radius) for storing the spatial coordinates in a grid
# 800 is enough for 20mln vertices
G = PROCESS_process.init_simulation_graph(300)

final_size = 20000
# can stop/start simulation for any given graph
G = PROCESS_process.continue_simulation(G, final_size)

#%% Data storage and reading is handled via EVAL_data_output
# Diverse evaluation schemes are bundled into
import EVAL_data_output
import EVAL_single_graph_evaluation

# store graph at specific position for creating an evaluation-directory
eval_dir = "/mnt/Prostate_Cancer_Partition/DRIVE_Prostate_2024_Sim_Files/Sim_Files_Stored/Clean_For_Publication/BLURP"
EVAL_data_output.save_graph_sql(G, eval_dir, fullpathname = True)

# load graph (lazy version for evaluation by default)
# G = EVAL_data_output.get_graph_sql(eval_dir, fullpathname = True)

# Load the graph and the eval-dir into an evaluation-object
E = EVAL_single_graph_evaluation.biopsy_evaluation_table(eval_dir)
E.get_graph(G)
    
# store phylogenetic information into a tree-db (for internal use) and an excel-file for readability
E.store_all_mutations(fullpathname= True)
# load tree into the eval-object
E.get_gen_tree_from_directory() 

# Biopsy-Eval, and plot of phylogenetic tree + spatial genotype distribution
# Per default, mutations with ccf<5% are ignored  
E.clonal_eval_2023(biopsy="4_quad_on_same_level", target_size=50)
E.spatial_trait_and_anc_plot(title = "Clones_Spatial", show = True, fullpathname = True)
