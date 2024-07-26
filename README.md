Tested with Python 3.8.16 and plotly 5.9.0.
(complete conda_env.yml in the repository)

This repository contains a pure python implementation of our prostate-cancer model.
This stochastic module can be used to grow three-dimensional tree-structures, as well as running a contact process on them - biologically speaking: placing mutations and letting the genotypes compete within the network.

On top of that, modules for all sorts of easy-to-use storage, loading, graphical and Excel-output are given.

For a ready-to-use example that explains the main steps, including storage and graphical output, follow the lines along the prepared main.py.


------------------------------------------------------------------------------------


IMPORTANT: Include INPUT_rates first into the main and define a type-of-rate (see example).
	   Changing the rate-functions in a running session can be done, but is not recommended.
	   Better store the graph and use a new instance.


For the graphs with 2*10⁷ nodes we used in the publication and with the given mutation rates,
a single sim can take up to 30GB ram and depending on the CPU up to ~three hours.
We recommend to test smaller examples first.


The modules are sorted into 3 different groups:

1) PROCESS-modules are relevant for the simulations and should not be touched
2) EVAL-modules concern graphical output, storage and evaluations
3) 	Via INPUT_rates, the entire behavior of a simulation can be defined (rates, types of events).
	The main.py should be used for all user-defined scripting

The graph-class defined in PROCESS_structures is the central object.
It is a collection of vertices and some structures/information.
A vertex is a structure that contains all locally relevant information:
    position (in R^3), neighbors connected via edges, current genotype, rates of process at v

The graph changes via a sequence of random events. At the moment, three events are implemented:
PROCESS_growth, PROCESS_mutation, PROCESS_competition (PROCESS_radial_growth is not tested properly).

The Gillespie-simulation is implemented in PROCESS_process, the underlying calculations are done in the PROCESS_event_manager
I recommend to leave these modules at peace :)

All evaluation-modules are bundled into the EVAL_single_graph_evaluation module.
The central object here is called "table_of_evaluation". After loading a graph into this table,
you can evaluate it in all possible senses and produce some graphical output. See again the main for an example workflow


The data IN/OUT is handled via the EVAL_data_output module. Graphs and ancestral trees can be stored
as SQL-databases. For working with stored data again, use the get_graph and get_tree functions from EVAL_data_output.


