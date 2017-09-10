# conformism
Conformism model based on an Montecarlo (metropolis) simulation with a 2D-Ising model

Author : M. Le Verge

Object : this is an algorithm dedicated to the study of the social conformism on individuals of a society, based on a modified 2D Ising Model, and using Montecarlo method with Metropolis algorithm.

Description of files :
- conformism.py : main script code. Two classes are defined : conformism2D and conformismfull. The first is a square lattice montecarlo script, the second one is based on a fully connected graph with montecarlo script.
- run_conformism.py : running script, which call conformism.py. Running run_conformism.py will produces graphs and plots of Energy, Magnetization, Heat Capacity and Susceptibility of the system, vs the desired parameters.

Sources :
