parbrain
========
Parbrain stands for "parallel brain". This code generates a binary tree, which globally couples multiple neurovascular units (NVUs). The binary tree represents a vascular tree, which perfuses the cerebral tissue. It is laid out in such a way that it spatially fills a 2-dimensional tissue slice consisting of NVU blocks. The code computes the blood flow in the branches of the tree and the pressure in the nodes. The NVUs are located at the leaves of the tree and can (through smooth muscle activation/relaxation) control the vessel radii. 


Inputs to the Model
===================
There are multiple inputs to the model coming from different "ends":

From the vasculature:
---------------------
* pressure in the root vessel
* PLC production flux 

From the tissue:
----------------
* neuronal K+ release/Na+ uptake - Ostby input
* neuronal Glu release -> activates NO production


How to Compile
==============
Configure the project with CMake and specify the *out-of-source* build directory.

After the the project has been configured, it can be opened and compiled with the target IDE, or in the case of Unix Makefile configuration simply run make in the build directory. 

Requirements:
-------------
* VTK
* CXsparse
* MPI


How to Run
==========
`mpirun -np <number_of_processors> <program> <number_of_levels> <subtree_size>`

Where `<program>` is "parBrainSim" located in the build directory.
The subtree size is usually 3. This outputs a directory with multiple data files.
	

How to View Binary Data
================
`python <script> <data_directory> <number_of_levels>`

Where `<script>` is "printData.py" located in the util subdirectory.
Gives you info, coordinates of the tissue block, state variables, flow and pressure for all timesteps.


Viewing in Paraview
==================

`<program> <data_directory>`

Where `<program>` is "ConvertBinToVtu" located in the build directory and <data_directory> is the directory generated by parBrainSim.

This will generate sets of vtu/vtp files in <data_directory> for the tissue blocks, H tree lines and tubes for each timestep. Open in Paraview for visualisation.


Required modules Power7:
------------------------
module load suitesparse/4.2.1_adv

module load vtk
