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
Simply run *make* in the main directory.

How to Run
==========
`mpirun -np <number_of_processors> ./simulate <number_of_levels> <subtree_size>`





Required modules Power7:
------------------------
module load suitesparse/4.2.1_adv
module load vtk
