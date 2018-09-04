parbrain - NVU 2.1 with extracellular diffusion, astrocytic gap junctions and cortical curvature
========
Parbrain stands for "parallel brain". This code generates a binary tree, which globally couples multiple neurovascular units (NVUs). The binary tree represents a vascular tree, which perfuses the cerebral tissue. It is laid out in such a way that it spatially fills a 2-dimensional tissue slice consisting of NVU blocks. The code computes the blood flow in the branches of the tree and the pressure in the nodes. The NVUs are located at the leaves of the tree and can (through smooth muscle activation/relaxation) control the vessel radii. 

NVU 2.1 contains the following pathways:
* Neuronal dynamics stimulated by a current input to the soma
* K+ pathway via K+ release from the neuron into the SC
* Nitric oxide pathway via glutamate release from the neuron into the SC
* Astrocytic calcium pathway via glutamate release from the neuron into the SC
* TRPV4 calcium channel on the astrocytic endfoot

The NVU model contains extracellular electrodiffusion of K+ and Na+ and an astrocytic K+ gap junction network throughout the tissue slice. Both normal neurovascular coupling and pathological conditions (cortical spreading depression) can be simulated. 

The model can also simulate the spatially varied Gaussian curvature of the cortex which affects the extracellular diffusion rate. This is an optional input and can be turned off in the constants.h file.

For further information on the NVU model refer to the following papers:
* Mathias, E., Kenny, A., Plank, M. J., & David, T. (2018). Integrated models of neurovascular coupling and BOLD signals: Responses for varying neural activations. NeuroImage, 174(March), 69–86. http://doi.org/10.1016/j.neuroimage.2018.03.010
* Kenny, A., Plank, M. J., & David, T. (2018). The role of astrocytic calcium and TRPV4 channels in neurovascular coupling. Journal of Computational Neuroscience, 44(1), 97–114. http://doi.org/10.1007/s10827-017-0671-7
* Dormanns, K., Brown, R. G. G., & David, T. (2016). The role of nitric oxide in neurovascular coupling. Journal of Theoretical Biology, 394, 1–17. http://doi.org/10.1016/j.jtbi.2016.01.009
* Dormanns, K., van Disseldorp, E. M. J., Brown, R. G., & David, T. (2015). Neurovascular coupling and the influence of luminal agonists via the endothelium. Journal of Theoretical Biology, 364, 49–70. http://doi.org/10.1016/j.jtbi.2014.08.029
* Farr, H., & David, T. (2011). Models of neurovascular coupling via potassium and EET signalling. Journal of Theoretical Biology, 286(1), 13–23. http://doi.org/10.1016/j.jtbi.2011.07.006

And for further information on the parallel implementation refer to the following papers:
* Kenny, A., Zakkaroff, C., Plank, M. J., & David, T. (2018). Massively parallel simulations of neurovascular coupling with extracellular diffusion. Journal of Computational Science, 24, 116–124. http://doi.org/10.1016/j.jocs.2017.07.001
* Dormanns, K., Brown, R. G., & David, T. (2015). Neurovascular coupling: a parallel implementation. Frontiers in Computational Neuroscience, 9, 109. http://doi.org/10.3389/fncom.2015.00109

Inputs to the Model
===================
There are multiple inputs to the model coming from different "ends":

From the tissue:
----------------
A current input to the neuron (soma) stimulates the following:
* neuronal K+ release/Na+ uptake
* neuronal Glu release -> activates NO production and astrocytic calcium release

From the vasculature:
---------------------
* pressure in the root vessel
* PLC production flux 


How to Compile
==============
To compile using Makefile configuration simply run `make` in the build directory (*must be done separately for both the main program in the main directory and the ConvertToVtu file in the util/bintovtu directory*). 

Alternatively, configure the project with CMake and specify the *out-of-source* build directory.
After the the project has been configured, it can be opened and compiled with the target IDE. 

Requirements:
-------------
* VTK (for bin_to_vtu script which outputs vtu files for viewing in Paraview)
* CXsparse (for sparse matrices)
* MPI (for core to core communication)

How to Run
==========
`mpirun -np <number_of_processors> <directory>/parBrainSim <number_of_levels> <final_time> <output per second> <curvature map (if using)>`

Where the program "parBrainSim" is located in the build directory. `<number_of_levels>` is the vascular tree size ,`<final_time>` is the final output time, and `<output per second>` is the number of times per second the data is saved. These are optional arguments. By default these are set in the constants.h file. If the parameter `CURVATURE_SWITCH` is set to 1 in constants.h then the name of the curvature map (in either txt or csv format), <curvature map (if using)>, is required.

This program outputs a directory `<data_directory>` with multiple data files (info.dat, tissueBlocks.dat, pressure.dat and flow.dat).

Modifying the Model
===================
Any changes to model parameters can be made in the constants.h file in the src directory, then re-compile the code.	

How to View Binary Data
================
`python <directory>/printData.py <data_directory>`

Where "printData.py" is located in the util subdirectory.
Prints out to the terminal: info, coordinates of the tissue block, state variables, flow and pressure for all timesteps.


Viewing in Paraview
==================
The program to convert the raw data into vtu files for paraview is optimised with OpenMP and is made with a separate makefile in the util/binToVtu directory. A bash wrapper is used to run the program with OpenMP. 

`<directory>/RunConvertBinToVtu.bash <Num of threads> [-c <CPU affinity mask> (optional)] <Data directory> <Final time> <Output per sec>`

Where "RunConvertBinToVtu.bash" is located in the util/binToVtu subdirectory. `<Num of threads>` is the number of threads dependent on the hardware (14 for brats01 machine), `<CPU affinity mask>` is an optional argument ("0-13" for brats01 machine), `<data_directory>` is the directory containing specific simulation data generated by parBrainSim, `<final_time>` is the final output time, and `<output per second>` is the number of times per second the data is saved to file.

This will generate sets of vtu/vtp files in `<data_directory>` for the tissue blocks, H tree lines and tubes for each timestep. Open in Paraview for visualisation.


Required modules Power7:
------------------------
module load suitesparse/4.2.1_adv

module load vtk
