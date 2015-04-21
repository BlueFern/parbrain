Visualisation for parbrain
==========================
Converts binary output files to .vtu / .vtp files that can be openend with paraView.

on power7:
module load vtk

LDFLAGS="-Wl,-rpath=$VTKLIB -Wl,-as-needed" cmake -DCMAKE_BUILD_TYPE=Release ../binToVtu

make

./ConvertBinToVtu <folder name> <final time>

