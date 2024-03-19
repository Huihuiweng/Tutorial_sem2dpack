
This code is used to simulate 2D and 2.5D dynamic rupture simulations, with slip-weakening, time-weakening, rate-and-state, and thermal pressurization friction laws. This tutorial aims to guide you to simulate the models and post process the results.

The structure of the tutorial document:
code              -- contains the original sem2dpack code
simulations       -- contains the model that can be run directly'
scripts_templates -- contains matlab and python scripts
post-process      -- choose the script you like and save it here


Current directory is defined as ${workdir}

Installation:
1. cd ${workdir}/code/sem2dpack/SRC
2. vi Makefile   ### Modify the executable dir and the type of compiler in Makefile
3. make

Running model:
1. cd ${workdir}/simulations
2. vi Par.inp    ### Modify the model parameters, refer to the manual of sem2dpack at the end of this file
3. ${path_to_bin}/sem2dsolve Par.inp

Post process results:
1. cd ${workdir}/post-process

for python version:
2. cp ../scripts_templates/python/* .
3. python Example_figure.py
for matlab version:
2. cp ../scripts_templates/matlab/* .  
3. open Example_figure.m in MATLAB and run the script directly.


# Manual for running sem2dpack see ${workdir}code/sem2dpack/doc/doc.txt 
