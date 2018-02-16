These are the set of scripts used to run umbrella sampling with OpenMM.

**Warning - they have passed through a few hands, are some what Frankenstein-ish and not necessarily well annotated.

restrOpenMM_GPU.py
This script used to take in a series of mol2 files and corresponding inpcrd and prmtop files starting from different angles of restraint and then run simulations restraining around the starting angle for each. I've set it to load the same prmtop and inpcrd file for every angle and then loop over a range of angles of restraint.  It does however, still look for a starting mol2 directory (should be addressed at some point).  This script is executed in the slurm_gpu.sh script with line: python restrOpenMM_GPU.py -i mol2files/ -f gaff -g mol2files/ > minGAFF.out

trajCompnent.tcl
To run this, copy script to same directory where GAFF/ is located and enter this line of code: for k in {0..350..10}; do echo $k; cd GAFF/$k; vmd -dispdev none -e ../../trajComponent.tcl -args 5 ../../vacDivaline.prmtop nvt01.dcd; cd ../../; done
This script computes the dihedral angles that were generated during the simulations (you can modify which atom indices that it computes the dihedral for, if you'd like).  The number of equilibration frames to be discarded can be specified following -args (it's 5 in the example above)

overlap.py
This script generates a graph of the dihedral angles for each simulation window and puts it in the GAFF/ directory

centers.dat
This file needs to be manually updated with the actual angles of restraint and corresponding force constatns.  It is used by the umbrella-sampling1.py script.

umbrella-sampling1.py
This script computes the PMF from your dihedral distributions. On line 24, you can modify the number of bins used to generate the PMF data.

Steps:
1) Set desired force restraint(s) on lines 241-250 of restrOpenMM_GPU.py as well as the range of restraint angles (line 195). Modify tthe path to the system input files (191 and 192).
2) Submit for sampling using gpu_umbrella.sh script
3) Once job completes, you will have a folder called "GAFF" containing other folders named by the angles you restrained around.
4) You should also have an output log file minGAFF.out which gives information about whether the simulations failed or not, and what force constants were used for each angle.
5) Navigate to folder containing GAFF/ directory (likely where you submitted job) and run this line of code to compute the dihedral angles for each simulation: 
    To run this, enter this line of code: for k in {0..350..10}; do echo $k; cd GAFF/$k; vmd -dispdev none -e ../../trajComponent.tcl -args 5 ../../vacDivaline.prmtop nvt01.dcd; cd ../../; done
6) Next run overlapy.py.  It will generate an overlap plot of all the sampled dihedral angles.
7) Update the centers.dat file to contain the appropriate restraint angles and force constants
8) Next run umbrella-sampling1.py and save the output to a log file: ie - python umbrella-sampling1.py > 15bin_PMF.log
8) Plot the bin (x-values) and f (y-values) at the end of the log file generated in step 8.


