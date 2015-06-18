# subnetwork_generator
Beard Lab Subnetwork Generator Program, February 2015

Clone the repository in a new directory

      git clone git@github.com:Beard-Group/subnetwork_generator.git
      
Build the program by navigating to subnet_generator/src/ and invoking make

      make
      
Run the program by running the automatic test script, "profile.sh".
This script is set up to run the program, and plots the predicted time course
profiles along with the actual data. You can tell the program how many time 
profiles to plot by giving a numerical argument to profile.sh on the command line.
The default is 3.
For example:

      ./profile.sh 4
    

You may need to modify paths in the Makefile to be able to compile the program.


---------------------------------------------------------------------------------------------
Inputs are read through src/io/* files. The input files being (for number of genes N_gene):

1. t_d: time points where data is available
2. x_d: value of expression at these timepoints. It has N_gene rows and size(t_d) columns.
3. n_ka: vector file with number of activators for each gene. Size: N_gene x 1.
4. n_kd: vector file with number of inhibitors for each gene. Size: N_gene x 1.
5. kavec: index of gene which acts as an activator. Index starts at 0 for the sample input files. Size is sum of all elements of n_ka vector.
6. kavec: index of gene which acts as an inhibitor. Index starts at 0 for the sample input files.Size is sum of all elements of n_kd vector.
7. kaval: Value of activating effect for the corresponding element in kavec. Size is sum of all elements of n_ka vector.
8. kdval: Value of inhibiting effect for the corresponding element in kdvec. Size is sum of all elements of n_kd vector.
-----------------------------------------------------------------------------------------------
