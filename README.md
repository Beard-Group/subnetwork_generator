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
