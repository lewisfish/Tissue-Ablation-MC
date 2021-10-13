#        Tissue-Ablation-MC

[![DOI](https://zenodo.org/badge/86356329.svg)](https://zenodo.org/badge/latestdoi/86356329)


This repository contains the source codes for the 3D grid Cartesian monte carlo radiation transfer codes developed as part of my phd.
The codes simulate the transfer of light through tissue/phantoms and calculates the heat induced due to the laser. Heats effect on tissue is also modelled with tissue ablation etc.

[Paper in Lasers in Surgery and Medicine 2020](https://doi.org/10.1002/LSM.23335)

[Thesis chapter on how code operates](https://github.com/lewisfish/Tissue-Ablation-MC/blob/master/chapter-finished.pdf)


![Heat Sim](https://github.com/lewisfish/Tissue-Ablation-MC/raw/master/Heat_3D.gif)


#### The FORTRAN source files are:
            
            ran2.f               random number generator
            constants.f90        contains the various constants used in the simulation
            photon_vars.f90      contains the various photons variables
            iarray.f90           contains the arrays variables names
            opt_prop.f90         tracks optical properties of the system
            subs.f90             inits arrays and directory paths 
            ch_opt.f90           sets rhokap and albedo arrays with correct values
            gridset.f90          sets up the various grids
            inttau2.f90          generates a tau value and integrates through the grid
            sourceph.f90         photon source subroutine
            writer.f90           writes out results
            3dFD.f90             Heat simulation
            mcpolar.f90          main program. calls heat sim and does MC sim

#### Input parameters are in:

	input.params

##### The file that compiles the code and creates the executable file 'mcgrid' is:

   Only been tested on linux so far. Works with intel and gfortran compilers.
   Also been tested on computing clusters. http://www-solar.mcs.st-and.ac.uk/~herbert/cluster/
	
	install.sh
	
	This can be run by ./install.sh
	May have to change permissions first in order to execute the script.
	This can be done by using sudo chmod +755 instsll.sh on linux
  
  ./install.sh -m 
  
   just complies the code
  
  ./install -n (number of cores)
  
   complies and runs the code on (number of cores). If left with nothing the default is 1 core.

#### To do/implement
   
   - Speed up if possible
   - More...
