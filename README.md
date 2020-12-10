# PROFASI_plugin_SBM

The software in this repository is meant to be used with the PROFASI 
software package for protein simulations. 

For information on PROFASI:

http://cbbp.thep.lu.se/activities/profasi/
http://people.binf.ku.dk/wb/profasi/

The code in this repository uses the PROFASI plugin facility and 
implements a structure-based potential energy term, which can be 
used to modify the base energy function implemented in PROFASI. 

The structure-based potential energy term is described in the 
following publication: 

Stefan Wallin and Bahman Seifi, "The C-terminal domain of 
transcription factor RfaH: Folding, fold switching and energy 
landscape", Biopolymers (2021). 

Installation: 

Place the directory BiasEnergy/ in the plugins/ directory within 
the PROFASI source tree. Re-compile PROFASI, as described in 
the general documentation. 

Usage: 

As an example of usage, see the plugin_test/ directory. 

