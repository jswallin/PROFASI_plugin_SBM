# PROFASI_plugin_SBM

The software in this repository is meant to be used with the PROFASI software package for protein simulations. For information on PROFASI:

http://cbbp.thep.lu.se/activities/profasi/
http://people.binf.ku.dk/wb/profasi/

This repository uses the PROFASI plugin facility. It implements an structure-based potential energy term, as described in the following publication: 

Stefan Wallin and Bahman Seifi, "The C-terminal domain of transcription factor RfaH: Folding, fold switching and energy landscape", Biopolymers (2021). 

Installation: 

Place the directory BiasEnergy/ in the plugins/ directory within the PROFASI source tree. Re-compile PROFASI, as described in PROFASI documentation. 

Use: 

As an example of usage, see the plugin_test/ directory. 
