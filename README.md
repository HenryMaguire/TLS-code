
# Background

This repository will contain the notes and code for investigating a quantum Two-Level System (TLS), coupled strongly to a phonon bath and coupled weakly to the ambient, incoherent electromagnetic field.
There is:
- a central reaction coordinate Liouvillian builder called UD_liouv.py. It creates the master equation for the strongly coupled vibrations.
- a module which deals with all of the different types of incoherent driving called driving_liouv.py.
- a module called superdriving_liouv.py which does the same as driving_liouv.py but for a superohmic spectral density (this is a cheap hack)
- a module with several different types of checks. The main feature is one which attempts to determine which systems are likely to be susceptible to non-secular effects.
- a plotting module which takes in the other two and plots graphs of the dynamics, coherences and emission spectra. This could be extended into an IPython notebook as well.
- a mathematica notebook for calculating the Franck-Condon overlap factors of the various vibrational states
- a module for calculating the exact dynamics of the Independent-Boson Model (currently does not work).
- a directory with all of the accompanying notes and figures for the investigation, read Vibronic_incoherent_notes.pdf to get some more physical insight.

# Requirements

All the python files are written in Python 2.7. The modules will need to be at least:
- Qutip >=3.1.0
- Numpy
- Scipy
- Matplotlib

# Getting started
- Clone the repo and install the python dependencies
- Open ME_plotting.py, choose some parameters and run it.
- Check in the Notes/Images folder for default plots of dynamics and spectra, or alternatively use the data files in DATA to plot your own.

# Bugs


# To do:
- Make sure all the code is sideways compatible after the two large merges. This means making sure the API is consistent between all the function calls in the notebooks.
- Have two notebooks to start with. One for underdamped and one for overdamped spectral densities. Focus initially on the underdamped one.
- Perhaps have a third notebook for the rigorous RWA comparisons. Clear the notebook these currently exist in of this work.
- Move all general code from the jupyter notebooks into the scripts. Also move all legacy shitcode into some legacy file.
- Make it so there's no OD/UD prefixes on the code files but one codebase is used for each case with some default parameter passed instead.

