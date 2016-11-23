
# Background

This repository will contain the notes and code for investigating a quantum Two-Level System (TLS), coupled strongly to a phonon bath and coupled weakly to the ambient, incoherent electromagnetic field.
There is:
- a central reaction coordinate Liouvillian builder called UD_liouv.py. It creates the master equation for the strongly coupled vibrations.
- a module which deals with all of the different types of incoherent driving called driving_liouv.py.
- a module with several different types of checks. The main feature is one which attempts to determine which systems are likely to be susceptible to non-secular effects.
- a plotting module which takes in the other two and plots graphs of the dynamics, coherences and emission spectra. This could be extended into an IPython notebook as well.
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
- In driving_liouv the nrwa Liouvillian does not work. I think it won't converge for a flat spectrum.
- It seems like when you have a small RC frequency the coherences diverge and the dynamics exhibit weird revivals and chaos.

# To do:
- Try to find the "smoking gun" for needing to use the non-secular treatment. E.g does the structure of a manifold and occupation number dictate exactly when secular works?
- Locate each manifold in the NE_checking.plot_manifold() and plot them a different colour.
- I still haven't accounted for the resonance effects in the Secular master equation (FRET).
- Choose a better optical field spectral density which allows the Principal parts to converge for FRET and renormalisation terms.
- Put all of the globally run code in ME_plotting into an IPython notebook with a full discussion of what's going on.
- Save data of a few relevant parameter regimes to file and load them into a pandas dataframe for easy analysis using matplotlib and seaborn.
- Maybe create a TLS class, which has methods to analyse different parts of the dynamics/steady states.
