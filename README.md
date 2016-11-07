This repository will contain the notes and code for investigating a Quantum Two-Level System, coupled strongly to phonon bath and coupled weakly to the ambient, incoherent electromagnetic field.  The main
There is:
- a central reaction coordinate Liouvillian builder called UD_liouv.py. It creates the master equation for the strongly coupled vibrations.
- a module which deals with all of the different types of incoherent driving called driving_liov.py.
- a module with several different types of checks. Convergence, non-secularity etc.
- a plotting module which takes in the other two and plots beautiful graphs. This could be extended into an ipython notebook as well.
- a directory with all of the accompanying notes for the investigation.

TODO:
- The spectral density in the nrwa case seems to be incompatible with the a flat spectrum. The coth term diverges in this case. I think I will need to have an ohmic spectral density with some cutoff.
