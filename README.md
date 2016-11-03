This directory will be for dealing with a TLS strongly coupled to vibrations.
There is:
- a central reaction coordinate Liouvillian builder called UD_liouv.py.
- a module which deals with all of the different types of incoherent driving called driving_liov.py
- a module with several different types of checks. Convergence, non-secularity etc. This is going to be left for a while.
- a plotting module which takes in the other two and plots beautiful graphs. This could be extended into an ipython notebook as well.


TODO:
The spectral density in the nrwa case seems to be incompatible with the a flat spectrum. The coth term diverges in this case. I think I will need to have an ohmic spectral density with some cutoff.
