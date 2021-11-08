# DFT-Clustering
This code allows one to take a DFT calculation using the computatational platform StoBe of NEXAFS spectra. The result is a set of peaks that make up a Building Block model using a tensor based formalism that can be used to carry out simultaneous fits on angle resolved NEXAFS for the extraction of the molecular tilt angle.

The algorithm works by defining the following parameters:
1. An energy cutoff that defines what is the maximum DFT transition energy to consider
2. An oscillator strength threshold that filters transitions that do not have a sufficiently high intensity from subsequent steps
3. A peak overlap threshold that determines whether two transitions can be clustered together depending on the overlap area between them. 

The code runs on IGOR however a Python implementation may be developed in the near future. Also, the code takes StoBe output files as input, however, as long as the computational platform provides transition energies, transition intensities and the components of the transition dipole moment, then the loading function can be modified to accomodate other platforms. 
