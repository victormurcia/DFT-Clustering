# DFT-Clustering
This program takes the output of DFT calculations carried out in the the computational platform StoBe for the simulation of NEXAFS to generate an optical tensor model derived from first principles. The result is a set of peaks that make up an optical model using a tensor based formalism that can be used to carry out simultaneous fits on angle resolved NEXAFS for the extraction of the molecular tilt angle. These peaks can then be used to obtain optical constants that could be subsequently used in the analysis of Resonant Soft X-Ray Scattering (R-SoXS) and Resonant X-Ray Reflectivity (XRR).

<p align="center">
  <img src="images/ovps.png" />
</p>

The algorithm works by defining the following parameters:
1. An energy cutoff that defines what is the maximum DFT transition energy to consider
2. An oscillator strength threshold that filters transitions that do not have a sufficiently high intensity from subsequent steps
3. A peak overlap threshold that determines whether two transitions can be clustered together depending on the overlap area between them. 

These peak overlaps can then be subsequently used to generate a more compact set of peaks (i.e. transition clusters) that are representative of all the transitions initially calculated by the TP-DFT. The transition clusters are then combined with angle-resolved NEXAFS measurements in order to generate a quantitatively accurate optical model derived from first principle calculations.

<p align="center">
  <img src="images/znpc bb fits for xrr both.png" width="500" height="500">
</p>

Finally, the transition clusters that comprise the optical model can be used to identify the chemical, energetic and orientational character of the various NEXAFS features in addition to allow these NEXAFS features to be associated with specific MOs calculated from the TP-DFT.

<p align="center">
  <img src="images/dft bb to mo cl8.png"  width="500" height="500">
</p>

The code runs on IGOR 8, however a Python implementation may be developed in the future. Also, the code takes StoBe output files as input, however, as long as the computational platform provides transition energies, transition intensities and the components of the transition dipole moment, then the loading function can be modified to accomodate other platforms. 

<p align="center">
  <img src="images/gui.png" />
</p>

The accompanying python files are there to:
1. Facilitate the procedural generation of the .run files for a Transition Potential calculation carried out in StoBe 
2. Extract the Mulliken Population Analysis from the StoBe output files that can be subsequently loaded into IGOR to aid in the chemical characterization of NEXAFS transitions. 
 
 Setup
1. Navigate to the following directory: Documents > Wavemetrics > Igor Pro 8 User Files
2. Place the file clusteringPanel v1.ipf in the folder named "Igor Procedures"
3. Navigate to the folder Documents > Wavemetrics > Igor Pro 8 User Files > User Procedures  
4. Place the contents of the folder "DFT_Clustering" inside the "User Procedures" directory
5. Open an IGOR instance. There should be a tab titled Macros there. Within the dropdown menu in Macros there should be an option titled "Clustering Algorithm" which will load in the control panel for the algorithm.

Bug reporting/Help
For any concerns email me at victor.murcia@wsu.edu
