Author: Jyotika Bahuguna - j.bahuguna@fz-juelich.de


This is the code to show the presence of DTT in a spiking neural network for Scenario II.
=========================================================================================

script_Fig4.py : This is the main script that calls two scripts. 1) To generate the firing rates (generateFig4.py) 2) To plot the firing rates
		with rasters as shown in Fig4.

generateFig4.py: This is the NEST simulation script that runs the simulation of striatal network for different cortical rates. NEST version used:2.2.2


Sim_raster.py: This plots the firing rates as shown in Fig4. Python version used: 2.7 and correspondingly compatible matplotlib, scipy ,numpy 
etc libraries (Lower versions of python or matplotlib might give errors for some functions ,eg, subplot2grid.)

Note: There is another version of generateFig4.py - "generateFig4_lossy.py". This was a a customized synapse type - "lossy" designed to 
efficiently generate Between and Within correlations (Refer to Methods -"Generation of B and W") by Alejandro Bujan (afbujan@gmail.com) and 
Susanne Kunkel (s.kunkel@fz-juelich.de). Unfortunately this is not yet a part of the main NEST code. In future, if this becomes a part of NEST,
generateFig4_lossy can be used to generate Figure 4 by uncommenting a line in script_Fig4.py. This version of code is also required to generate
figures with different values of B and W correlations.





