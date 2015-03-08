import sys
import params # Parameters for the neurons
import numpy as np
import itertools
import shutil
import os
import generateFig4 as genF4  # Main NEST code to run simulations
import generateFig4_lossy as genF4_los # NEST code with lossy synapse
import Sim_raster as S_r # Plots the results

pars = params.get_parameters() # 	
postfix = "Fig4_"
#genF4_los.run(postfix,pars) # Uncomment this line and comment the line below for using the lossy synapse version to generate DTT
genF4.run(postfix,pars) 
S_r.drawFig4()
			
