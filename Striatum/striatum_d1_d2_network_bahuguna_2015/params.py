#  Author: Jyotika Bahuguna: -j.bahuguna@fz-juelich.de
from NeuroTools.parameters import ParameterSet
from NeuroTools.parameters import ParameterRange
from NeuroTools.parameters import ParameterTable
from NeuroTools.parameters import ParameterSpace
import NeuroTools.signals as signal

import numpy,shelve,pylab,os


def get_parameters():

	p = ParameterSpace({})
	# Parameters for neuronal features
	p.outpath = '.'
	p.vm = -80.
	
	p.th1 = -45.
	p.th2 = -54.
	p.th3 = -45.
	p.tau_synE1 = 0.3
	p.tau_synE2 = 0.3
	p.tau_synI1 = 2.
	p.tau_synI2 = 2.
	p.E_ex = 0.
	p.E_in1 = -64.
	p.E_in2 = -76.
	p.ie = 0.
	p.cm1 = 192.   	# For MSN (Gertler 2008)
	p.gL11 = 8. # From Gertler, D1 and D2 1/(124.4Mohm) and 1/(154.83Mohm)
	p.gL12 = 6.
	p.cm = 200.		# For MSN (Wolf 2005)	
	p.cm_fsi = 500.		# For MSN (Wolf 2005)	
	p.gL1 = 12.5		# For MSN (Wolf 2005)	
	p.cm2 = 157.
	p.gL2 = 25.
	p.tref = 2.
	p.vi1low = -80.
	p.vi2low = -80.
	p.vi3low = -80.
	p.vi1hi = -45.
	p.vi2hi = -54.
	p.vi3hi = -45.
	
	# Parameters for running
	p.timestep = 0.1
	p.min_delay = 0.1
	p.max_delay = 50.
	p.runtime = 500.
	
	p.num01 = 150 # Neuron population in the cortex that recieve correlated input 
	p.num02 = 1
	p.num1 = 4000 # Pair of neurons in MSN ( inhibitory ) which recieve input from cortex  The ratio of cortex::MSn is 10:1
	p.numAll = p.num1/2
	p.numFSI = 80
	p.p_copy = 0.03
	p.nc21 = 10
	p.prob11 = 0.23
	p.delay11 = 1.
	p.delay12 = 4.
	p.delay21 = 1.
	p.delay22 = 0.5
	p.j01 = 3.1
	p.j02 = 0.55
	minInh = -0.5
		
	p.jd1d1 = minInh
	p.jd2d2 = minInh*2.0
	p.jd1d2 = minInh*2.0*1.21
	p.jd2d1 = minInh*2.0*1.21*1.32
	
	p.jfsi = -2.6
	
	p.j21 = -1.0  # New value , IPSP Koos1999
	

	p.Rate = numpy.arange(50.,4560.,500.) 

	return p	
	
	
