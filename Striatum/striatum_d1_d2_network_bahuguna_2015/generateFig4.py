#  Author: Jyotika Bahuguna: -j.bahuguna@fz-juelich.de

import numpy as np
import time
import pylab as pl
import matplotlib.cm as cm
import matplotlib
import scipy
import params 
import scipy.spatial.distance as dist
import gc
from matplotlib.pyplot import *
import pdb


matplotlib.rcParams['legend.fontsize'] = 8
p = params.get_parameters() 

# These will be the dictionaries which will be saved as a .pickle files. 
spike_senders = dict()
spike_times = dict()
rates = dict()

rates["d1"] =[]
rates["d2"] =[]
rates["fsi"]=[]
rates["ctx"]=[]

spike_senders["d1"]=[]
spike_senders["d2"]=[]
spike_senders["fsi"]=[]
spike_senders["ctx"]=[]

spike_times["d1"]=[]
spike_times["d2"]=[]
spike_times["fsi"]=[]
spike_times["ctx"]=[]


def run(prefix,pars):		
	w1 = 3.72  # Jc1
	w2 = 3.0   # Jc2
	w3 = 4.5   # Jfsictx
	np.set_printoptions(threshold='nan')			
	import nest

	for rate in pars['Rate']:
		nest.ResetKernel() # This is needed to destroy the network and run a fresh version for every "rate"
	
		nest.SetKernelStatus({"resolution":p.timestep,"overwrite_files": True,"local_num_threads":8 }) 
		# D1 population
		std1 = nest.Create('iaf_cond_alpha',p.numAll,params={"V_reset":p.vm, "V_th": p.th1, "tau_syn_ex":p.tau_synE1, "tau_syn_in":p.tau_synI1, "E_L":p.vm, "E_ex":p.E_ex, "E_in":p.E_in1, "I_e":p.ie, "C_m":p.cm, "g_L":p.gL1, "t_ref":p.tref})
		# D2 population
		std2 = nest.Create('iaf_cond_alpha',p.numAll,params={"V_reset":p.vm, "V_th": p.th1, "tau_syn_ex":p.tau_synE1, "tau_syn_in":p.tau_synI1, "E_L":p.vm, "E_ex":p.E_ex, "E_in":p.E_in1, "I_e":p.ie, "C_m":p.cm, "g_L":p.gL1, "t_ref":p.tref})
		# FSI population
		stfsi = nest.Create('iaf_cond_alpha',int(p.numFSI),params={"V_reset":p.vm, "V_th": p.th2, "tau_syn_ex":p.tau_synE1, "tau_syn_in":p.tau_synI1, "E_L":p.vm, "E_ex":p.E_ex, "E_in":p.E_in2, "I_e":p.ie, "C_m":p.cm_fsi, "g_L":p.gL2, "t_ref":p.tref})
		# Heterogenity in D1,D2 and FSI by using a distribution of resting potentials instead of a single value
		vid1 = np.random.uniform(p.vi1low,p.vi1hi,p.numAll)
		vidl = vid1.tolist()
		nest.SetStatus(std1,params="V_m",val=vidl)
		nest.SetStatus(std2,params="V_m",val=vidl)
		vifsi = np.random.uniform(p.vi2low,p.vi2hi,p.numFSI)
		vifsil = vifsi.tolist()
		nest.SetStatus(stfsi,params="V_m",val=vifsil)
		
		import random
		
		# Conn probabilities
		# d1 - d1 13%
		# d1 - d2 3%
		# d2 - d1 13.5%
		# d2 - d2 18%
					
		d1d1 = 0.13
		d1d2 = 0.03
		d2d1 = 0.135
		d2d2 = 0.18		
		
		d1fsi = 0.27
		d2fsi = 0.18
		# Connect the populations
		nest.RandomDivergentConnect(std1, std1, int(d1d1*p.numAll), weight=float(p.jd1d1), delay=p.delay21,options={'allow_autapses':False,'allow_multapses':False})
		nest.RandomDivergentConnect(std1, std2, int(d1d2*p.numAll), weight=float(p.jd1d2), delay=p.delay21,options={'allow_autapses':False,'allow_multapses':False})
		nest.RandomDivergentConnect(std2, std1, int(d2d1*p.numAll), weight=float(p.jd2d1), delay=p.delay21,options={'allow_autapses':False,'allow_multapses':False})
		nest.RandomDivergentConnect(std2, std2, int(d2d2*p.numAll), weight=float(p.jd2d2), delay=p.delay21,options={'allow_autapses':False,'allow_multapses':False})
		# FSI connections are faster 
		nest.RandomDivergentConnect(stfsi,std1, int(d1fsi*p.numAll),weight=float(p.jfsi),delay=p.delay22,options={'allow_autapses':False,'allow_multapses':False})
		nest.RandomDivergentConnect(stfsi,std2, int(d2fsi*p.numAll),weight=float(p.jfsi),delay=p.delay22,options={'allow_autapses':False,'allow_multapses':False})
		
		# Background input to D1,D2 and FSI
		noise_d1 = nest.Create('poisson_generator',1,{'rate':2500.})
		noise_d2 = nest.Create('poisson_generator',1,{'rate':2500.})
		noise_fsi = nest.Create('poisson_generator',1,{'rate':2500.})	

		#Connect background input			
		nest.DivergentConnect(noise_d1,std1,weight=w2,delay=1.0,model='static_synapse')
		nest.DivergentConnect(noise_d2,std2,weight=w2,delay=1.0,model='static_synapse')
		nest.DivergentConnect(noise_fsi,stfsi,weight=w2,delay=1.0,model='static_synapse')
		
		print 'Connect the cortex neuronal population to first subpopulation of MSN neurons'
	
		# Connect spike detectors to D1,D2 and FSIs to record spikes	
		detect_d1 = nest.Create("spike_detector")
		nest.SetStatus(detect_d1,{"withgid" : True , "withtime" : True })
		nest.ConvergentConnect(std1,detect_d1)

		detect_d2 = nest.Create("spike_detector")
		nest.SetStatus(detect_d2,{"withgid" : True , "withtime" : True })
		nest.ConvergentConnect(std2,detect_d2)
		
		detect_fsi = nest.Create("spike_detector")
		nest.SetStatus(detect_fsi,{"withgid" : True , "withtime" : True })
		nest.ConvergentConnect(stfsi,detect_fsi)
		
		# p.num01=150 is the size of every cortical input pool that converges to a single D1/D2 neuron.
		# So to generate a input population rate of say 10Hz, each neuron should spike with a rate of 10/150.
		pg_d1 = nest.Create('poisson_generator',1,{'rate':(rate/float(p.num01))})
		st01_list=[]

		pg_d2 = nest.Create('poisson_generator',1,{'rate':(rate/float(p.num01))})
		pg_fsi = nest.Create('poisson_generator',1,{'rate':rate})
		nest.DivergentConnect(pg_fsi,stfsi,weight=w3,delay=1.0)


		# parrot neurons represent the input pre-synaptic population. This could be done without parrots, but then no recording of
		# spikes possible since spike detector cannot be connected directly to poisson_generator (NEST limitation)
		for i in xrange((p.num1)):
			st01_list.append(nest.Create('parrot_neuron',p.num01))
		
		for i in xrange(p.numFSI):
			st01_list.append(nest.Create('parrot_neuron',int(p.num01)))
		
		# Connect poisson generators to the cortical pools
		nest.DivergentConnect(pg_d1,np.array(st01_list[:p.num1/2]).flatten().tolist(),weight=p.j02,delay=1.0,model='static_synapse')		
		nest.DivergentConnect(pg_d2,np.array(st01_list[p.num1/2:p.num1]).flatten().tolist(),weight=p.j02,delay=1.0,model='static_synapse')		
			
		# Record from 2 (example) cortical pools, to check if indeed MSNs are receiving the intended cortical input.
		detect_cortex_d1 = nest.Create("spike_detector")
		nest.SetStatus(detect_cortex_d1,{"withgid" : True , "withtime" : True })
		nest.ConvergentConnect(st01_list[0],detect_cortex_d1)

		detect_cortex_d2 = nest.Create("spike_detector")
		nest.SetStatus(detect_cortex_d2,{"withgid" : True , "withtime" : True })
		nest.ConvergentConnect(st01_list[p.num1/2+2],detect_cortex_d2)

		# Convergent connect from each pool to MSNs
		for i,msn in enumerate(std1):
			nest.ConvergentConnect(st01_list[i],[msn],weight=[w1],delay=1.0 ,model='static_synapse')

		for i,msn in enumerate(std2):
			nest.ConvergentConnect(st01_list[i+p.num1/2],[msn],weight=[w2],delay=1.0 ,model='static_synapse')		
		
		# Simulate
		nest.Simulate(p.runtime)

		# Read from spike detectors			
		dSD_d1 = nest.GetStatus(detect_d1)[0]
		evs_d1 = dSD_d1['events']['senders']
		ts_d1 = dSD_d1['events'] ['times']

		dSD_d2 = nest.GetStatus(detect_d2)[0]
		evs_d2 = dSD_d2['events']['senders']
		ts_d2 = dSD_d2['events'] ['times']

		dSD_fsi = nest.GetStatus(detect_fsi)[0]
		evs_fsi = dSD_fsi['events']['senders']
		ts_fsi = dSD_fsi['events'] ['times']

		dSD_cortex_d1 = nest.GetStatus(detect_cortex_d1)[0]
		evs_cortex_d1 = dSD_cortex_d1['events']['senders']
		ts_cortex_d1 = dSD_cortex_d1['events'] ['times']

		dSD_cortex_d2 = nest.GetStatus(detect_cortex_d2)[0]
		evs_cortex_d2 = dSD_cortex_d2['events']['senders']
		ts_cortex_d2 = dSD_cortex_d2['events'] ['times']
		
		# Store spike times and IDs in dictionaries
		spike_senders['d1'].append(evs_d1)
		spike_senders['d2'].append(evs_d2)
		spike_senders['ctx'].append(evs_cortex_d1)
		spike_senders['fsi'].append(evs_fsi)

		spike_times['d1'].append(ts_d1)
		spike_times['d2'].append(ts_d2)
		spike_times['ctx'].append(ts_cortex_d1)
		spike_times['fsi'].append(ts_fsi)

		# Calculate the mean firing rate of D1,D2and FSI populations
		secs = float(p.runtime)/1000.
		rateD1 = (len(ts_d1)/secs)/(float(p.num1/2))
		rateD2 = (len(ts_d2)/secs)/(float(p.num1/2))
		rateFSI = (len(ts_fsi)/secs)/(float(p.numFSI))
		rateCort_d1 = (len(ts_cortex_d1)/secs)/float(p.num01)
		rateCort_d2 = (len(ts_cortex_d2)/secs)/float(p.num01)

		rates['d1'].append(rateD1)
		rates['d2'].append(rateD2)
		rates['ctx'].append(rateCort_d1)
		rates['fsi'].append(rateFSI)

	# Dump in pickle files for future analysis
	import pickle
	print rates
	name = 'Senders_'+prefix+'.pickle'
	pickle.dump(spike_senders,open(name,"w"))
	
	name = 'Times_'+prefix+'.pickle' 
	pickle.dump(spike_times,open(name,"w"))

	name = 'Rates_'+prefix+'.pickle' 
	pickle.dump(rates,open(name,"w"))
		


