#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from plotting_functions import *
from sel_list import *

stdev = np.std
sqrt = np.sqrt
nullfmt = NullFormatter()

file1 = sys.argv[1]
time_offset = int(sys.argv[2])		# units of ns 

sys = []
sys.append(['Apo','steelblue',','])
sys.append(['ATP','cadetblue',','])
sys.append(['ssRNA','turquoise',','])
sys.append(['ssRNA+ATP','forestgreen',','])
sys.append(['ssRNA+ADP+Pi','limegreen',','])
sys.append(['ssRNA+ADP','orangered',','])
sys.append(['ssRNA+Pi','crimson',','])

nSys = len(sys)
nSel = len(sel)

# ----------------------------------------
# MAIN PROGRAM:

datalist1 = np.loadtxt(file1)
nSteps = len(datalist1)/nSys
print 'Number of steps analyzed for each system is = %d' %(nSteps)

time = np.zeros(nSteps)
for i in range(nSteps):
	time[i] = time_offset + (i*0.002)		# units of time in ns; each frame is separated by 0.002 ns 

for i in range(nSel):
	selection = sel[i][0]
	
	fig = plt.figure(figsize=(12,10))

	y_min = int(np.min(datalist1[:,i]))
	y_max = int(np.max(datalist1[:,i])+1)

	x = 0
	y = 0
	for j in range(nSys):
		system = sys[j][0]
		y += nSteps
		
		temp = fig.add_subplot(nSys,1,j+1)
		plt.locator_params(axis='y',nbins=6)
		plt.locator_params(axis='x',nbins=13)
		
		temp.plot(time[:],datalist1[x:y,i],c=sys[j][1],marker=sys[j][2],linestyle='None',label=sys[j][0])

		x += nSteps
		temp.set_xlim((time_offset,1500))
		temp.set_ylim((y_min,y_max))
		temp.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
		temp.legend(fontsize=14,loc=2,frameon=False,markerfirst=False,numpoints=None,scatterpoints=None,bbox_to_anchor=(-0.015,1.15),handlelength=0)

		if j != nSys-1:
			temp.set_xticklabels([])

	fig.text(0.5,0.04,'Time (ns)',ha='center',va='center',fontsize=20)
	fig.text(0.06,0.5,'RMSD ($\AA$)',ha='center',va='center',rotation='vertical',fontsize=20)
	
	plt.savefig('%02d.%s.comparison.png' %(i,selection),dpi=200)
	plt.close()
	fig.clf()

