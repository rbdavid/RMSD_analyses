#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:



# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d import Axes3D
import sel_list
from plotting_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

file1 = sys.argv[1]

frame = []
frame.append(['Apo',650000,'steelblue','.'])
frame.append(['ATP',650000,'cadetblue','.'])
frame.append(['ssRNA',650000,'turquoise','.'])
frame.append(['ssRNA+ATP',650000,'forestgreen','.'])
frame.append(['ssRNA+ADP+Pi',650000,'limegreen','.'])
frame.append(['ssRNA+ADP',650000,'orangered','.'])
frame.append(['ssRNA+Pi',650000,'crimson','.'])
#frame.append([,,,])

nSys = len(frame)
nSel = len(sel_list.sel)

legend_list = []
for i in range(nSys):
	legend_list.append(frame[i][0])

flush = sys.stdout.flush

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
        print '%s' %(string)
        flush()

# ----------------------------------------
# MAIN PROGRAM:

data1 = np.loadtxt(file1)

#for i in range(nSel-1):
#	for j in range(i+1,nSel):
#		x = 0
#		y = 0
#		for k in range(nSys):
#			y += frame[k][1]
#			plt.plot(data1[x:y,i],data1[x:y,j],c=frame[k][2],marker=frame[k][3],ls='None',ms=0.2)
#			x += frame[k][1]
#		plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#		plt.xlabel('RMSD data for %s' %(sel_list.sel[i][0]))
#		plt.ylabel('RMSD data for %s' %(sel_list.sel[j][0]))
#		leg = plt.legend(legend_list,bbox_to_anchor=(-0.05, 1.03, 1.1, .100),fontsize='10',loc=3,ncol=4,mode="expand",borderaxespad=0.,markerscale=100,numpoints=1)
#		plt.savefig('%02d.%02d.png' %(i,j))
#		plt.close()

for i in range(nSel):
	events, edges, patches = plt.hist([
		data1[0:650000,i],
		data1[650000:1300000,i],
		data1[1300000:1950000,i],
		data1[1950000:2600000,i],
		data1[2600000:3250000,i],
		data1[3250000:3900000,i],
		data1[3900000:4550000,i]],
		bins=100, histtype='bar',
		color=[frame[0][2],frame[1][2],frame[2][2],frame[3][2],frame[4][2],frame[5][2],frame[6][2]],stacked=True)
	
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel('RMSD data for %s' %(sel_list.sel[i][0]))
	plt.ylabel('Frequency')
	plt.xlim((min(data1[:,i]),max(data1[:,i])))
	leg = plt.legend(legend_list,bbox_to_anchor=(-0.05, 1.03, 1.1, .100),fontsize='10',loc=3,ncol=4,mode="expand",borderaxespad=0.,markerscale=100,numpoints=1)
	plt.savefig('%02d.hist1d.png' %(i),dpi=200)
	plt.close()

