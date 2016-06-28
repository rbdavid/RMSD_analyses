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
#from fn_plotting import *

# ----------------------------------------
# VARIABLE DECLARATION

num_eigens = int(sys.argv[1])

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

legend_list = []
for i in range(nSys):
	legend_list.append(frame[i][0])

flush = sys.stdout.flush
changedir = os.chdir

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
        print '%s' %(string)
        flush()

# ----------------------------------------
# MAIN PROGRAM:

eigens = []
for i in range(num_eigens):
	eigens.append('%02d.projection.dat' %(i))

# PLOTTING EIGENVALUES...
eigvalues = np.loadtxt('eigenvalues.dat')
nSteps = len(eigvalues)
print nSteps
cumulative_eigval = np.zeros(nSteps)
total_eigval = 0
for j in range(nSteps):
	total_eigval += eigvalues[j]
	cumulative_eigval[j] = total_eigval

out1 = open('cumulative_eigenvalues.dat','w')
for j in range(nSteps):
	out1.write('%f   %f\n' %(eigvalues[j]/total_eigval,cumulative_eigval[j]))
out1.close()

plt.plot(cumulative_eigval[:10]/total_eigval,'k.')
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.xlabel('Eigen Number')
plt.ylim((0,1.01))
plt.ylabel('% of Total Eigenvalues')
plt.savefig('eigenvalues.png')
plt.close()

for j in range(num_eigens-1):
	proj0 = np.loadtxt('%02d.projection.dat' %(j))
	for k in range(j+1, num_eigens):
		proj1 = np.loadtxt('%02d.projection.dat' %(k))
		x = 0
		y = 0
		for l in range(nSys):
			y += frame[l][1]
			plt.plot(proj0[x:y],proj1[x:y],c=frame[l][2],marker=frame[l][3],ls='none',markersize=3)
			x += frame[l][1]
		plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
		leg = plt.legend(legend_list, bbox_to_anchor=(-0.05, 1.02, 1.1, 0.030), fontsize='10', loc=3, ncol=4, mode="expand", borderaxespad=0., markerscale=3,numpoints=1)
		plt.xlabel('Projection onto Eigenvector %s' %(j+1))
		plt.ylabel('Projection onto Eigenvector %s' %(k+1))
		plt.savefig('projection_%02d_%02d.png' %(j,k),dpi=300)
		plt.close()


proj0 = np.loadtxt('00.projection.dat')
proj1 = np.loadtxt('01.projection.dat')
proj2 = np.loadtxt('02.projection.dat')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = 0
y = 0
for j in range(nSys):
	y += frame[j][1]
	ax.scatter(proj0[x:y],proj1[x:y],proj2[x:y],c=frame[j][2],marker=frame[j][3],linewidths=0.2)
	x += frame[j][1]
leg = plt.legend(legend_list, bbox_to_anchor=(-0.05, 1.03, 1.1, .100), fontsize='10', loc=3, ncol=4, mode="expand", borderaxespad=0., markerscale=3,scatterpoints=1)
ax.set_xlabel('Projection onto Eigenvector 1')
ax.set_ylabel('Projection onto Eigenvector 2')
ax.set_zlabel('Projection onto Eigenvector 3')
plt.savefig('projection_00_01_02.png')
plt.close()

