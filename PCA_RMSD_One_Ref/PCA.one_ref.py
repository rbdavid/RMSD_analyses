#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:


# ----------------------------------------
# PREAMBLE:

import numpy as np
from numpy.linalg import *
import sys
import os

# ----------------------------------------
# VARIABLE DECLARATION

zeros = np.zeros
flush = sys.stdout.flush
eigen = np.linalg.eig
dot_prod = np.dot
changedir = os.chdir
sums = np.sum
sqrt = np.sqrt

system = sys.argv[1]

# ----------------------------------------
# FUNCTIONS:

def ffprint(string):
        print '%s' %(string)
        flush()

# ----------------------------------------
# MAIN:

datalist1 = np.loadtxt('%s.rmsd.dat' %(system))

centered_matrix = datalist1 - np.mean(datalist1,axis=0)          # average RMSD value of columns (axis=0)

xTx = dot_prod(centered_matrix.T, centered_matrix)
eigval, eigvec = eigen(xTx)
idx = eigval.argsort()[::-1]
eigval = eigval[idx]

nVec = len(eigval)
cumulative_eigval = np.zeros(nVec)
total_eigval = 0
for i in range(nVec):
	total_eigval += eigval[i]
	cumulative_eigval[i] = total_eigval

out1 = open('eigenvalues.dat','w')
out2 = open('eigenvectors.dat','w')
for i in range(nVec):
	out1.write('%f   %f   %f   %f\n' %(eigval[i],eigval[i]/total_eigval,cumulative_eigval[i],cumulative_eigval[i]/total_eigval))
	for j in range(len(eigvec[:,0])):
		out2.write('%f   ' %(eigvec[j,i]))
	out2.write('\n')
out1.close()
out2.close()

for i in range(nVec):
	out = open('%02d.projection.dat' %(i), 'w')
	a_b = dot_prod(centered_matrix, eigvec[:,idx[i]])
	for j in range(len(a_b)):
		out.write('%f \n' %(a_b[j]))
	out.close()

