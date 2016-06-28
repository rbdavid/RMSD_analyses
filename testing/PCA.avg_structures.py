#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:


# ----------------------------------------
# PREAMBLE:

import numpy as np
from numpy.linalg import *
import sys
import os
import sel_list

# ----------------------------------------
# VARIABLE DECLARATION

zeros = np.zeros
flush = sys.stdout.flush
eigen = np.linalg.eig
dot_prod = np.dot
changedir = os.chdir
sums = np.sum
sqrt = np.sqrt

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
        print '%s' %(string)
        flush()

# ----------------------------------------
# MAIN PROGRAM:

nSel = len(sel_list.sel)

datalist1 = np.loadtxt('AMBER_ssrna_atp.rmsd.dat')
#nSteps = len(datalist1)
#
#datalist2 = zeros((nSteps,nSel),dtype=np.float32)
#
#datalist2[:,0] = datalist1[:,4]
#datalist2[:,1] = datalist1[:,5]
#datalist2[:,2] = datalist1[:,6]
#datalist2[:,3] = datalist1[:,7]
#datalist2[:,4] = datalist1[:,8]
#datalist2[:,5] = datalist1[:,9]
#datalist2[:,6] = datalist1[:,10]
#datalist2[:,7] = datalist1[:,11]
#datalist2[:,8] = datalist1[:,12]
#datalist2[:,9] = datalist1[:,13]
#datalist2[:,10] = datalist1[:,14]
#datalist2[:,11] = datalist1[:,15]
#datalist2[:,12] = datalist1[:,16]
#datalist2[:,13] = datalist1[:,17]
#datalist2[:,14] = datalist1[:,18]
#datalist2[:,15] = datalist1[:,19]
#
#print datalist2[0]

centered_matrix = datalist1 - np.mean(datalist1,axis=0)          # average RMSD value of columns (axis=0)

xTx = dot_prod(centered_matrix.T, centered_matrix)
eigval, eigvec = eigen(xTx)
idx = eigval.argsort()[::-1]
eigval = eigval[idx]

out1 = open('eigenvalues.dat','w')
out2 = open('eigenvectors.dat','w')
for j in range(len(eigval)):
	out1.write('%f \n' %(eigval[j]))
	for k in range(len(eigvec[:,0])):
		out2.write('%f   ' %(eigvec[k,j]))
	out2.write('\n')
out1.close()
out2.close()

for j in range(len(eigval)):
	out = open('%02d.projection.dat' %(j), 'w')
	a_b = dot_prod(centered_matrix, eigvec[:,idx[j]])
	for k in range(len(a_b)):
		out.write('%f \n' %(a_b[k]))
	out.close()

