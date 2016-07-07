#!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
##!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:


# ----------------------------------------
# PREAMBLE:

import numpy as np
from numpy.linalg import *
from sel_list import *
import sys
import os

# ----------------------------------------
# VARIABLE DECLARATION

rmsd_file = sys.argv[1]
important_eigens = int(sys.argv[1])

square = np.square
abs = np.absolute
sums = np.sum
sqrt = np.sqrt
zeros = np.zeros
flush = sys.stdout.flush
eigen = np.linalg.eig
dot_prod = np.dot
changedir = os.chdir

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
        print '%s' %(string)
        flush()

# ----------------------------------------
# MAIN PROGRAM:

nSel = len(sel)

datalist1 = np.loadtxt(rmsd_file)

nSteps = int(len(datalist))

for i in range(nSel):
	ffprint('Beginning PCA analysis of %s selection' %(sel[i][0]))
	changedir('%s' %(sel[i][0]))
	datalist1 = np.loadtxt('%s.rmsd_matrix.dat' %(sel_list.sel[i][0]))

	nSteps = int(datalist1[-1][1]+1)
	sq_matrix = zeros((nSteps,nSteps),dtype=float)

	for j in range(len(datalist1)):
		sq_matrix[int(datalist1[j][0]),int(datalist1[j][1])] = float(datalist1[j][2])
		sq_matrix[int(datalist1[j][1]),int(datalist1[j][0])] = float(datalist1[j][2])

	centered_matrix = sq_matrix - np.mean(sq_matrix,axis=0)

	out = open('euclid_distances.dat','w')
	for j in range(nSteps):
		eu_dist = sqrt(sums(square(centered_matrix[j,:])))
		out.write('%d    %f\n' %(j,eu_dist))
	out.close()

	xTx = dot_prod(centered_matrix.T, centered_matrix)
	eigval, eigvec = eigen(xTx)
	idx = eigval.argsort()[::-1]
	eigval = eigval[idx]
	# Eigvec is still unorganized... Columns correspond to the eigvectors...
	
	out = open('eigenvalues.dat','w')
	for j in range(len(eigval)):
		out.write('%f \n' %(eigval[j]))
	out.close()

	for j in range(important_eigens):
		out = open('%02d.projection.dat' %(j), 'w')
		a_b = dot_prod(centered_matrix,eigvec[:,idx[j]])
		for k in range(len(a_b)):
			out.write('%f \n' %(a_b[k]))
		out.close()

	changedir('..')

	ffprint('Finishing PCA analysis of %s selection' %(sel_list.sel[i][0]))

