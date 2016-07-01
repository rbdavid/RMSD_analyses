#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:


# ----------------------------------------
# PREAMBLE:

import MDAnalysis
from MDAnalysis.analysis.align import *
import sys
import os
from sel_list import *
from res_list import *
from distance_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

system_loc = sys.argv[1]
pdb_file = sys.argv[2]
system = sys.argv[3]
start = int(sys.argv[4])
end = int(sys.argv[5])

flush = sys.stdout.flush

alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'

nRef = len(ref)
nSel = len(sel)

if sel[0][1] != alignment:
	print 'First atom selection needs to correspond with the alignment landmark... which it does not, you idiot'
	sys.exit()

# ----------------------------------------
# FUNCTIONS:

def ffprint(string):
        print '%s' %(string)
        flush()

def summary(nSteps):
	sum_file = open('%s.rmsd.summary' %(system),'w')
	sum_file.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
	sum_file.write('To recreate this analysis, run this line in terminal:\n')
	for i in range(len(sys.argv)):
		sum_file.write('%s ' %(sys.argv[i]))
	sum_file.write('\n\n')
	sum_file.write('output is written to:\n')
	sum_file.write('	%s.output\n' %(system))
	sum_file.write('	%s.rmsd.dat\n' %(system))
	sum_file.write('\nTotal number of steps analyzed: %d\n' %(nSteps))
	sum_file.write('\nAtom selections analyzed:\n')
	for i in range(nSel):
		sum_file.write('	%02d   %s   %s\n' %(i,sel[i][0],sel[i][1]))
	sum_file.write('\nSystems used as reference structures:\n')
	for i in range(nSys):
		sum_file.write('	%s, %s%s%s\n' %(ref[i][0],system_loc,ref[i][0],ref_name))
	sum_file.close()

# ----------------------------------------
# MAIN:

out1 = open('%s.output' %(system),'w',1) 

# LOOP THROUGH ALL REFERENCE STRUCTURES AND SAVE COORDINATES OF THE ATOM SELECTIONS
pos_list['']*nRef
out_list = []
for i in range(nRef):
	temp_file = '%s%s/%s' %(system_loc,ref[i][1],ref[i][2])
	out1.write(' %02d -- Reference structure: %s\n' %(temp_file))
	temp_ref = MDAnalysis.Universe(temp_file)
	temp_all = temp_ref.select_atoms('all')
	temp_align = temp_ref.select_atoms(alignment)
	temp_all.translate(-temp_align.center_of_mass())

	pos0 = []
	# SAVE REFERENCE COORDINATES FOR ALL SELECTIONS...
	for j in range(nSel):
		pos0.append(temp_ref.select_atoms(sel[j][1]).positions)

	pos_list[i] = pos0
	temp_file = open('%02d_ref.rmsd.dat' %(i))
	out_list.append(temp_file)

out1.write('Finished collecting the reference structure data\n')

# INITIALIZING UNIVERSE AND SAVING ATOM SELECTIONS...
u = MDAnalysis.Universe(pdb_file)
u_all = u.select_atoms('all')
u_align = u.select_atoms(alignment)
u_selection_list = []
for i in range(nSel):
	temp_sel = u.select_atoms(sel[i][1])
	u_selection_list.append([int(temp_sel.n_atoms),temp_sel])
	out1.write('%s corresponds to %s atom selection\n' %(sel[i][0],u_selection_list[i][1]))

# BEGINNING TO ANALYZE TRAJECTORIES
nSteps = 0
while start <= end:
	u.load_new('%s%s/Trajectories/production.%s/production.%s.dcd' %(system_loc,system,start,start))
	nSteps += len(u.trajectory)
	for ts in u.trajectory:
		u_all.translate(-u_align.center_of_mass())
		for i in range(nRef):
			R, rmsd = rotation_matrix(u_align.positions,pos_list[i][0])
			u_all.rotate(R)
			for j in range(nSel):
				rmsd = RMSD(u_selection_list[j][1],pos_list[i][j],u_selection_list[j][0])
				out_list[i].write('%f   ' %(rmsd))
			out_list[i].write('\n')
	out1.write('Finished analyzing trajectory %02d\n' %(start))
	start += 1
out1.write('Analyzed %d steps.' %(nSteps))

out1.close()
for i in range(nRef):
	out_list[i].close()

summary(nSteps)

