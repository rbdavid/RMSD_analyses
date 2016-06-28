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
from distance_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

ref_file = sys.argv[1]			# pointer to the pdb file to be used as the reference structure
traj_loc = sys.argv[2]			# pointer to the trajectory positions (or really the position where all systems are stored); look at this variable's use in the script
number = int(sys.argv[3])		# integer identifying which system is used as the reference structure; use python indexing for the ref_list variable

flush = sys.stdout.flush

alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'
important = 'protein'

ref_list = []
ref_list.append(['AMBER_apo', 21, 150])	
ref_list.append(['AMBER_atp', 21, 150])		
ref_list.append(['AMBER_ssrna', 21, 150])		
ref_list.append(['AMBER_ssrna_atp', 21, 150])	
ref_list.append(['AMBER_ssrna_adp_pi', 21, 150])	
ref_list.append(['AMBER_ssrna_adp', 21, 150])	
ref_list.append(['AMBER_ssrna_pi', 21, 150])		

nSys = len(ref_list)
nSel = len(sel)

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
	sum_file.write('	%s.output\n' %(ref_list[number][0]))
	sum_file.write('	%s.rmsd.dat\n' %(ref_list[number][0]))
	sum_file.write('\nTotal number of steps analyzed: %d\n' %(nSteps))
	sum_file.write('\nAtom selections analyzed:\n')
	for i in range(nSel):
		sum_file.write('	%02d   %s   %s\n' %(i,sel[i][0],sel[i][1]))
	sum_file.write('\nSystems analyzed:\n')
	for i in range(nSys):
		sum_file.write('	%s, Trajectories %03d to %03d\n' %(ref_list[i][0],ref_list[i][1],ref_list[i][2]))
	sum_file.close()

# ----------------------------------------
# MAIN:

out1 = open('%s.output' %(ref_list[number][0]),'w',1) 
out2 = open('%s.rmsd.dat' %(ref_list[number][0]),'w')

out1.write('Reference structure: %s\n' %(ref_file))
ref = MDAnalysis.Universe(ref_file)
ref_all = ref.select_atoms('all')
ref_align = ref.select_atoms(alignment)
ref_all.translate(-ref_align.center_of_mass())
pos0 = ref_align.positions

# SAVE COORDINATES FOR ALL SELECTIONS...
pos_list = []
for i in range(nSel):
	selection = sel[i][1]
	temp_sel = ref.select_atoms(selection)
	temp_pos = temp_sel.positions
	pos_list.append(temp_pos)

out1.write('Finished collecting the reference structure data\n')
nSteps = 0
# INITIALIZING UNIVERSES, LOADING TRAJECTORIES IN, ANALYZING, ETC...
for i in range(nSys):
	out1.write('Loading in Trajectories from %s\n' %(ref_list[i][0]))
	u = MDAnalysis.Universe('%s%s/truncated.pdb' %(traj_loc,ref_list[i][0]))
	
	u_align = u.select_atoms(alignment)
	u_important = u.select_atoms(important)

	u_selection_list = []
	for a in range(nSel):
		selection = sel[a][1]
		temp_sel = u.select_atoms(selection)
		u_selection_list.append([int(temp_sel.n_atoms),temp_sel])
		out1.write('%s corresponds to %s atom selection\n' %(sel[a][0],u_selection_list[a][1]))

	count = 0
	out1.write('Beginning trajectory analysis from system %s\n' %(ref_list[i][0]))
	a = ref_list[i][1]
	while a <= ref_list[i][2]:
		u.load_new('%s%s/Trajectories/production.%s/production.%s.dcd' %(traj_loc,ref_list[i][0],a,a))
		nSteps += len(u.trajectory)
		for ts in u.trajectory:
			u_important.translate(-u_align.center_of_mass())
			R, rmsd = rotation_matrix(u_align.positions,pos0)
			u_important.rotate(R)
			for j in range(nSel):
				temp_pos = u_selection_list[j][1].positions
				rmsd = RMSD(temp_pos,pos_list[j],u_selection_list[j][0])
				out2.write('%f   ' %(rmsd))
			out2.write('\n')

		out1.write('Finished analyzing Trajectory: %s%s/Trajectories/production.%s/production.%s.dcd\n' %(traj_loc,ref_list[i][0],a,a))
		a +=1
	out1.write('Analyzed %d frames from system %s\n' %(count,ref_list[i][0]))

out1.close()
out2.close()
summary(nSteps)

