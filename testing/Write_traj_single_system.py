#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import MDAnalysis
from MDAnalysis.analysis.align import *
import sys

# ----------------------------------------
# VARIABLE DECLARATION

pdb = sys.argv[1]		# pdb file of protein system; needs to have the same # of atoms as the trajectories do
traj_loc = sys.argv[2]		# pointer to the directories of the systems
start = int(sys.argv[3])	# integer of first production run to analyze
end = int(sys.argv[4])		# integer of last production run to analyze
grab = int(sys.argv[5])		# step integer between grabbing frames
system = sys.argv[6]		# system descriptor
out = sys.argv[7]		# output file name w/ extension

selection = 'protein or nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname MG or resname PHX'
alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'

flush = sys.stdout.flush

# ----------------------------------------
# FUNCTIONS:

def ffprint(string):
        print '%s' %(string)
        flush()

def summary(nSteps):
	sum_file = open('%s.traj_writing.summary' %(system),'w')
	sum_file.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
	sum_file.write('To recreate this analysis, run this line in terminal:\n')
	for i in range(len(sys.argv)):
		sum_file.write('%s ' %(sys.argv[i]))
	sum_file.write('\n\n')
	sum_file.write('Output is written to:\n')
	sum_file.write('	ref_structure.pdb\n')
	sum_file.write('	%s\n' %(out))
	sum_file.write('\nTotal number of frames written to a new trajectory: %d\n' %(nSteps/grab))
	sum_file.write('\nAtom selection written to trajectory:\n')
	sum_file.write('	%s\n' %(selection))
	sum_file.write('\n Frames were aligned to:\n')
	sum_file.write('	%s\n' %(alignment))
	sum_file.close()

# ----------------------------------------
# MAIN:

u = MDAnalysis.Universe('%s' %(pdb))		
u_align = u.select_atoms(alignment)
u_sel = u.select_atoms(selection)

u_sel.translate(-u_align.center_of_mass())
u_sel.write('ref_structure.pdb')

nAtoms = len(u_sel.atoms)
pos0 = u_align.positions

nSteps = 0
with MDAnalysis.Writer('%s' %(out), nAtoms) as W:

	while start <= end:
		ffprint('Beginning to analyzing trajectory %02d.' %(start))
		u.load_new('%s/production.%s/production.%s.dcd' %(traj_loc,start,start))
		nFrames = len(u.trajectory)
		nSteps += nFrames
		for k in range(nFrames):
			if k%(grab) == 0:
				u.trajectory[k]
				u_sel.translate(-u_align.center_of_mass())
				R, rmsd = rotation_matrix(u_align.positions,pos0)
				u_sel.rotate(R)
				W.write(u_sel)

		start += 1 
W.close()

summary(nSteps)

