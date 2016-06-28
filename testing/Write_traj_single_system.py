#!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
##!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
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

selection = 'protein'
alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'

flush = sys.stdout.flush

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
        print '%s' %(string)
        flush()

# ----------------------------------------
# MAIN PROGRAM:

u = MDAnalysis.Universe('%s' %(pdb))		
u_backbone = u.select_atoms('backbone')
u_align = u.select_atoms(alignment)
u_sel = u.select_atoms(selection)

u_sel.translate(-u_backbone.center_of_mass())
u_sel.write('ref_structure.pdb')

nAtoms = len(u_sel.atoms)
pos0 = u_align.positions

nSteps = 0
with MDAnalysis.Writer('%s' %(out), nAtoms) as W:

	while start <= end:
		u.load_new('%s/production.%s/production.%s.dcd' %(traj_loc,start,start))
		nSteps += len(u.trajectory)
		nFrames = len(u.trajectory)
		
		for k in range(nFrames):
			if k%(grab) == 0:
				u.trajectory[k]
				u_sel.translate(-u_backbone.center_of_mass())
				R, rmsd = rotation_matrix(u_align.positions,pos0)
				u_sel.rotate(R)
				W.write(u_sel)

		start += 1 
W.close()

