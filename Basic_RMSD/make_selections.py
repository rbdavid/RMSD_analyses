
import MDAnalysis
import sys

# ----------------------------------------
# HOMEMADE ANALYSIS SELECTIONS
selections = []
selections.append(['aligned_CAs','protein and name CA and (resid 14:19 44:49 67:69 84:88 106:110 136:141 159:163 184:188 208:212 230:234 247:252 295:299)'])
selections.append(['aligned_betas','protein and (resid 14:19 44:49 67:69 84:88 106:110 136:141 159:163 184:188 208:212 230:234 247:252 295:299) and not name H*'])
selections.append(['full_protein','protein and not name H*'])
selections.append(['full_backbone','backbone'])
selections.append(['protein-10','protein and not resid 0:13 and not name H*'])
selections.append(['backbone-10','backbone and not resid 0:13'])
selections.append(['motif_1','protein and resid 19:30 and not name H*'])
selections.append(['motif_1a','protein and resid 49:53 and not name H*'])
selections.append(['a2_1','protein and resid 54:63 and not name H*'])
selections.append(['motif_1c','protein and resid 71:74 and not name H*'])
selections.append(['motif_1b','protein and resid 87:93 and not name H*'])
selections.append(['motif_2','protein and resid 111:118 and not name H*'])
selections.append(['motif_3','protein and resid 142:144 and not name H*'])
selections.append(['post_m_3','protein and resid 145:156 and not name H*'])
selections.append(['motif_4','protein and resid 190:196 and not name H*'])
selections.append(['motif_4a','protein and resid 213:218 221 and not name H*'])
selections.append(['motif_5','protein and resid 234:242 and not name H*'])
selections.append(['motif_6','protein and resid 281:290 and not name H*'])

# ----------------------------------------
# STANDARD RESIDUES:
nucleic = ['A5','A3','A','G5','G3','G','C5','C3','C','T5','T3','T','U5','U3','U']
triphosphate = ['atp','adp','PHX']
other = ['MG']

# ----------------------------------------
# HOMEMADE ATOM SELECTION STRINGS FOR THE STANDARD RESIDUES:
sugar = "name C5' C4' O4' C1' C3' C2' O2' " + " C5* C4* O4* C1* C3* O3* C2* O2* "		# NO HYDROGENS; DOES NOT INCLUDE THE O5' atom (which I will include in the phosphate atom selection string...; the atoms with * are found in triphosphates;
sugar_5= sugar + " O5'"		# NO HYDROGENS
sugar_3= sugar + " O3' "	# NO HYDROGENS
base = 'name N9 C8 N7 C5 C6 N6 N1 C2 N3 C4 O6 N4 C2 O2 O4'	# NO HYDROGENS; selection string that will select all appropriate atoms for any of the nucleic residues...
a_phos = 'name O5* O2A O1A PA O3A'
b_phos = 'name PB O1B O2B O3B'
g_phos = 'name PG O1G O2G O3G'
inorg_phos = 'name P O1 O2 O3 O4'	# NO HYDROGENS

# ----------------------------------------
# FUNCTION USED TO MAKE ANY OF THE HOMEMADE ATOM SELECTIONS FOR THE STANDARD RESIDUES

def make_selections(analysis_universe,ref_universe,resname,resid,output_file,selection_list,nAtoms,ref_pos,count):
	"""A function that takes in a residue name and creates a non-standard MDAnalysis atom selection
	
	Usage: make_selection(........)
	
	Arguments:
		analysis_universe: MDAnalysis Universe object to be used as the analysis universe.
		reference_universe: MDAnalysis Universe object to be used as the reference universe.
		resname: string of the residue name;
		resid: int of the residue ID number;
		output_file: file object that is to be written to;
        """

	# ----------------------------------------
	# CREATING THE NUCLEIC SELECTIONS
	if resname in nucleic:
		# CREATING THE SLECTION FOR THE BASE OF NUCLEIC RESIDUES
		sel_string = 'resname %s and resid %d and %s' %(resname,resid,base)
		u_temp = analysis_universe.select_atoms(sel_string)
		selection_list.append(u_temp)
		nAtoms.append(u_temp.n_atoms)
		ref_temp = ref_universe.select_atoms(sel_string)
		ref_pos.append(ref_temp.positions)
		if u_temp.n_atoms != ref_temp.n_atoms:
			print 'Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string)
			sys.exit()
		output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
		count +=1

		# CREATING THE SLECTION FOR THE SUGAR OF NUCLEIC RESIDUES
		if resname in ['A5','G5','C5','T5','C5']:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar_5)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			if u_temp.n_atoms != ref_temp.n_atoms:
				print 'Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string)
				sys.exit()
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1
			return

		elif resname in ['A3','U3','C3','G3']:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar_3)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			if u_temp.n_atoms != ref_temp.n_atoms:
				print 'Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string)
                                sys.exit()
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1

		else:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			if u_temp.n_atoms != ref_temp.n_atoms:
				print 'Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string)
				sys.exit()
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1

		# CREATING THE SLECTION FOR THE PHOSPHATE OF NUCLEIC RESIDUES
		sel_string = "(resname %s and resid %s and name P OP1 OP2 O5') or (resid %s and name O3')" %(resname,resid,analysis_universe.residues[resid-1].resid) 
		u_temp = analysis_universe.select_atoms(sel_string)
		selection_list.append(u_temp)
		nAtoms.append(u_temp.n_atoms)
		ref_temp = ref_universe.select_atoms(sel_string)
		ref_pos.append(ref_temp.positions)
		output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
		count += 1
		return

	# ----------------------------------------
	# CREATING THE TRIPHOSPHATE ATOM SELECTIONS
	elif resname in triphosphate:
		if resname in ['atp','adp']:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,base)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			if u_temp.n_atoms != ref_temp.n_atoms:
				print 'Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string)
				sys.exit()
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			if u_temp.n_atoms != ref_temp.n_atoms:
				print 'Number of atoms do not match for selection %02d, %s, %s' %(count,resname,sel_string)
				sys.exit()
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1

		if resname == 'atp':
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,a_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,b_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,g_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1
			return

		elif resname == 'adp':
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,a_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,b_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1
			return

		elif resname == 'PHX':
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,inorg_phos)
			u_temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(u_temp)
			nAtoms.append(u_temp.n_atoms)
			ref_temp = ref_universe.select_atoms(sel_string)
			ref_pos.append(ref_temp.positions)
			output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
			count +=1
			return

	# ----------------------------------------
	# CREATING ANY REMAINING ATOM SELECTIONS...
	elif resname in other:
		sel_string = 'resname %s and resid %d' %(resname,resid)
		u_temp = analysis_universe.select_atoms(sel_string)
		selection_list.append(u_temp)
		nAtoms.append(u_temp.n_atoms)
		ref_temp = ref_universe.select_atoms(sel_string)
		ref_pos.append(ref_temp.positions)
		output_file.write('%02d   %s   %s\n' %(count,resname,sel_string))
		count +=1
		return

# ----------------------------------------


