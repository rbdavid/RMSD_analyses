#!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
##!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

sel = []
sel.append(['aligned_CAs','protein and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307) and name CA'])
sel.append(['aligned_betas','protein and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307) and not name H*'])
sel.append(['full_protein','protein and not name H*'])
sel.append(['full_backbone','backbone'])
sel.append(['protein-19','protein and not resid 0:19 and not name H*'])
sel.append(['backbone-19','backbone and not resid 0:19'])
sel.append(['motif_1','protein and resid 25:35 and not name H*'])
sel.append(['motif_1a','protein and resid 55:60 and not name H*'])
sel.append(['Lb3b4_1','protein and resid 76:89 and not name H*'])
sel.append(['motif_1b','protein and resid 94:100 and not name H*'])
sel.append(['motif_2','protein and resid 117:124 and not name H*'])
sel.append(['motif_3','protein and resid 148:150 and not name H*'])
sel.append(['post_m_3','protein and resid 151:158 and not name H*'])
sel.append(['motif_4','protein and resid 195:202 and not name H*'])
sel.append(['motif_4a','protein and resid 219:235 and not name H*'])
sel.append(['motif_5','protein and resid 241:250 and not name H*'])
sel.append(['motif_6','protein and resid 286:297 and not name H*'])
sel.append(['beta_wedges','protein and resid 260:280 and not name H*'])
sel.append(['d3_gate','protein and resid 434:442 and not name H*'])
sel.append(['a2_1','protein and resid 61:69 and not name H*'])

sel.append(['nucleic','(nucleic or resname A5 or resname A3 or resname U5) and not name H*'])
#sel.append(['ATP','resname atp and not name H*'])
#sel.append(['ADP','resname adp and not name H*'])
#sel.append(['PHX','resname PHX and not name H*'])
#sel.append(['MG','resname MG'])

#sel.append(['',''])

