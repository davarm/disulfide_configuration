import os
from amino_acids_sorted_native_frequency import plotting_amino_acid_frequncy
#------------------------------------
#Store in secondary structures in a counter list
#------------------------------------
secondary_structure=[]
secondary_structure_list=['Turn,Coil',
'Coil,Bridge',
'Coil,Helix',
'Strand,Coil',
'Helix,Bridge',
'Bridge,Turn',
'Strand,Bridge',	
'Turn,Turn',
'Turn,Helix',
'Strand,Helix',
'Bridge,Bridge',
'Strand,Strand',
'Strand,Turn',
'Helix,Helix',
'Coil,Coil']


#Make one large dictionary,secondary structure combination as key:
secondary_structure_dict={}
get=open('./../'+'angles_secondary_structure.txt','r')
for line in get:
	line=line.strip('\n')
	lines=line.split(',')
	cys1_ss=lines[-2:]
	cys2_ss=lines[-1:]+ lines[-2:-1]
	cys1_ss=",".join(cys1_ss)
	cys2_ss=",".join(cys2_ss)
	
	#print cys1_ss
	#print cys2_ss
	#---
	#Have to store in lists based on secondary structure
	#Will assign the line based to a dictionary:
	
	#Secondary structure list will determine the order of dictionaries

	if cys1_ss in secondary_structure_list:
		#print 'YES',cys1_ss
		try:
			secondary_structure_dict[cys1_ss].append(line)
		except KeyError:
			secondary_structure_dict[cys1_ss]=[line]
	
	if cys1_ss not in secondary_structure_list:
		#print 'NO',cys1_ss
		try:
			secondary_structure_dict[cys2_ss].append(line)
		except KeyError:
			secondary_structure_dict[cys2_ss]=[line]


for key in secondary_structure_dict:
	print key
	plotting_amino_acid_frequncy(secondary_structure_dict[key],key) 