from bokeh.plotting import figure, output_notebook, show
from bokeh.charts import Bar, output_file, show
from bokeh.charts.attributes import cat
from numpy import cos, linspace
from collections import Counter,OrderedDict

import pandas as pd
import os
#x_file = open(direct+"/5_1.txt", "r")
# f=open(os.path.join(sub_dir,file))
#------------------------------------
#Store in secondary structures in a counter list
#------------------------------------

"""Dictionary to store the amino acids and their corresponding van der Walls volume"""
vdw_volume={}
vdw_volume['A']=67
vdw_volume['R']=148
vdw_volume['N']=96
vdw_volume['D']=91
vdw_volume['C']=86
vdw_volume['c']=86
vdw_volume['Q']=114
vdw_volume['E']=109
vdw_volume['G']=48
vdw_volume['H']=118
vdw_volume['I']=124
vdw_volume['L']=124
vdw_volume['K']=135
vdw_volume['M']=124
vdw_volume['F']=135
vdw_volume['P']=90
vdw_volume['S']=73
vdw_volume['T']=93
vdw_volume['W']=163
vdw_volume['Y']=141
vdw_volume['V']=105

vdw_volume_list=['67','148','96','91','86','114','109','48','118','124','135','90','73','93','163','141','105']


vdw_volume_list=sorted(vdw_volume_list, key=lambda x: float(x))
vdw_volume_list=[str(_) for _ in vdw_volume_list]
print vdw_volume_list
before_residue_list=[]
after_residue_list=[]
get=open('./../'+'angles_secondary_structure.txt','r')
for line in get:
	line=line.strip('\n')
	lines=line.split(',')
	cys1_ss=lines[-2:]
	cys2_ss=lines[-1:]+ lines[-2:-1]
	cys1_ss=",".join(cys1_ss)
	# cys1_ss=cys1_ss.strip('\n')


	cys2_ss=",".join(cys2_ss)
	# cys2_ss=cys2_ss.strip('\n')


	#------------------------------------#------------------------------------
	#Calculating the VDW volume of neighbouring amino acids
	#------------------------------------#---------------------------------
	"""Defining the neighbouring amino acids of the two Cys residues of a disulfide bond"""
	neighbouring_aminoacids1=[lines[37]] + [lines[47]] 
	neighbouring_aminoacids2=[lines[57]] + [lines[67].strip('\n')]
	#Van der walls volumes for residues before first Cys
	cys1_before=[(vdw_volume[neighbouring_aminoacids1[0]])]
	cys1_before=str(cys1_before)
	cys1_after=[(vdw_volume[neighbouring_aminoacids1[1]])]
	cys1_after=str(cys1_after)
	#Van der walls volumes for residues before second Cys
	cys2_before=[(vdw_volume[neighbouring_aminoacids2[0]])]
	cys2_before=str(cys2_before)
	cys2_after=[(vdw_volume[neighbouring_aminoacids2[1]])]
	cys2_after=str(cys2_after)
	

	before_residue_list.append(cys1_before)
	before_residue_list.append(cys2_before)

	after_residue_list.append(cys1_after)
	after_residue_list.append(cys2_after)

total=len(before_residue_list)
counts_before=Counter(before_residue_list)
counts_after=Counter(after_residue_list)
for key, value in counts_before.items():
	x=float(value) / float(total)
	counts_before[key] = x
for key, value in counts_after.items():
    counts_after[key] = float(value) / float(total)


y_list_before=[]
y_list_after=[]
for value in vdw_volume_list:
	y_list_before.append((counts_before['['+value+']']))
	y_list_after.append(counts_after['['+value+']'])

print len(vdw_volume_list)
print len(y_list_before)
print len(y_list_after)
x_range=range(len(vdw_volume_list))
print x_range
p = figure(plot_width=400, plot_height=400)
p.line(vdw_volume_list,y_list_before,label=cat(columns='Secondary Structure',sort=False))

#p.multi_line([x_range,x_range],[y_list_before,y_list_after],color=["firebrick", "navy"], alpha=[0.8, 0.3], line_width=4)

show(p)