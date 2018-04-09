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


amino_acids=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

#------------------------------------#------------------------------------
#NATIVE AMINO ACID FREQUENCY
#------------------------------------#------------------------------------
native_amino_freq=[0.0892,0.0566,0.0395,0.0542,0.0126,0.0384,0.0615,0.0716,0.0222,0.057,0.0983,0.0506,0.024,0.0393,0.0487,0.0681,0.0559,0.013,0.0295,0.0681]


print len(native_amino_freq)
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
	#Appending amino acid length
	#------------------------------------#---------------------------------
	"""Defining the neighbouring amino acids of the two Cys residues of a disulfide bond"""
	neighbouring_aminoacids1=[lines[37]] + [lines[47]] 
	neighbouring_aminoacids2=[lines[57]] + [lines[67].strip('\n')]
	#Van der walls volumes for residues before first Cys
	cys1_before=(neighbouring_aminoacids1[0])
	cys1_after=neighbouring_aminoacids1[1]
	#Van der walls volumes for residues before second Cys
	cys2_before=neighbouring_aminoacids2[0]
	cys2_before=str(cys2_before)
	cys2_after=neighbouring_aminoacids2[1]
	

	before_residue_list.append(cys1_before)
	before_residue_list.append(cys2_before)

	after_residue_list.append(cys1_after)
	after_residue_list.append(cys2_after)

#print before_residue_list
#print after_residue_list
before_residue_list=['C' if x=='c' else x for x in before_residue_list]
after_residue_list=['C' if x=='c' else x for x in after_residue_list]


#------------------------------------#------------------------------------
#The total is the length of the residue list, or the total number of neighbouring residues 
#Used to then divide the number of each occurence to find the frequency
#------------------------------------#------------------------------------

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
for value in amino_acids:
	y_list_before.append((counts_before[value]))
	y_list_after.append(counts_after[value])

print len(y_list_before)
x_range=range(1,21)
print len(x_range)
print x_range
p = figure(plot_width=600, plot_height=400,x_range=amino_acids)
#p.multi_line([x_range,x_range,x_range],[y_list_before,y_list_after,native_amino_freq],color=["darkmagenta", "mediumblue","black"], alpha=[0.9, 0.9,0.5], line_width=3)
#p.multi_line([x_range,x_range],[y_list_before,y_list_after],color=["darkmagenta", "mediumblue"], alpha=[0.9, 0.9], line_width=3, legend=["before",'after'])

p.line(x_range,y_list_before,color="darkmagenta",alpha=0.9,line_width=3,legend='before')
p.line(x_range,y_list_after,color="mediumblue",alpha=0.9,line_width=3,legend='after')
p.line(x_range,native_amino_freq,color="black",alpha=0.9,line_dash=[6,3],line_width=3,legend='native')
#p.legend.orientation = "left_top"
p.legend.location = "top_right"
p.legend.orientation="horizontal"
p.yaxis.axis_label = "Frequency"
p.xaxis.axis_label = "Amino Acid"
p.xaxis.major_label_text_font_size='12pt'
p.yaxis.major_label_text_font_size='12pt'
p.yaxis.axis_label_text_font_size='13pt'
p.xaxis.axis_label_text_font_size='13pt'
show(p)