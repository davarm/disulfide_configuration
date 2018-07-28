import os
from bokeh.models import Jitter
from bokeh.plotting import figure, show, output_file
import math
from bokeh.charts.utils import cycle_colors


colors=['Black','Red','Lime','Blue','Yellow','Cyan','Magenta','Silver','Gray','Maroon','Olive','Green','Purple','Teal','Navy']

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
	
	cys_distance=math.log(int(lines[2])-int(lines[4]))
	#print cys1_ss
	#print cys2_ss
	#---
	#Have to store in lists based on secondary structure
	#Will assign the line based to a dictionary:
	
	#Secondary structure list will determine the order of dictionaries

	if cys1_ss in secondary_structure_list:
		#print 'YES',cys1_ss
		try:
			secondary_structure_dict[cys1_ss].append(cys_distance)
		except KeyError:
			secondary_structure_dict[cys1_ss]=[cys_distance]
	
	if cys1_ss not in secondary_structure_list:
		#print 'NO',cys1_ss
		try:
			secondary_structure_dict[cys2_ss].append(cys_distance)
		except KeyError:
			secondary_structure_dict[cys2_ss]=[cys_distance]



p = figure(plot_width=1000, plot_height=400, x_range=secondary_structure_list, y_range=(0,8),
           title="X1 Probabilities",)

#palette = all_palettes['Category20c'][15]
for i,value in enumerate(secondary_structure_list):
	y1=secondary_structure_dict[value]

	p.circle(x={'value': i+1, 'transform': Jitter(width=0.3)}, y=y1,
         color=colors[i], alpha=0.1,size=5)



p.yaxis.axis_label = "DISH Probability Score"
#p.xaxis.axis_label = "Amino Acid"
p.xaxis.major_label_text_font_size='12pt'
p.yaxis.major_label_text_font_size='12pt'
p.yaxis.axis_label_text_font_size='13pt'
p.xaxis.axis_label_text_font_size='13pt'
p.xaxis.major_label_orientation='vertical'


#print palette
show(p)


















