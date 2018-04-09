from bokeh.plotting import figure, output_notebook, show
from bokeh.charts import Bar, output_file, show
from bokeh.charts.attributes import cat
from numpy import cos, linspace
from collections import Counter,OrderedDict

import pandas as pd
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
'Helix,Helix']


get=open('angles_secondary_structure.txt','r')
for line in get:
	line=line.strip('\n')
	lines=line.split(',')
	cys1_ss=lines[-2:]
	cys2_ss=lines[-1:]+ lines[-2:-1]
	cys1_ss=",".join(cys1_ss)
	# cys1_ss=cys1_ss.strip('\n')


	cys2_ss=",".join(cys2_ss)
	# cys2_ss=cys2_ss.strip('\n')
	
	if cys1_ss in secondary_structure_list:
		secondary_structure.append(cys1_ss)
	else:
		secondary_structure.append(cys2_ss)

total=len(secondary_structure)
print total
counts=Counter(secondary_structure)
for key, value in counts.items():
    counts[key] = float(value) / float(total)
print counts
#counts=sorted(counts, key=counts.get, reverse=True)
#counts=sorted(counts.items(), key=lambda item: item[1])
print counts
df=pd.DataFrame(counts.items(), columns=['Secondary Structure', 'Frequency'])
print df
df = df[['Frequency','Secondary Structure']]
print df
df=df.sort(['Frequency','Secondary Structure'], ascending=False)
#df=df.sort(inplace=True)
print df
#------------------------------------
#Making x,y list
#------------------------------------
x=[]
y=[]
print ''
i=1
for key in counts:
	
	x.append(i)
	frequency=float(counts[key])/float(total)
	y.append(frequency)
	i=i+1

print x
y=sorted(y, key=float, reverse=True)


#x = linspace(-6, 6, 100)
#y = cos(x)
#p = figure(width=500, height=500,background_fill_color='lightgrey')
#p.line(x, y, color="blue", alpha=0.5,line_width=2,)
#p.square(x, y, color="blue", alpha=1,)

#------------------------------------
# Bar graph
#------------------------------------
p = Bar(df,values='Frequency',label=cat(columns='Secondary Structure',sort=False),legend=None,color='blue',ylabel='Frequency (%)')
p.xaxis.major_label_text_font_size='10pt'
p.yaxis.major_label_text_font_size='10pt'
p.yaxis.axis_label_text_font_size='11pt'
p.xaxis.axis_label_text_font_size='11pt'
show(p)