import numpy as np

from bokeh.models import Jitter
from bokeh.plotting import figure, show, output_file


#------------------------------------#------------------------------------
#Analysing the distribution of DISH probabilities for True and False connections
#------------------------------------#------------------------------------

#Database previously generated that provides the DISH probability and angle for all possible connections of a peptide from experimental database, used as inputs for SVM detection
true_x1_prob=[]
false_x1_prob=[]
true_x2_prob=[]
false_x2_prob=[]
true_product_prob=[]
false_product_prob=[]



get=open('angles_probs.txt','r')
for line in get:
	line=line.strip('\n')
	lines=line.split(',')

	#------------------------------------
	#Sort into true and false groups
	#------------------------------------

	if lines[4]=='true':
		true_x1_prob.append(float(lines[1]))
		true_x2_prob.append(float(lines[3]))
		true_product_prob.append(float(lines[1])*float(lines[3]))
	
	if lines[4]=='false':
		false_x1_prob.append(float(lines[1]))
		false_x2_prob.append(float(lines[3]))
		false_product_prob.append(float(lines[1])*float(lines[3]))

test=['True','False']
p = figure(plot_width=500, plot_height=400, x_range=test, y_range=(0.49,1.01),
           title="X1 Probabilities",)

y1 = true_x1_prob
y2 = false_x1_prob

p.circle(x={'value': 1, 'transform': Jitter(width=0.3)}, y=y1,
         color="navy", alpha=0.6,size=5)

p.circle(x={'value': 2, 'transform': Jitter(width=0.3)}, y=y2,
         color="firebrick", alpha=0.6,size=5)


p.yaxis.axis_label = "DISH Probability Score"
#p.xaxis.axis_label = "Amino Acid"
p.xaxis.major_label_text_font_size='12pt'
p.yaxis.major_label_text_font_size='12pt'
p.yaxis.axis_label_text_font_size='13pt'
p.xaxis.axis_label_text_font_size='13pt'


output_file("x1_probabilities.html")

show(p)
