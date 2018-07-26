import pandas as pd
import numpy as np

df = pd.read_csv('xx.csv')
print len(df)

ss_dict      = {}
ss_dict['H'] = 'a-Helix'
ss_dict['B'] = 'B-Bridge'
ss_dict['E'] = 'Strand'
ss_dict['G'] = '3-Helix'
ss_dict['I'] = '5-Helix'
ss_dict['T'] = 'Turn'
ss_dict['S'] = 'Bend'
ss_dict['C'] = 'Coil'

def classify_ss(ss):
    try:
        ss_category = ss_dict[ss]
    except KeyError:
        ss_category = ss_dict['C']
    return(ss_category)

df['Cys1_SS_cat'] = df['Cys1_SS'].apply(classify_ss)
df['Cys2_SS_cat'] = df['Cys2_SS'].apply(classify_ss)

import collections
import operator

############################
# Round all five cystine X angles to nearest 5 degress
#############################


config = df[['Cys1_x1','Cys1_x2','x3','Cys2_x2','Cys2_x1']]
config['Cys1_x1'] = config['Cys1_x1'].map(lambda Cys1_x1: int(5 * round(float(Cys1_x1)/5)))
config['Cys2_x1'] = config['Cys2_x1'].map(lambda Cys2_x1: int(5 * round(float(Cys2_x1)/5)))

config['Cys1_x2'] = config['Cys1_x2'].map(lambda Cys1_x2: int(5 * round(float(Cys1_x2)/5)))
config['Cys2_x2'] = config['Cys2_x2'].map(lambda Cys2_x2: int(5 * round(float(Cys2_x2)/5)))

config['x3'] = config['x3'].map(lambda x3: int(5 * round(float(x3)/5)))
print config['x3']
#config = config.loc[config['x3'].isin(range(0,180))]
#########################################
## Looking at total X angle distribution so combine X1 / X1B and X2 / X2B
#########################################
x1 = config['Cys1_x1'].tolist() + config['Cys2_x1'].tolist()
x2 = config['Cys1_x2'].tolist() + config['Cys2_x2'].tolist()
x3 = config['x3'].tolist()
##########################################################
## Use collections to counter the frequency of each angle
###########################################################

#### For range betwen -180,180, count frequency of number so that includes missing numbers##

x1 = [180 if x == -180 else x for x in x1]
x2 = [180 if x == -180 else x for x in x2]
x2 = [180 if x == -180 else x for x in x2]
print x1

x1_frequency = collections.Counter(x1)
x2_frequency = collections.Counter(x2)
x3_frequency = collections.Counter(x3)

x1_frequency[-180] = x1_frequency[180]
x2_frequency[-180] = x2_frequency[180]

########################################
# Sort dictionary in orderered fashion (smallest to largest)
###########################################
ordered_x1_frequency = sorted(x1_frequency.items(), key=operator.itemgetter(0))
ordered_x2_frequency = sorted(x2_frequency.items(), key=operator.itemgetter(0))
ordered_x3_frequency = sorted(x3_frequency.items(), key=operator.itemgetter(0))

print ordered_x1_frequency

################################
# Define each axis
################################
x1_axis = [(_[0]) for _ in ordered_x1_frequency]
x2_axis = [(_[0]) for _ in ordered_x2_frequency]
x3_axis = [(_[0]) for _ in ordered_x3_frequency]

#########################################################################
# Convert the number of angles to a frequency by dividing by total angles
#########################################################################
x1_frequency = [(float(_[1]) / float(len(x1))) for _ in ordered_x1_frequency]
x2_frequency = [(float(_[1]) / float(len(x2))) for _ in ordered_x2_frequency]
x3_frequency = [(float(_[1]) / float(len(x3))) for _ in ordered_x3_frequency]



##########################################
# Start to generate line graphs'
#########################################
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np



#######################################
# Create X1 graph
#######################################
plt.plot(x1_axis,x1_frequency, color = 'darkblue',linewidth = 2)

ax = plt.subplot(111, projection='polar')
x1_axis = [x+180 for x in x1_axis]
x2_axis = [x+180 for x in x2_axis]
x3_axis = [x+180 for x in x3_axis]
import numpy as np
x1_axis = [np.deg2rad(x) for x in x1_axis]
x2_axis = [np.deg2rad(x) for x in x2_axis]
x3_axis = [np.deg2rad(x) for x in x3_axis]

print x1_axis
print x1_frequency
ax.plot(x1_axis, x1_frequency, linewidth =2, color = 'blue',label=r"$\chi$1")
ax.plot(x2_axis, x2_frequency,linewidth =2, color = 'red',label=r"$\chi$2")
ax.plot(x3_axis, x3_frequency,linewidth =2, color = 'green',label=r"$\chi$3")
ax.set_rmax(0.14)
# ax.rlabels([0.02,0.06,0.10,0.14])
ax.set_yticklabels(['','',0.04,'','0.08','','0.12'])
# less radial ticks
ax.set_xticklabels([r"180$\degree$", r"-135$\degree$",  r"-90$\degree$", r"-45$\degree$",  r"0$\degree$",  r"45$\degree$",  r"90$\degree$", r"135$\degree$"])
ax.set_rlabel_position(190)  # get radial labels away from plotted line
ax.grid(True)
box = ax.get_position()
ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5))

ax.set_title("Frequency of cystine dihedral angles", va='bottom')
plt.show()
