import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
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

x1_frequency = collections.Counter(x1)
x2_frequency = collections.Counter(x2)
x3_frequency = collections.Counter(x3)

########################################
# Sort dictionary in orderered fashion (smallest to largest)
###########################################
ordered_x1_frequency = sorted(x1_frequency.items(), key=operator.itemgetter(0))
ordered_x2_frequency = sorted(x2_frequency.items(), key=operator.itemgetter(0))
ordered_x3_frequency = sorted(x3_frequency.items(), key=operator.itemgetter(0))


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

import matplotlib.pyplot as plt

import numpy as np


import itertools
import os,sys
import shutil
import numpy
import math 

##########################################################
# GENERATE A LIST OF ALL 90 THEORETICAL CONFIGURATIONS
# X1 and X2 angles can be Gauche+, Gauche- or Trans
# X3 angles can be +90 or -90
# Over 90 theoretical configurations
# Store in a list as 'configurations'
#########################################################
a = 60
b = -60
c = 180
d = -90
e = 90

import string
exclude        = set(string.punctuation)
dihedralsx     = numpy.empty(shape=(5,1))
possible       = [[a,b,c],[a,b,c],[d,e],[a,b,c],[a,b,c]]
configurations = list(itertools.product(*possible))
unique         = []
for value in configurations:
      if value[::-1] in unique:
            continue
      unique.append(value)
configurations = unique

###############################################################
### Round all of the dihedral X angles to suitable conformation
# X1 angles are +/- 30 degrees of ,+60,-60 or 180
#                          +60 = +30 to +90
#                          -60 = -90 to -30 
#                          180 = 150 to 210
#
# X2 angles show more distribution (refer to previous figure) Therefore are:
#                          +60 = +30 to + 120
#                          -60 = -120 to -30
#                          180 =  150 to 210
#
# X3 angles are either +90 (Right Handed) or -90 (Left Handed)
#                          +90 = +60 to 120
#                          -90 = -120 to -60
################################################################

# df = df.loc[df['seqlength'] > 200]
dihedrals = df[['Cys1_x1','Cys1_x2','x3','Cys2_x2','Cys1_x1']]

#################################################################
# Define the rounding functions and apply to angles
###################################################################
def x1_rounded(x1):
    x1 = float(x1)
    if (x1 <=  90) & (x1 >= 30):
        x1= 60
    if (x1 >= -90)  & (x1 <= -30):
        x1 = -60   
    if (x1 <=  180) & (x1 >= 150):
        x1= 180
    if (x1 >= -180) & (x1 <= -150):
        x1=180
    return(x1)

def x2_rounded(x2):
    x2 = float(x2)
    if (x2 <=  130) & (x2 >= 30):
        x2= 60
    if (x2 >= -130) & (x2 <= -30):
        x2 = -60   
    if (x2 <=  180) & (x2 >= 150):
        x2= 180
    if (x2 >= -180) & (x2 <= -150):
        x2=180
    return(x2)

def x3_rounded(x3):
    x3 = float(x3)
    if (x3  <=  130)  & (x3 >= 60):
        x3   = 90
    if (x3  >=  -130) & (x3 <= -60):
        x3   = -90
    return(x3)
        
df['Cys1_x1']  = df['Cys1_x1'].apply(x1_rounded)
df['Cys2_x1']  = df['Cys2_x1'].apply(x1_rounded)
df['Cys1_x2']  = df['Cys1_x2'].apply(x2_rounded)
df['Cys2_x2']  = df['Cys2_x2'].apply(x2_rounded)
df['x3'     ]  = df['x3'    ].apply(x3_rounded)

###########################################################################################
# Count the number of each configuration
# Have to consider assymetical configuraitons (x1,x2,x3,x2',x1' ==  x1',x2',x3,x2,1)
##############################################################################################

def configuration_count(config):
    config         = list(config)
    forward_config = len(df.loc[(df['Cys1_x1'] == float(config[0])) & (df['Cys1_x2'] == float(config[1]))  & (df['x3'] == float(config[2])) & (df['Cys2_x2'] == float(config[3]))& (df['Cys2_x1'] == float(config[4])) ]) #& df['x2'] == float(config[1])])
    
    ###########################################################################################
    # If configuration IS NOT symmetrical, have to consider the reverse order (Cys2 - Cys1)
    #############################################################################################
    if config != config[::-1]:
        reverse_config =  len(df.loc[(df['Cys2_x1'] == float(config[0])) & (df['Cys2_x2'] == float(config[1]))  & (df['x3'] == float(config[2])) & (df['Cys1_x2'] == float(config[3]))& (df['Cys1_x1'] == float(config[4])) ])
        forward_config = reverse_config + forward_config
    return(forward_config)

#####################################################################
# Separate configurations based on the x3 angle (either +90 or -90)
# Store in dictionary
########################################################################
plus_x3_configuration_dict  = {}
minus_x3_configuration_dict = {}
for config in configurations:
    if list(config)[2] == 90:
        plus_x3_configuration_dict[config] = configuration_count(config)
    if list(config)[2] == -90:
        minus_x3_configuration_dict[config] = configuration_count(config)
        
plus_x3_ordered_configuration  = sorted(plus_x3_configuration_dict.items(),  key=operator.itemgetter(0))
minus_x3_ordered_configuration = sorted(minus_x3_configuration_dict.items(), key=operator.itemgetter(0))

total_lefthanded = sum((item) for item in minus_x3_configuration_dict.values())
total_righhanded = sum((item) for item in plus_x3_configuration_dict.values())
total_disulfides = total_lefthanded + total_righhanded
print 'Total Number of Disulfides:',len(df)
print 'Total Number of Disulfides in Defined Configurations for Structural Analysis:',total_disulfides

##########################################################################################
# The re_write columns is for the 'reverse configuraiton' dataframes
# where the config = x1b,x2b,x3,x2,x1
# By reversing the column names can easily then append to the forward dataframe for analysis
############################################################################################
def rewrite_columns(dataframe):
    dataframe.columns = dataframe.columns.str.replace("Cys2", "Cys3")
    dataframe.columns = dataframe.columns.str.replace("Cys1", "Cys2")
    dataframe.columns = dataframe.columns.str.replace("Cys3", "Cys1")
    return (dataframe)

###################################################################
# Function to return the pandas dataframe for the desired configuration
####################################################################

configuration_dataframe_dict = {}
def configuration_dataframe_return(config):
    config = list(config)
    forward_config = pd.DataFrame(columns = list(df))
    reverse_config = pd.DataFrame(columns = list(df))
    ##############################################################################################################
    # Forward Configuration:
    # Dihedrals are labelled as X1,X2,X3,X2b,X1b:
    # The forward configuration is when this order matches the configuration
    # For example: Config == -180,-60,-90,-60,-60: 
    # Search the database and get any cystine residues where we observe that configuration
    # However the reverse of the CONFIG: -60,-60,-90,-60,-180 is also the same:
    # Therefore have to search where X1b, X2B, X3, X2, X1 == Config
    #
    # IF THE CONFIGURATION IS IDENTICAL: -60,-60,-90,-60,-60 we do not search for the reverse configuration as then
    # the dataframe would be doubled
    ###############################################################################################################
    forward_config = (df.loc[
                        (df['Cys1_x1'] == float(config[0])) & 
                        (df['Cys1_x2'] == float(config[1])) & 
                        (df['x3'     ] == float(config[2])) & 
                        (df['Cys2_x2'] == float(config[3])) & 
                        (df['Cys2_x1'] == float(config[4])) 
                       ])
    
    ###################################################################################
    # If the configuration IS NOT symmetrical, then search for the reverse configuration
    # x1b,x2b,x3,x2,x1
    ###################################################################################
    if config != config[::-1]:
        reverse_config = (df.loc[
                             (df['Cys2_x1'] == float(config[0])) & 
                             (df['Cys2_x2'] == float(config[1])) & 
                             (df['x3'     ] == float(config[2])) & 
                             (df['Cys1_x2'] == float(config[3])) & 
                             (df['Cys1_x1'] == float(config[4])) 
                            ])
        ###############################################################################################
        # Apply the rewrite_columns function so Cys1 and Cys2 properties match with that of forward_config
        ################################################################################################
        reverse_config = rewrite_columns(reverse_config)

        #################################################################################
        # Append the two dataframes together
        #################################################################################
        forward_config = pd.concat([forward_config, reverse_config])
    
    return(forward_config)

for config in configurations:
    configuration_dataframe_dict[config] = configuration_dataframe_return(config)


#################################################################################
# TAKE SECONDARY STRUCTURE:::::
#################################################################################

from collections import Counter
import itertools

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
vdw_volume['X']= 5

import os
from bokeh.models import Jitter
from bokeh.plotting import figure, show, output_file
import math
# from bokeh.charts.utils import cycle_colors

#MATPLOTLIB
def generate_ss_bar_graph(cys1_b_vdw,cys1_a_vdw):
    fig,ax = plt.subplots()
    i=0
    x_axis_labels = []
    new_axis = []
    for key in cys1_b_vdw:
        b_vdw_mean = np.mean(cys1_b_vdw[key])
        b_vdw_std  = np.std(cys1_b_vdw[key])
        a_vdw_mean = np.mean(cys1_a_vdw[key])
        a_vdw_std  = np.std(cys1_a_vdw[key])
        x_axis_labels.append(key)
        new_axis.append(key)
        x_axis_labels.append('')

        if i ==0:
            ax.errorbar(i-0.2, b_vdw_mean, yerr=b_vdw_std, fmt='o', color = 'blue', capsize=3, ms =5, label = "i-1 residue")
            ax.errorbar(i+0.2, a_vdw_mean, yerr=a_vdw_std, fmt='o', color = 'red',capsize=3, ms =5, label = "i+1 residue")
        if i !=0:
            ax.errorbar(i-0.2, b_vdw_mean, yerr=b_vdw_std, fmt='o', color = 'blue', capsize=3, ms =5)
            ax.errorbar(i+0.2, a_vdw_mean, yerr=a_vdw_std, fmt='o', color = 'red',capsize=3, ms =5, )
        #for _ in distance_dict[key]:
        #    plt.scatter(i, _,color = 'blue')
        i = i+2
    ax.legend(loc='best')
    plt.ylim(0,180)
    # xticks_pos = [x for x in range(0,len(x_axis_labels))]
    xticks_pos = (np.arange(0, len(x_axis_labels), step=2))
    print xticks_pos

    ax.set_xticks(xticks_pos)
    ax.set_xticklabels(new_axis,rotation=90)
    plt.ylabel("VdW radius")
    plt.xlabel("Configuration")
    plt.title("Cys1")
    plt.savefig('Cys1_vdw_radi.png', dpi=300, bbox_inches='tight')
    plt.show()


    return()

distance_dict = {}  
config_list = [] 
cys1_b_vdw = {}
cys1_a_vdw = {}
cys2_b_vdw = {}
cys2_a_vdw = {}
def vdw_radi(aa):
    vdw_radi_v = vdw_volume[aa]
    return[vdw_radi_v]

for config in configurations:
    config_dataframe = configuration_dataframe_return(config)
    # config_dataframe = config_dataframe.loc[config_dataframe['chain1'] == config_dataframe['chain2']]
#    # Test if dictionary works first
#    config_dataframe = configuration_dataframe_dict[config]
#    Cys1 b res
#    
    config_dataframe['Cys1 b res vdw'] = config_dataframe['Cys1 b res'].apply(vdw_radi)
    config_dataframe['Cys2 b res vdw'] = config_dataframe['Cys2 b res'].apply(vdw_radi)
    config_dataframe['Cys1 a res vdw'] = config_dataframe['Cys1 a res'].apply(vdw_radi)
    config_dataframe['Cys2 a res vdw'] = config_dataframe['Cys2 a res'].apply(vdw_radi)



    if len(config_dataframe) > 100: 
            config_list.append(str(config))
            cys1_b_vdw[config] = config_dataframe['Cys1 b res vdw'].tolist()
            cys1_a_vdw[config] = config_dataframe['Cys1 a res vdw'].tolist()
            cys2_b_vdw[config] = config_dataframe['Cys2 b res vdw'].tolist()
            cys2_a_vdw[config] = config_dataframe['Cys2 a res vdw'].tolist()
            #print config_dataframe
#config_list = ['1','2''1','2''1','2''1','2''1','2''1','2''1','2''1','2''1','2']
generate_ss_bar_graph(cys1_b_vdw,cys1_a_vdw)
# print distance_dict
#            
##            #