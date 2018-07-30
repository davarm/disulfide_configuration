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


from matplotlib.pyplot import cm
import matplotlib as mpl
# mpl.rcParams['figure.figsize'] = (20.0, 10.0)
aminos=['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','X']
# aminos=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','X','Z']
amino_dict ={}
"""Assigning amino acid a number (1 to 20)"""
for k,value in enumerate(aminos):
    amino_dict[value]=k
# amino_dict['C']=amino_dict['c']



########################################################
# Frequency of Amino Acids By Biochemical Property
###########################################################
from matplotlib.pyplot import cm
import matplotlib as mpl
aminos=['A','R','N','D','c','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','X']
aa_group_dict = {}

aa_group_dict = (dict.fromkeys(['D', 'E'], 'Acid'))
aa_group_dict.update(dict.fromkeys(['K', 'R','H'], 'Base'))
aa_group_dict.update(dict.fromkeys(['P'], 'Pro'))
aa_group_dict.update(dict.fromkeys(['G'], 'Gly'))
aa_group_dict.update(dict.fromkeys(['F', 'W','Y'], 'Hydro.'))
aa_group_dict.update(dict.fromkeys(['S', 'T','C','c','N','Q'], 'Polar'))
aa_group_dict.update(dict.fromkeys(['A', 'V','L','M','I'], 'Aliphatic'))
aa_group_dict.update(dict.fromkeys(['X'], 'Ter'))

group_list = []
for k,value in enumerate(aa_group_dict.values()):
    if value not in group_list:
        group_list.append(value)
amino_dict = {}
for k,value in enumerate(group_list):
    amino_dict[value]=k

def aa_group(aa):
    aa_biochem_group = aa_group_dict[aa]
    return aa_biochem_group

"""Assigning amino acid a number (1 to 20)"""

df['Cys1 b res group'] = df['Cys1 b res'].apply(aa_group)
df['Cys2 b res group'] = df['Cys2 b res'].apply(aa_group)
df['Cys1 a res group'] = df['Cys1 a res'].apply(aa_group)
df['Cys2 a res group'] = df['Cys2 a res'].apply(aa_group)   

def amino_heatmaps(amino_array,amino_array_cys2,config):
        # plt.figure(figsize=(18, 6)) 
        fig = plt.subplot(1,1,1)
        
        amino_array=amino_array/np.sum(amino_array)
        amino_array_cys2=amino_array_cys2/np.sum(amino_array_cys2)
        
        for axis in [fig.xaxis, fig.yaxis]:
            axis.set(ticks=np.arange(0.5, len(group_list)), ticklabels=group_list)
        
        
        im = fig.pcolor(amino_array, cmap='nipy_spectral', edgecolor='black', linestyle=':', lw=1)
        im.set_clim(vmin=0,vmax =1)
        plt.ylabel("Cys1(i-1)")
        plt.xlabel("Cys1(i+1)")
        plt.title("Cys1 "+str(config))
        #plt.colorbar(im, orientation ='horizontal')
        for axis in [fig.xaxis, fig.yaxis]:
            axis.set(ticks=np.arange(0.5, len(group_list)), ticklabels=group_list)
        plt.savefig(str(config)+'cys1_amino_properties.png', dpi=300, bbox_inches='tight')
        #plt.figure(figsize=(18, 8)) 
        

        fig = plt.subplot(1,1,1)
        im = fig.pcolor(amino_array_cys2, cmap='nipy_spectral', edgecolor='black', linestyle=':', lw=1)
        im.set_clim(vmin=0,vmax =0.5)
        # plt.colorbar(im)
        
        print config
        
        for axis in [fig.xaxis, fig.yaxis]:
            axis.set(ticks=np.arange(0.5, len(group_list)), ticklabels=group_list)
        plt.ylabel("Cys2(i-1)")
        plt.xlabel("Cys2(i+1)")
        plt.title("Cys2 "+str(config))
        plt.savefig(str(config)+'cys2_amino_properties.png', dpi=300, bbox_inches='tight')
        #plt.show()
        return()
#


def search_config(config):   
    config_dataframe = configuration_dataframe_return(config)
    if len(config_dataframe)> 100 :
        amino_array = np.zeros(shape=(8,8))
        amino_array_cys2 = np.zeros(shape=(8,8))
        for index,row in config_dataframe.iterrows():
            cys1_b = row['Cys1 b res group']
            cys1_a = row['Cys1 a res group']
            cys2_b = row['Cys2 b res group']
            cys2_a = row['Cys2 a res group']
            amino_array[amino_dict[cys1_b],amino_dict[cys1_a]] = amino_array[amino_dict[cys1_b],amino_dict[cys1_a]]+1
            amino_array_cys2[amino_dict[cys2_b],amino_dict[cys2_a]] = amino_array_cys2[amino_dict[cys2_b],amino_dict[cys2_a]]+1
        amino_heatmaps(amino_array,amino_array_cys2,config)
    return()

for config in configurations:
    print config
    search_config(config)  