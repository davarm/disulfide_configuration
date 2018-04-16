import pandas as pd
import numpy as np
import collections
import operator
df = pd.read_csv('disulfides.csv')
print len(df)


import itertools

import os,sys
import shutil
import numpy
import math 

### GENERATE A LIST OF ALL 90 THEORETICAL CONFIGURATIONS
a= 60
b= -60
c= 180
d= -90
e= 90

import string
exclude = set(string.punctuation)
dihedralsx=numpy.empty(shape=(5,1))

possible=[[a,b,c],[a,b,c],[d,e],[a,b,c],[a,b,c]]
configurations=list(itertools.product(*possible))
unique = []
for value in configurations:
      if value[::-1] in unique:
            continue
      unique.append(value)
configurations = unique

#print configuration_dict
#
### Round all of the dihedral X angles to suitable conformation
dihedrals = df[['x1','x2','x3','x2b','x1b']]

def x1_rounded(x1):
    x1=float(x1)
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
    x2=float(x2)
    if (x2 <=  120) & (x2 >= 30):
        x2= 60
    if (x2 >= -120) & (x2 <= -30):
        x2 = -60   
    if (x2 <=  180) & (x2 >= 150):
        x2= 180
    if (x2 >= -180) & (x2 <= -150):
        x2=180
    return(x2)

def x3_rounded(x3):
    x3=float(x3)
    if (x3  <=  120)  & (x3 >= 60):
        x3   = 90
    if (x3  >=  -120) & (x3 <= -60):
        x3   = -90
    return(x3)
        
df['x1' ]  = df['x1'].apply(x1_rounded)
df['x1b']  = df['x1b'].apply(x1_rounded)
df['x2' ]  = df['x2'].apply(x2_rounded)
df['x2b']  = df['x2b'].apply(x2_rounded)
df['x3' ]  = df['x3'].apply(x3_rounded)

# Count the number of each configuration
# Have to consider assymetical configuraitons (x1,x2,x3,x2',x1' ==  x1',x2',x3,x2,1)

def configuration_count(config):
    config = list(config)
    forward_config = len(df.loc[(df['x1'] == float(config[0])) & (df['x2'] == float(config[1]))  & (df['x3'] == float(config[2])) & (df['x2b'] == float(config[3]))& (df['x1b'] == float(config[4])) ]) #& df['x2'] == float(config[1])])
    if config != config[::-1]:
        reverse_config =  len(df.loc[(df['x1b'] == float(config[0])) & (df['x2b'] == float(config[1]))  & (df['x3'] == float(config[2])) & (df['x2'] == float(config[3]))& (df['x1'] == float(config[4])) ])
        forward_config = reverse_config + forward_config
    
    return(forward_config)

# Store in dictionary
plus_x3_configuration_dict = {}
minus_x3_configuration_dict = {}
for config in configurations:
    if list(config)[2] == 90:
        plus_x3_configuration_dict[config] = configuration_count(config)
    if list(config)[2] == -90:
        minus_x3_configuration_dict[config] = configuration_count(config)
        
plus_x3_ordered_configuration = sorted(plus_x3_configuration_dict.items(), key=operator.itemgetter(0))
minus_x3_ordered_configuration = sorted(minus_x3_configuration_dict.items(), key=operator.itemgetter(0))
print plus_x3_ordered_configuration


import math
import os
import sys
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import PDB
from matplotlib import colors


rama_preferences = {
            "General": {
            "file": "./rama_data/pref_general.data",
            "cmap": colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
            "bounds": [0, 0.0005, 0.02, 1],
        },
           "Second": {
            "file": "./rama_data/pref_general.data",
            "cmap": colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
            "bounds": [0, 0.0005, 0.02, 1],
        }
    }
rama_pref_values = {}
for key, val in rama_preferences.items():
        rama_pref_values[key] = np.full((360, 360), 0, dtype=np.float64)
        with open(val["file"]) as fn:
            for line in fn:
                if not line.startswith("#"):
                    # Preference file has values for every second position only
                    rama_pref_values[key][int(float(line.split()[1])) + 180][int(float(line.split()[0])) + 180] = float(
                        line.split()[2])
                    rama_pref_values[key][int(float(line.split()[1])) + 179][int(float(line.split()[0])) + 179] = float(
                        line.split()[2])
                    rama_pref_values[key][int(float(line.split()[1])) + 179][int(float(line.split()[0])) + 180] = float(
                        line.split()[2])
                    rama_pref_values[key][int(float(line.split()[1])) + 180][int(float(line.split()[0])) + 179] = float(
                        line.split()[2])
normals = {}
outliers = {}
for key, val in rama_preferences.items():
        normals[key] = {"x": [], "y": []}
        outliers[key] = {"x": [], "y": []}   
        
def ramachandran(cys1_phi,cys1_psi,cys2_phi,cys2_psi,config):
    """
    Main calculation and plotting definition
    :param file_name_list: List of PDB files to plot
    :return: Nothing
    """
    for idx, (key, val) in enumerate(sorted(rama_preferences.items(), key=lambda x: x[0].lower())):
        print idx
        plt.subplot(1, 2, idx + 1)
        #plt.subplot(1,2,1)
        plt.imshow(rama_pref_values[key], cmap=rama_preferences[key]["cmap"],
                   norm=colors.BoundaryNorm(rama_preferences[key]["bounds"], rama_preferences[key]["cmap"].N),
                   extent=(-180, 180, 180, -180))
        if idx == 0:
            plt.scatter(cys1_phi, cys1_psi,s = 6)
        if idx ==1:
            plt.scatter(cys2_phi, cys2_psi,s = 6, color = 'red')
        plt.xlim([-180, 180])
        plt.ylim([-180, 180])
        plt.plot([-180, 180], [0, 0], color="black")
        plt.plot([0, 0], [-180, 180], color="black")
        plt.locator_params(axis='x', nbins=7)
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$\psi$')
        plt.grid()
        plt.title(config,size=10)
        plt.tight_layout()
    plt.show()
        
        #plt.subplot(1,2,2)
        #plt.imshow(rama_pref_values[key], cmap=rama_preferences[key]["cmap"],
        #           norm=colors.BoundaryNorm(rama_preferences[key]["bounds"], rama_preferences[key]["cmap"].N),
        #           extent=(-180, 180, 180, -180))
    #
        #plt.scatter(cys2_phi, cys2_psi, color = 'red', s= 6)
        #plt.xlim([-180, 180])
        #plt.ylim([-180, 180])
        #plt.plot([-180, 180], [0, 0], color="black")
        #plt.plot([0, 0], [-180, 180], color="black")
        #plt.locator_params(axis='x', nbins=7)
        #plt.xlabel(r'$\phi$')
        #plt.ylabel(r'$\psi$')
        #plt.grid()
        #plt.title(config,size=10)
        #plt.tight_layout()
        #plt.show()
#
    
    ## Iterate though configurations:
def configuration_count(config):
    print config
    reverse_config = []
    forward_config = []
    config = list(config)
    cys1_phi_list = []
    cys2_phi_list = []
    cys1_psi_list = []
    cys2_psi_list = []   
    
    forward_config =    (df.loc[(df['x1'] == float(config[0])) & 
                             (df['x2'] == float(config[1])) & 
                             (df['x3'] == float(config[2])) & 
                             (df['x2b'] == float(config[3]))& 
                             (df['x1b'] == float(config[4])) 
                            ])
    

    if config == config[::-1]:
    
        cys1_phi_list = forward_config['phi'].tolist() 
        cys1_psi_list = forward_config['psi'].tolist() 
        cys2_phi_list = forward_config['phi_x'].tolist()
        cys2_psi_list = forward_config['psi_x'].tolist()


    
    if config != config[::-1]:
        reverse_config = (df.loc[(df['x1b'] == float(config[0])) & 
                             (df['x2b'] == float(config[1])) & 
                             (df['x3'] == float(config[2])) & 
                             (df['x2'] == float(config[3]))& 
                             (df['x1'] == float(config[4])) 
                            ])  
    
    
    	cys1_phi_list = forward_config['phi'].tolist() +  reverse_config['phi_x'].tolist()
    	cys1_psi_list = forward_config['psi'].tolist() +  reverse_config['psi_x'].tolist()
    	
    	cys2_phi_list = forward_config['phi_x'].tolist() +  reverse_config['phi'].tolist()
    	cys2_psi_list = forward_config['psi_x'].tolist() +  reverse_config['psi'].tolist()
    

        
        
    if len(cys1_phi_list)> 200:
        ramachandran(cys1_phi_list,
                     cys1_psi_list,
                     cys2_phi_list,
                     cys2_psi_list,
                    config)

    return()


for config in configurations:
      config,(configuration_count(config))