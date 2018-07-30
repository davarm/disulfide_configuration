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
        
#df['Cys1_x1']  = df['Cys1_x1'].apply(x1_rounded)
#df['Cys2_x1']  = df['Cys2_x1'].apply(x1_rounded)
#df['Cys1_x2']  = df['Cys1_x2'].apply(x2_rounded)
#df['Cys2_x2']  = df['Cys2_x2'].apply(x2_rounded)
#df['x3'     ]  = df['x3'    ].apply(x3_rounded)

###########################################################################################
# Count the number of each configuration
# Have to consider assymetical configuraitons (x1,x2,x3,x2',x1' ==  x1',x2',x3,x2,1)
##############################################################################################

#def configuration_count(config):
#    config         = list(config)
#    forward_config = len(df.loc[(df['Cys1_x1'] == float(config[0])) & (df['Cys1_x2'] == float(config[1]))  & (df['x3'] == float(config#[2])) & (df['Cys2_x2'] == float(config[3]))& (df['Cys2_x1'] == float(config[4])) ]) #& df['x2'] == float(config[1])])
#    
#    ###########################################################################################
#    # If configuration IS NOT symmetrical, have to consider the reverse order (Cys2 - Cys1)
#    #############################################################################################
#    if config != config[::-1]:
#        reverse_config =  len(df.loc[(df['Cys2_x1'] == float(config[0])) & (df['Cys2_x2'] == float(config[1]))  & (df['x3'] == float(#config[2])) & (df['Cys1_x2'] == float(config[3]))& (df['Cys1_x1'] == float(config[4])) ])
#        forward_config = reverse_config + forward_config
#    return(forward_config)
#
######################################################################
## Separate configurations based on the x3 angle (either +90 or -90)
## Store in dictionary
#########################################################################
#plus_x3_configuration_dict  = {}
#minus_x3_configuration_dict = {}
#for config in configurations:
#    if list(config)[2] == 90:
#        plus_x3_configuration_dict[config] = configuration_count(config)
#    if list(config)[2] == -90:
#        minus_x3_configuration_dict[config] = configuration_count(config)
#        
#plus_x3_ordered_configuration  = sorted(plus_x3_configuration_dict.items(),  key=operator.itemgetter(0))
#minus_x3_ordered_configuration = sorted(minus_x3_configuration_dict.items(), key=operator.itemgetter(0))
#
#total_lefthanded = sum((item) for item in minus_x3_configuration_dict.values())
#total_righhanded = sum((item) for item in plus_x3_configuration_dict.values())
#total_disulfides = total_lefthanded + total_righhanded
#print 'Total Number of Disulfides:',len(df)
#print 'Total Number of Disulfides in Defined Configurations for Structural Analysis:',total_disulfides
#
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

#configuration_dataframe_dict = {}
#def configuration_dataframe_return(config):
#    config = list(config)
#    forward_config = pd.DataFrame(columns = list(df))
#    reverse_config = pd.DataFrame(columns = list(df))
#    ##############################################################################################################
#    # Forward Configuration:
#    # Dihedrals are labelled as X1,X2,X3,X2b,X1b:
#    # The forward configuration is when this order matches the configuration
#    # For example: Config == -180,-60,-90,-60,-60: 
#    # Search the database and get any cystine residues where we observe that configuration
#    # However the reverse of the CONFIG: -60,-60,-90,-60,-180 is also the same:
#    # Therefore have to search where X1b, X2B, X3, X2, X1 == Config
#    #
#    # IF THE CONFIGURATION IS IDENTICAL: -60,-60,-90,-60,-60 we do not search for the reverse configuration as then
#    # the dataframe would be doubled
#    ###############################################################################################################
#    forward_config = (df.loc[
#                        (df['Cys1_x1'] == float(config[0])) & 
#                        (df['Cys1_x2'] == float(config[1])) & 
#                        (df['x3'     ] == float(config[2])) & 
#                        (df['Cys2_x2'] == float(config[3])) & 
#                        (df['Cys2_x1'] == float(config[4])) 
#                       ])
#    
#    ###################################################################################
#    # If the configuration IS NOT symmetrical, then search for the reverse configuration
#    # x1b,x2b,x3,x2,x1
#    ###################################################################################
#    if config != config[::-1]:
#        reverse_config = (df.loc[
#                             (df['Cys2_x1'] == float(config[0])) & 
#                             (df['Cys2_x2'] == float(config[1])) & 
#                             (df['x3'     ] == float(config[2])) & 
#                             (df['Cys1_x2'] == float(config[3])) & 
#                             (df['Cys1_x1'] == float(config[4])) 
#                            ])
#        ###############################################################################################
#        # Apply the rewrite_columns function so Cys1 and Cys2 properties match with that of forward_config
#        ################################################################################################
#        reverse_config = rewrite_columns(reverse_config)
#
#        #################################################################################
#        # Append the two dataframes together
#        #################################################################################
#        forward_config = pd.concat([forward_config, reverse_config])
#    
#    return(forward_config)
#
#for config in configurations:
#    configuration_dataframe_dict[config] = configuration_dataframe_return(config)


#################################################################################
# TAKE SECONDARY STRUCTURE:::::
#################################################################################

from collections import Counter
import itertools


secondary_structure = [
'a-Helix',
'B-Bridge',
'Strand',
'3-Helix',
'5-Helix',
'Turn',
'Bend',
'Coil']

ss_possible = itertools.combinations(secondary_structure,2)
unique_ss = []
for _ in ss_possible:
   if  _[::-1] not in unique_ss:
        unique_ss.append(_)
from collections import OrderedDict

unique_ss.append(('a-Helix','a-Helix'))
unique_ss.append(('B-Bridge','B-Bridge'))
unique_ss.append(('Strand','Strand'))
unique_ss.append(('3-Helix','3-Helix'))
unique_ss.append(('5-Helix','5-Helix'))
unique_ss.append(('Turn','Turn'))
unique_ss.append(('Bend','Bend'))
unique_ss.append(('Coil','Coil'))

def generate_circle_plot(x1,x2,x3,ss):
    print 'xxxxxxxxxxxxxxxxx'
    print x1
    x1 = [180 if x == -180 else x for x in x1]
    x2 = [180 if x == -180 else x for x in x2]
    x2 = [180 if x == -180 else x for x in x2]
    
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
    ax.set_rmax(0.2)
    # ax.rlabels([0.02,0.06,0.10,0.14])
    ax.set_yticklabels(['','',0.05,'','0.1','','0.15','','0.20'])
    # ax.set_yticklabels(['','',0.1,'','0.08','','0.12'])
    # less radial ticks
    ax.set_xticklabels([r"180$\degree$", r"-135$\degree$",  r"-90$\degree$", r"-45$\degree$",  r"0$\degree$",  r"   45$\degree$",  r"90$\degree$", r"135$\degree$"])
    ax.set_rlabel_position(190)  # get radial labels away from plotted line
    ax.grid(True)
    box = ax.get_position()
    ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    
    ax.set_title(ss, va='bottom')
    plt.savefig(str(ss)+'.png', dpi=100, bbox_inches='tight')
    plt.close()
    # plt.show()
    return()







df['Cys1_x1'] = df['Cys1_x1'].map(lambda Cys1_x1: int(5 * round(float(Cys1_x1)/5)))
df['Cys2_x1'] = df['Cys2_x1'].map(lambda Cys2_x1: int(5 * round(float(Cys2_x1)/5)))
df['Cys1_x2'] = df['Cys1_x2'].map(lambda Cys1_x2: int(5 * round(float(Cys1_x2)/5)))
df['Cys2_x2'] = df['Cys2_x2'].map(lambda Cys2_x2: int(5 * round(float(Cys2_x2)/5)))
df['x3'] = df['x3'].map(lambda x3: int(5 * round(float(x3)/5)))

df['combined'] = df[['Cys1_SS_cat','Cys2_SS_cat']].apply(lambda x: ','.join(x), axis=1)
for ss in unique_ss:
    ss_config = []
    x1_list = []
    x2_list = []
    x3_list = []
    ss_forward         = ss[0]+','+ss[1]
    ss_reverse = ss[1]+','+ss[0]
    print ss, ss_reverse
    ss_config = df.loc[(df['combined'] == ss_forward) | (df['combined'] == ss_reverse) ]

    #x1_list.append(ss_config['Cys1_x1'].tolist())
    #x1_list.append(ss_config['Cys2_x1'].tolist())
    x1_list= (ss_config['Cys2_x1'].tolist()) +(ss_config['Cys2_x1'].tolist()) 
    x2_list= (ss_config['Cys2_x2'].tolist()) +(ss_config['Cys2_x2'].tolist()) 
    #x2_list.append(ss_config['Cys1_x2'].tolist())
    #x2_list.append(ss_config['Cys1_x2'].tolist())

    x3_list = (ss_config['x3'].tolist())

    if len(x3_list) > 100:
        generate_circle_plot(x1_list,x2_list,x3_list,ss)
    #x1_list = []
    #x2_list = []
    #x3_list = []
    #complete_total = 0
    #ss_forward         = ss[0]+','+ss[1]
    #ss_reverse = ss[1]+','+ss[0]
    ##print ss_forward
    #total_dict ={}
    #for config in configurations:
    #    # print config
    #    config_dataframe = configuration_dataframe_return(config)
    #    ss_list = config_dataframe[['Cys1_SS_cat','Cys2_SS_cat','Cys1_x1','Cys1_x2','x3','Cys2_x1','Cys2_x2']]
    #    ss_list['combined'] = ss_list[['Cys1_SS_cat','Cys2_SS_cat']].apply(lambda x: ','.join(x), axis=1)
    #    # print ss_list
    #    if len(ss_list) > 100:
    #        #print 'YES'
    #        #total = ss_list['combined'].tolist().count(ss_forward) +ss_list['combined'].tolist().count(ss_reverse)
    #        #complete_total = complete_total+total
    #        #total_dict[config] = total
    #        #print config,ss_forward, (float(total)/float(len(ss_list))*100)
    #        # print config,ss_forward, float(total)/flen(ss)
    #total_dict=collections.OrderedDict(sorted(total_dict.items()))
    ## print ss,complete_total
    #if complete_total < 10:
    #   print ss
    #   #generate_ss_bar_graph(total_dict,complete_total,ss)
#
#fig, ax = plt.subplots()
#image = np.random.uniform(size=(10, 10))
#
#fig.subplots_adjust(hspace=0.3, wspace=0.05)

#def generate_ss_bar_graph(frequency_dict,config_total,config):
#  #y_pos = np.arange(len(ss_x_axis_top3))
#    #plt.xticks(y_pos,ss_x_axis, rotation = 90)
#    
#    i= 0
#    for key in frequency_dict:
#        plt.bar(i,(float(key[1])/(float(config_total))))
#        i =i+1
#    #for key in frequency_dict:
#    #    plt.bar(ss_x_axis_top3.index(key[0]),(float(key[1])/(float(config_total))), color = 'blue')
#    #plt.xticks(y_pos,ss_x_axis_top3, rotation = 90)
#    #plt.ylim(0,1)
#    plt.title(config)
#    plt.ylabel("Frequency")
#    plt.xlabel('Cys1-Cys2 Secondary Structure')
#    plt.show()

 




  #config_dataframe = configuration_dataframe_return(config)
  #ss_list = config_dataframe[['Cys1_SS_cat','Cys2_SS_cat']]
  #ss_list['combined'] = ss_list[['Cys1_SS_cat','Cys2_SS_cat']].apply(lambda x: ','.join(x), axis=1)
  #if len(ss_list) > 3000:
  #    common_ss = Counter(ss_list['combined'].tolist())
  #    print common_ss
# #          common_ss = common_ss.most_common(3)
# #          generate_ss_bar_graph(common_ss, len(ss_list),config)
##            #