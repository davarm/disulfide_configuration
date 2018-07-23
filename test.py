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
def custom_round(x, base=5):
    x = int(base * round(float(x)/base))
    
    ###############################################################################
    # To took at the distribution of angles at 180 degrees: We classify -150 to -180 as 180 - 210
    ###############################################################################
    if x < -150:
        x = x + 360
    return x


config = df[['Cys1_x1','Cys1 x2','x3','Cys2 x2','Cys2 x1']]

config['Cys1_x1'] = config['Cys1_x1'].map(lambda Cys1_x1: int(5 * round(float(Cys1_x1)/5)))
config['Cys2_x1'] = config['Cys2_x1'].map(lambda Cys2_x1: int(5 * round(float(Cys2_x1)/5)))

config['Cys1_x2'] = config['Cys1_x2'].map(lambda Cys1_x2: int(5 * round(float(Cys1_x2)/5)))
config['Cys2_x2'] = config['Cys2_x2'].map(lambda Cys2_x2: int(5 * round(float(Cys2_x2)/5)))

config['x3'] = config['x3'].map(lambda x3: int(5 * round(float(x3)/5)))
print config