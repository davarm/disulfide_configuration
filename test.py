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

def vdw_radi(aa):
    vdw_radi_v = vdw_volume[aa]
    return[vdw_radi_v]

import pandas as pd
df = pd.read_csv('xx.csv')

df['Cys1 b res vdw'] = df['Cys1 b res'].apply(vdw_radi)
df['Cys2 b res vdw'] = df['Cys2 b res'].apply(vdw_radi)
df['Cys1 a res vdw'] = df['Cys1 a res'].apply(vdw_radi)
df['Cys2 a res vdw'] = df['Cys2 a res'].apply(vdw_radi)

df.to_csv('xx.csv',index=False)
# print distance_dict