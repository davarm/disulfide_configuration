# -*- coding: utf-8 -*-
"""
Created on Thu May 21 13:52:42 2015

@author: clarkelab
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:20:34 2015

@author: clarkelab
"""

# -*- coding: utf-8 -*-
"""
Same as conformations script, just writing out the CB shifts instead of angles
"""




import itertools

import os,sys
import shutil
import numpy
import math 
a=60
b=-60
c=180
d=-90
e=90
get2=open('all.txt','w')
import string
exclude = set(string.punctuation)

path='all'


dihedralsx=numpy.empty(shape=(5,1))

possible=[[a,b,c],[a,b,c],[d,e],[a,b,c],[a,b,c]]
configurations=list(itertools.product(*possible))
unique = []
for value in configurations:
      if value[::-1] in unique:
            continue
      unique.append(value)
configurations = unique

if os.path.exists(path):
      shutil.rmtree(path)
os.makedirs(path)



for k in configurations:
    #dir_path = os.path.join(self.feed, self.address)
    
    get=open(os.path.join(path,'%s.txt' %str(k)),'w+')
   
    
    with open("angles.txt") as fd:
        for i, line in enumerate(fd):

            lines=line.strip(",").split(",")
            dihedrals=(lines[17:22])
            dihedrals = [float(i) for i in dihedrals]
            cb1=lines[9]
            cb2=lines[15]
            #shifts = [float(i) for i in shifts]
            pdb=str(lines[0:3])
            pdb = ''.join(ch for ch in pdb if ch not in exclude)
            pdb=pdb.replace(" ","")
            rb1=lines[37]
            ra1=lines[47]
            rb2=lines[57]
            ra2=lines[67]
            if len(ra2)==2:
                ra2=ra2[0]
            
            dihedralsx[0]=dihedrals[0]
            dihedralsx[1]=dihedrals[1]
            dihedralsx[2]=dihedrals[2]
            dihedralsx[3]=dihedrals[3]
            dihedralsx[4]=dihedrals[4]
            phi=lines[22]
            psi=lines[23]
           
            
            

            
            phix=lines[25]
            psix=lines[26]

            """ UPDATED GAUCHE RANGES SO X2 IS BETWEEN 180-120 AND 120-30"""     
            dihedralsx[2:3][(dihedralsx[2:3]  <=  120) & (dihedralsx[2:3] >= 60)] = 90
            dihedralsx[2:3][(dihedralsx[2:3]  >=  -120) & (dihedralsx[2:3] <= -60)] = -90
              
         
            """x2 angles"""
        
            dihedralsx[1:4:2][(dihedralsx[1:4:2]  <=  120) & (dihedralsx[1:4:2] >= 30)] = 60
        
            dihedralsx[1:4:2][(dihedralsx[1:4:2]  >= -120) & (dihedralsx[1:4:2] <= -30)] = -60
        
            dihedralsx[1:4:2][(dihedralsx[1:4:2]  <=  180) & (dihedralsx[1:4:2] > 150)] = 180
        
            dihedralsx[1:4:2][(dihedralsx[1:4:2]  >= -180) & (dihedralsx[1:4:2] < -150)] = 180
        
            """ x1 angles """
              
            dihedralsx[0:5:4][(dihedralsx[0:5:4]  <=  90) & (dihedralsx[0:5:4] >= 30)] = 60
        
            dihedralsx[0:5:4][(dihedralsx[0:5:4]  >= -90) & (dihedralsx[0:5:4] <= -30)] = -60
        
            dihedralsx[0:5:4][(dihedralsx[0:5:4]  <=  180) & (dihedralsx[0:5:4] >= 150)] = 180
        
            dihedralsx[0:5:4][(dihedralsx[0:5:4]  >= -180) & (dihedralsx[0:5:4] <= -150)] = 180
            


    
            if dihedralsx[0] == k[0]:
                if dihedralsx[1] == k[1]:
                    if dihedralsx[2] == k[2]:
                        if dihedralsx[3] == k[3]:
                            if dihedralsx[4] == k[4]:
                                
                                get.write(rb1+','+ra1+','+rb2+','+ra2)
                                get.write('\n')

            if dihedralsx[0] == dihedralsx[4] and dihedralsx[1] == dihedralsx[3]:
                continue

            if dihedralsx[4] == k[0]:
                if dihedralsx[3] == k[1]:
                    if dihedralsx[2] == k[2]:
                        if dihedralsx[1] == k[3]:
                            if dihedralsx[0] == k[4]:
                                #get.write(shifts)
                                get.write(rb2+','+ra2+','+rb1+','+ra1)
                                get.write('\n')
                                #get.write('\n')

for filename in os.listdir(path):
    with open(os.path.join(path,'%s' %filename), "r") as f:
            bonds = sum(1 for _ in f)
            if bonds != 0:
                config=filename[:-5]
                config=config.strip('(')
                configa=config.replace(" ","")
                configa=[]
                print configa,bonds
                #for line in f:
                 #   configa.append(line[0])
                #print [[x,configa.count(x)] for x in set(configa)]
                
                #get2.write('\n')
          
