
path='all'
import itertools
import pylab
import os,sys
import shutil
import numpy as np
import math 
import matplotlib.pyplot as plt
amino_dict={}
aminos=['A','R','N','D','c','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
#x_axis=range(len(aminos))
z=len(aminos)
x_axis=(np.arange(0,z *2,2))
onex_axis=range(len(aminos))
twox_axis=range(len(aminos))
threex_axis=range(len(aminos))
fourx_axis=range(len(aminos))
for k,value in enumerate(x_axis):
    x_axis[k]=value+1
    onex_axis[k]=value
    twox_axis[k]=value+0.35
    threex_axis[k]=value+0.7
    fourx_axis[k]=value+1.05
for k,value in enumerate(aminos):
    amino_dict[aminos[k]]=k
resid_dict={}
for filename in os.listdir(path):
    x_label=[]
    xx=[]
    y=[]
    config=filename[:-5]
    config=config.strip('(')
    configa=config.replace(" ","")
    print configa
    first=[]
    second=[]
    after_count=[]
    third=[]
    fourth=[]
    get= open(os.path.join(path,'%s' %filename), "r")
    #bonds = sum(1 for _ in get)
    #rint bonds
   # if bonds >= 15:
    #    continue
    i=0
    for line in get:
          #  print line
            i=i+1
            first.append(line[0])
            second.append(line[2])
            third.append(line[4])
            fourth.append(line[6])
    
    xx= [[x,first.count(x)] for x in set(first)]
    after_count=[[x,second.count(x)] for x in set(second)]
    before_2_count=[[x,third.count(x)] for x in set(third)]
    after_2_count=[[x,fourth.count(x)] for x in set(fourth)]
    
    #print xx
        #print
   # print after_count 
    temp_dict={}
    after_dict={}
    before_2_dict={}
    after_2_dict={}
    
    present_aminos=[]
    after_aminos=[]
    before_2_aminos=[]
    after_2_aminos=[]
    
    for value in xx:
            temp_dict[value[0]]=value[1]
            present_aminos.append(value[0])
            
    for value in after_count:
            after_dict[value[0]]=value[1]
            after_aminos.append(value[0])
    
    for value in before_2_count:
            before_2_dict[value[0]]=value[1]
            before_2_aminos.append(value[0])
    
    for value in after_2_count:
            after_2_dict[value[0]]=value[1]
            after_2_aminos.append(value[0])
            
            
    y_values=[]
    after_y_values=[]
    before_2_values=[]
    after_2_values=[]
    
    
    for value in aminos:
        try:
            y_values.append(temp_dict[value])
        except KeyError:
            y_values.append(0)
    
    for value in aminos:
        try:
            after_y_values.append(after_dict[value])
        except KeyError:
            after_y_values.append(0)
            
    for value in aminos:
        try:
            before_2_values.append(before_2_dict[value])
        except KeyError:
            before_2_values.append(0)
    
    for value in aminos:
        try:
            after_2_values.append(after_2_dict[value])
        except KeyError:
            after_2_values.append(0)

    
    

    #print y_values
    if i>=15:
        """Convert to fraction"""
        before_total=sum(y_values)
        after_total=sum(after_y_values)
        
        before_2_total=sum(before_2_values)
        after_2_total=sum(after_2_values)
    
        before_fraction=[]
        after_fraction=[]
        
        before_2_fraction=[]
        after_2_fraction=[]
        
        
        for k,value in enumerate(y_values):
            before_fraction.append(float(float(value)/float(before_total)))
        
        for k,value in enumerate(after_y_values):
            after_fraction.append(float(float(value)/float(after_total)))
        
        for k,value in enumerate(before_2_values):
            before_2_fraction.append(float(float(value)/float(before_2_total))*-1)
       
        for k,value in enumerate(after_2_values):
            after_2_fraction.append(float(float(value)/float(after_2_total))*-1)
        
        width=1
        
        
        plt.figure(figsize=(20,10))
        ax=plt.subplot(111)
        plt.axhline(y=0.00,xmin=0,xmax=1, hold=None, color='black')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        
        plt.bar(onex_axis,before_fraction,width=0.35)
        plt.bar(twox_axis,after_fraction,width=0.35,color='red')
        
        plt.bar(threex_axis,before_2_fraction,width=0.35)
        plt.bar(fourx_axis,after_2_fraction,width=0.35,color='red')
        
        
        plt.xticks(x_axis, aminos, ha='center')
        axes=plt.gca()
        axes.set_ylim([-0.8,0.8])
        axes.set_xlim([0,len(aminos)*2])
        
        """Reverse axis"""
        #plt.gca().invert_yaxis()
        
        plt.show()
        #fig.close()
     #   $plt.savefig(configa+'.png')
     
  