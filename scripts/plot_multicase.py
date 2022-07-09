#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 15:06:45 2022

@author: fbrient

# Goal:
# Plot domaine-averaged figures for several cases

"""
#import sys
import numpy as np
from matplotlib import pyplot as plt
import tools as tl
import makefigures as mf
from infoplot import *


styles=['-','--',':','-.','-','--',':','-.','-','--',':','-.']
# Dimensions
var1D  = ['vertical_levels','W_E_direction','S_N_direction'] #Z,Y,X

def main(dirin,dirout,infocas,var0,relx,rely):
    if __name__== "__main__" :
#      file0    = dirin+"{cas}/{simu}/sel_{simu}.1.{vtype}.TTTT.nc4"
      
      fig   = plt.figure()
      ax    = fig.add_subplot(111)
      names = infocas.keys()
      yaxis = 'z (m)' 
      
      for nn,name in enumerate(names):
        info  = infocas[name]
        print('name: ', name)
        for hh,hour in enumerate(info['hour']):
          nc = None
          if 'nc' in info.keys():
            nc = info['nc']   
          vtype = 'V0301'
          if 'vtyp' in info.keys():
            vtype = info['vtyp']
#          name1 = '--'
#          if 'name1' in info.keys():
#            name1 = info['name1']
      
        # Open netcdf
          Data,filenc,cas,simu = tl.opennc(dirin,info,name,hour,nc=nc,vtype=vtype)
          ZZ  = Data['vertical_levels']
          var = var0
          if 'svt' in info.keys():
            var  = tl.replacename(var0)
          print(var)
          
          if var in Data.variables.keys():   
            data = Data[var] #np.mean(,axis=(1,2,))
          else:
            data = tl.createnew(var,Data,var1D)
          if data is not None:
            data = np.mean(np.squeeze(data),axis=(1,2,))

          # Figure
          color = info['color']
          if data is not None:
            if relx: # relative to the surface value
                data -= data[0]
            if rely: # relative to inversion
                offset = 0.25
                idxzi,toppbl,grad = tl.findpbltop('THLM',Data,offset=offset)
                print('idxzi: ', idxzi)
                ZZ    = ZZ/ZZ[idxzi] 
                yaxis = 'z/zi (m)' 
            ax.plot(data,ZZ,linestyle=styles[hh],color=color,label=name)
            DZ     = 100
            ax.text(data[0],ZZ[0]-DZ,hour,color=color,fontsize=10)
            
      if relx:  
          ax.axvline(x=0, color='k', linewidth=1)
      if rely:
          ax.axhline(y=1, color='k', linewidth=0.5, linestyle='--')
            
      # remove duplicate in legend
      handles, labels = plt.gca().get_legend_handles_labels()
      by_label        = OrderedDict(zip(labels, handles))
      print(by_label)
      plt.legend(by_label.values(), by_label.keys())

      plt.xlabel(var0)
      plt.ylabel(yaxis)
      #plt.show()
      title   = 'prof_'+var0
      if relx:
          title+='_relx'
      if rely:
          title+='_rely'
          zmax = 2
          ax.set_ylim([0,zmax])
      xsize  = (8,10)
      mf.savefig(fig,ax,dirout,title=title,xsize=xsize)


dirin   = "../data/"
dirout  = "../data/d_unif/"
relx    = True 
rely    = True
vars   = ['THV'] #['NNV','THLM','RNPM','THV','W_VAR','W2','RHO','THT','REHU','THS2','THS1','MSE','CLDFR','RCT','SVT001','SVT002','SVT003','WINDSHEAR','TA','N']
for var in vars:
  main(dirin,dirout,infocas,var,relx,rely)