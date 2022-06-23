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
import seaborn as sns
from scipy import stats
from scipy.stats import moment
from scipy.stats import skew
from scipy.stats import skewtest
from scipy.stats import shapiro
from copy import deepcopy



styles=['-','--',':','-.','-','--',':','-.','-','--',':','-.']
# Dimensions
var1D  = ['vertical_levels','W_E_direction','S_N_direction'] #Z,Y,X

def main(dirin,dirout,infocas,var0,plotobj):
    if __name__== "__main__" :
      fig   = plt.figure()
      ax    = fig.add_subplot(111)
      names = infocas.keys()
      
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
          # Open netcdf
          Data,filenc,cas,simu = tl.opennc(dirin,info,name,hour,nc=nc,vtype=vtype)
          
          data1D,ALT,dz,nzyx,nxnynz,sizezyx = tl.gridsize(Data,var1D)
          
          var = var0
          #if 'svt' in info.keys():
          #  var  = tl.replacename(var0)
          print(var)
          
          
          if var in Data.variables.keys():   
            data = Data[var] #np.mean(,axis=(1,2,))
          else:
            data = tl.createnew(var,Data,var1D)
          if data is not None:
            # We assume that all grid have the same size
            # Domain size of interest
            offset = 0.25
            idxzi,toppbl,grad = tl.findpbltop('THLM',Data,offset=offset)
            print('Altitude of idxzi : ',data1D[var1D[0]][idxzi])
            data_distrib = data[1:idxzi,:,:]
            data_distrib = data_distrib.flatten()
            
          if var == 'WT':
            minb, maxb, binwidth = [-3., 3., 0.1]
          bins = np.arange(minb, maxb + binwidth, binwidth)
            
          # Figure
          color = info['color']
          style = styles[0]
          if 'style' in info.keys():
              style = info['style']
          if data is not None:
             # Statistics
             #k2, p = stats.normaltest(data_distrib)
             stat, p2 = shapiro(data_distrib)
             #print(k2,p)
             print(stat,p2)
             aa = [moment(data_distrib,moment=ij) for ij in range(5)]
             print(aa)
             sk = skew(data_distrib)
             print('Skewness: ',sk)
             # Test shewness
             bb = skewtest(data_distrib)
             print(bb)
             
             if plotobj:
                 plt.hist(data_distrib, color = color , density = False \
                      , edgecolor = 'black', bins = bins)
             else:
                #sns.set_theme(style="ticks")
                label = "{0} {1} (s={2})".format(name,hour,'{:.2f}'.format(sk))
                ax = sns.distplot(data_distrib,bins=bins
                          ,hist=False, kde=True, 
                          color = color, 
                          hist_kws={'edgecolor':'black'},
                          kde_kws={'linewidth': 3,'linestyle':style},
                          label=label
                          #palette=sns.color_palette('bright')[:3],
                          )
                l1 = ax.lines[0]
                x1 = l1.get_xydata()[:, 0]
                y1 = l1.get_xydata()[:, 1]
                #plt.fill_between(x1, y1, color=color, alpha=0.3)

             #sns.displot(data_distrib,color=color, kde=True)

             if plotobj == True:
                 # Creating mask with objects
                 thrs   = 2
                 thch   = str(thrs).zfill(2)
                 minchar = 'volume' #unit
                 if minchar == 'volume':
                     vmin   = 0.02 # by default
                     #vmin   = 0
                     #suffixmin = '_vol'+str(vmin)
                 sens = info['sens']
                 case = name
                 if 'cas' in info.keys():
                     case = info['cas']
                 nbplus    = tl.svttyp(case,sens) #1
                 # Object based on tracer concentration AND vertical velocity?
                 AddWT  = 1
                 typs, objtyp = tl.def_object(thrs,nbplus=nbplus,AddWT=AddWT)
                 # Compute objects once
                 dataobjs = {}
                 #nameobjs = []
                 #mask     = []
                 for typ in typs:
                     print(typ)
                     nameobj   = typ+'_'+thch
                     #nameobjs += [nameobj] #updr_SVT001_WT_02
                     #try:
                     dataobj   = Data[nameobj][:] #tl.removebounds(DATA[nameobj])
                         # new version (volume)
                     tmpmask = tl.do_unique(deepcopy(dataobj))*nxnynz
                         #print('tmpmask ',tmpmask[0,:,0])
                     dataobjs[nameobj],nbr  = tl.do_delete2(dataobj,tmpmask,vmin,rename=True)
                     zm                     = deepcopy(dataobjs[nameobj][1:idxzi,:,:])
                     zm[zm>0]               = 1
                     zm                     = np.ma.masked_where(zm==0,zm)
                     #mask.append(mask0)
                     #except:
                     #    dataobjs[nameobj] = None
    
                     #print(data_distrib)
                     #datamask = np.ma.masked_where(zm.flatten()==0,data_distrib)
                     #print(datamask)
                     
                     colorobj    = tl.findhatch(nameobj,typ='color')
                     plt.hist(data_distrib, color = colorobj, density = False \
                          , edgecolor = 'black', bins = bins, weights=zm.flatten())
                     #sns.distplot(data_distrib, hist=True,#, kde=True, 
                     #         color = 'red', 
                     #         hist_kws={'edgecolor':'black', 'weights': zm.flatten()},
                     #         #kde_kws={'linewidth': 3},
                     #         bins = bins)

            
            
      # remove duplicate in legend
      handles, labels = plt.gca().get_legend_handles_labels()
      by_label        = OrderedDict(zip(labels, handles))
      print(by_label)
      plt.legend(by_label.values(), by_label.keys())

      plt.xlabel(var0)
      ax.axvline(x=0,color='k', linewidth=1, linestyle='--')
      tl.adjust_spines(ax, ['left', 'bottom'])
      #plt.ylabel(yaxis)
      title   = 'distrib_'+var0
      xsize  = (10,5)
      mf.savefig(fig,ax,dirout,title=title,xsize=xsize)


dirin   = "../data/"
dirout  = "../data/d_unif/"
plotobj = False
vars   = ['WT'] #['NNV','THLM','RNPM','THV','W_VAR','W2','RHO','THT','REHU','THS2','THS1','MSE','CLDFR','RCT','SVT001','SVT002','SVT003','WINDSHEAR','TA','N']
for var in vars:
  main(dirin,dirout,infocas,var,plotobj)