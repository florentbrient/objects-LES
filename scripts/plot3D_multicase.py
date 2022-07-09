#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 15:06:45 2022

@author: fbrient

# Goal:
# Plot figures for several cases
# Use of the 3D field
# Use of objects
# Plot:
    # Distribution of 3D field below zi
    # 
"""
#import sys
import numpy as np
from matplotlib import pyplot as plt
import tools as tl
import toolstats as ts
import makefigures as mf
from infoplot import *
import seaborn as sns
from scipy import stats
from scipy.stats import moment
from scipy.stats import skew
from scipy.stats import skewtest
from scipy.stats import shapiro
from copy import deepcopy
import math
from skimage import measure

def getnames(name):
    info  = infocas[name]
    print('name: ', name)
    nc = None
    if 'nc' in info.keys():
        nc = info['nc']   
    vtype = 'V0301'
    if 'vtyp' in info.keys():
        vtype = info['vtyp']
    # Figure
    color = info['color']
    style = styles[0]
    if 'style' in info.keys():
        style = info['style']
    sens = info['sens']
    case = name
    if 'cas' in info.keys():
        case = info['cas']
    return info,nc,vtype,color,style,sens,case

def maxbin(var):
    print('maxbin ',var)
    if var == 'WT':
        minb, maxb, binwidth = [-3., 3., 0.1]
    return minb, maxb, binwidth

def orgaindex(field,name,gridsize=1,vmin2D=0,connectivity=2):
    OBJ           = measure.label(field, connectivity=connectivity)
    #objmask       = tl.do_unique(deepcopy(OBJ)) # use if delete used
    nbobj         = len(np.unique(OBJ))-1
    tmp           = np.nan
    #print(nbobj)
    if nbobj>0:
      objmask       = tl.do_unique(deepcopy(OBJ))*gridsize
      OBJ,nbr       = tl.do_delete2(OBJ,objmask,vmin2D*gridsize,rename=True)
      
      OBJ           = np.ma.masked_array(OBJ, mask=(OBJ==0))
      cloudprop     = measure.regionprops(OBJ)
      if name == 'eqdia':
        tmp         = np.asarray([prop.equivalent_diameter for prop in cloudprop])
        tmp         = tmp*gridsize
      elif name == 'ecc':
        tmp         = np.asarray([prop.eccentricity for prop in cloudprop])
      elif name == 'nd':
        # Minimum distance between objects (in pixels)
        mindist     = 0 #30
        if nbobj>2: # at least 3 objects
          tmp         = ts.neighbor_distance(cloudprop,mindist=mindist)
      else:
        sys.exit('Problem orgaindex')
    return tmp,nbobj

def plotdistrib(ax,data_distrib,bins
                ,color='k',style='-'
                ,label=None,weights=None):
    if weights is not None:
        print(data_distrib, color, bins, label, weights)
        plt.hist(data_distrib, color = color, density = True \
                 , edgecolor = 'black', bins = bins\
                 , linestyle = style\
                 , label = label\
                 , weights=weights)
    else:
        ax = sns.distplot(data_distrib,bins=bins,
                          hist=False, kde=True, 
                          color = color, 
                          hist_kws={'edgecolor':'black'},
                          kde_kws={'linewidth': 2,'linestyle':style},
                          label=label
                          #palette=sns.color_palette('bright')[:3],
                          )
        #l1 = ax.lines[0]
        #x1 = l1.get_xydata()[:, 0]
        #y1 = l1.get_xydata()[:, 1]
    return ax
                
def plotindex(ax,data,alt\
              ,color='k',style='-',width=3\
              ,label=None):
    ax.plot(data,alt,\
            color=color,linestyle=style,linewidth=width,\
            label=label)
    return ax
    

styles=['-','--',':','-.','-','--',':','-.','-','--',':','-.']
# Dimensions
var1D  = ['vertical_levels','W_E_direction','S_N_direction'] #Z,Y,X

def main(dirin,dirout,plottyp,var,plotobj,rely=False):
    if __name__== "__main__" :
      fig   = plt.figure()
      fts   = 12
      ax    = fig.add_subplot(111)
      names = infocas.keys()
      
      namesall = '_'.join(names)
      
      #if nameplot=='distrib':
      #    distrib,indexplot  = 1,0
      #else:
      #    distrib,indexplot  = 0,1
            
      # For distrib plot
      if plottyp == 1:
        minb, maxb, binwidth = maxbin(var)   
        bins = np.arange(minb, maxb + binwidth, binwidth)
      
      for nn,name in enumerate(names):
        # get information about the simulation
        info,nc,vtype,color,style,sens,case = getnames(name)
        
        for hh,hour in enumerate(info['hour']):
          # Open netcdf
          Data,filenc,cas,simu = tl.opennc(dirin,info,name,hour,nc=nc,vtype=vtype)
          
          data1D,ALT,dz,nzyx,nxnynz,sizezyx = tl.gridsize(Data,var1D)
          gridsize     = (nzyx[var1D[1]]*nzyx[var1D[2]])**0.5
          print(nzyx)
          
          varread = var
          if plottyp==2:
              varread = "WT"
          
          if varread in Data.variables.keys():   
            data = Data[varread] #np.mean(,axis=(1,2,))
          else:
            data = tl.createnew(varread,Data,var1D)
            
          if data is not None:
             # We assume that all grid have the same size
             # Domain size of interest
             offset = 0.25
             idxzi,toppbl,grad = tl.findpbltop('THLM',Data,offset=offset)
             print('Altitude of idxzi : ',ALT[idxzi])
             data_distrib = data
             NL           = data_distrib.shape[0]
             
             data1D       = data_distrib[1:idxzi,:,:].flatten() # remove surface
                          
             ####### Statistics
             #k2, p = stats.normaltest(data_distrib)
             stat, p2 = shapiro(data1D)
             #print(k2,p)
             print(stat,p2)
             aa = [moment(data1D,moment=ij) for ij in range(5)]
             print(aa)
             sk = skew(data1D)
             print('Skewness: ',sk)
             # Test shewness
             #bb = skewtest(data_distrib)
             #print(bb)
             
             if plottyp==1:
                 relative = False
                 xx       = data1D
                 if rely:
#                     relative = True
                     std  = np.nanstd(data1D)
                     mean = np.nanmean(data1D)
                     xx   = (data1D-mean)/std
                 
                 label = "{0} {1} (s={2})".format(name,hour,'{:.2f}'.format(sk))
                 ax = plotdistrib(ax,xx,bins,\
                                  style=styles[hh],\
                                  label=label,color=color)
             if plottyp==2:
                nbobj = np.zeros(NL)
                 
             if plotobj:
                 # Creating mask with objects
                 thrs   = 2; thch   = str(thrs).zfill(2)
                 AddWT  = 1 # Vertical velocity
                 minchar = 'volume' #unit
                 if minchar == 'volume':
                     vmin   = 0.02 # by default
                 nbplus    = tl.svttyp(case,sens) #1
                 # Object based on tracer concentration AND vertical velocity?
                 typs, objtyp = tl.def_object(thrs,nbplus=nbplus,AddWT=AddWT)
                 # Remove returning shells
                 typs   = [ij for ij in typs if not (('down_SVT001' in ij) or ('down_SVT004' in ij))]
                 print(typs)
                 #stop                       
                 # Compute objects once
                 dataobjs = {}
                 for typ in typs:
                     nameobj = typ+'_'+thch
                     dataobj = Data[nameobj][:] #tl.removebounds(DATA[nameobj])
                     # new version (volume)
                     tmpmask                = tl.do_unique(deepcopy(dataobj))*nxnynz
                     dataobjs[nameobj],nbr  = tl.do_delete2(dataobj,tmpmask,vmin,rename=True)
                     zm                     = deepcopy(dataobjs[nameobj])
                     zm[zm>0]               = 1
                     zm                     = np.ma.masked_where(zm==0,zm)
                     colorobj               = tl.findhatch(nameobj,typ='color')

                     if plottyp==1:
                        weights             = zm[1:idxzi,:,:].flatten() 
                        ax = plotdistrib(ax,data1D,bins,\
                                         color=colorobj,style=style,\
                                         weights=weights)
                     if plottyp==2:
                         result = np.zeros(NL)
                         # remove smaller objects than xx grid (xx=4)
                         xx           = 100 #10 #4
                         vmin2D       = xx #*gridsize**2. # grid size
                         for ij in range(NL):
                             field        = zm[ij,:,:]
                             #objmask      = tl.do_unique(deepcopy(field))*nxnynz
                             #field,nbr    = tl.do_delete2(field,objmask,vmin2D,rename=True)
                             #print(ij,field)
                             connectivity = 2
                             #tmp          = np.nan
                             #if nbr>0: 
                             tmp,nbobj    = orgaindex(field,var,gridsize=gridsize,
                                                      vmin2D=vmin2D,connectivity=connectivity)
                             result[ij]    = np.nanmean(tmp) #*(gridsize*1000.) # in meters
                         #print(NL,ALT,index)
                         
                         ZZ = ALT
                         if rely: # relative to inversion
                             ZZ    = ALT/ALT[idxzi] 
                         ax = plotindex(ax,result,ZZ\
                                        ,color=colorobj,style=styles[nn]
                                        ,label=name)
      
      if plottyp == 1:            
          # plot ideal Gaussian
          mu = 0; variance = 1
          sigma = math.sqrt(variance)
          ax.plot(bins, stats.norm.pdf(bins, mu, sigma),'ko'\
                  ,markersize=3, mfc='none')
      else:
          yaxis = 'Altitude (km)'
          if rely:
            yaxis = 'z/zi'
            ax.axhline(y=1, color='k', linewidth=0.5, linestyle='--')

      # remove duplicate in legend
      handles, labels = plt.gca().get_legend_handles_labels()
      #by_label        = OrderedDict(zip(labels, handles))
      by_label        = dict(zip(labels, handles))
      print(by_label)
      plt.legend(by_label.values(), by_label.keys())
      #
      ax.axvline(x=0,color='k', linewidth=1, linestyle='--')
      ax.axhline(y=0,color='k', linewidth=1, linestyle='-')
      tl.adjust_spines(ax, ['left', 'bottom'])
      
      if plottyp == 1:
          namex = var
          namey = 'Probability (-)'
          xsize  = (8,4)
          namefig = 'distrib'
      else:
          namex = var
          namey = yaxis
          xsize  = (5,8)
          namefig = 'prof'
          if rely:
              zmax = 1.5
              ax.set_ylim([0,zmax])
          
      suffix = ''   
      if rely:
          suffix='_rel'
     
      ax.set_xlabel(namex,fontsize=fts)
      ax.set_ylabel(namey,fontsize=fts)
      title   = namefig+'_'+var+'_'+namesall+suffix
      mf.savefig(fig,ax,dirout,title=title,xsize=xsize)
      
      
      #
      #xlev=findxlev(var)
      #ax.set_xlim(xlev)
      #ax.set_ylim([0,zmax])
      #titlefig = titlefig.format(cases='_'.join(keys))
      #plt.legend(loc=2, prop={'size': fts-5})


dirin   = "../data/"
dirout  = "../data/d_unif/"

# Which plot?
distrib  = 1 # 
index    = 0

plottyp  = 1

# With objects
plotobj = False

# variables
vars   = ['WT',] #['NNV','THLM','RNPM','THV','W_VAR','W2','RHO','THT','REHU','THS2','THS1','MSE','CLDFR','RCT','SVT001','SVT002','SVT003','WINDSHEAR','TA','N']

# for index 
if index == 1:
    plottyp  = 2 #index
    vars = ['eqdia','ecc','nd']
#index   = 'ecc'
#index   = 'nd'

rely    = True # normalized Gauss if distrib

#if distrib:
#  nameplot   = 'distrib'
#else:
#  nameplot   = index

for var in vars:
  main(dirin,dirout,plottyp,var,plotobj,rely=rely)