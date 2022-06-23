# -*- coding: utf-8 -*-

import sys
sys.path.append('/cnrm/tropics/user/brientf/MESONH/scripts/')
from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
import tools2 as tl
import makefigures as mf
# open infocas
from infoplot_home import *

styles=['-','--',':','-.','-','--',':','-.','-','--',':','-.']

def main(path,infocas,var0):
    if __name__== "__main__" :
      pathout   = '/Users/florentbrient/Dropbox/GitHub/scripts_LES/'+'d_unif/'
      tl.mkdir(pathout)
      file0    = path+"{cas}/{simu}/{simu}.1.{vtype}.TTTT.nc4"

      fig = plt.figure()
      ax  = fig.add_subplot(111)
      names     = infocas.keys()

      for nn,name in enumerate(names):
        info  = infocas[name]
        print name
        for hh,hour in enumerate(info['hour']):
          nc = None
          if 'nc' in info.keys():
            nc = info['nc']   
          vtype = 'V0301'
          if 'vtyp' in info.keys():
            vtype = info['vtyp']
          name1 = '--'
          if 'name1' in info.keys():
            name1 = info['name1']

          # Open netcdf
          Data,filenc,cas,simu = tl.opennc(path,info,name,hour,nc=nc,vtype=vtype,rmsimu=True)
          ZZ  = Data['ZHAT']
          var = var0
          if 'svt' in info.keys():
            var  = tl.replacename(var0)
          print var

          if var in Data.variables.keys():   
            data = Data[var] #np.mean(,axis=(1,2,))
          else:
            data = tl.createnew(var,Data)
          if data is not None:
            data = np.mean(np.squeeze(data),axis=(1,2,))

          # Figure
          color = info['color']
          if data is not None:
            ax.plot(data,ZZ,linestyle=styles[hh],color=color,label=name)
            DZ     = 100
            ax.text(data[0],ZZ[0]-DZ,hour,color=color,fontsize=10)
            
      # remove duplicate in legend
      handles, labels = plt.gca().get_legend_handles_labels()
      by_label = OrderedDict(zip(labels, handles))
      print by_label
      plt.legend(by_label.values(), by_label.keys())

      plt.xlabel(var0)
      plt.ylabel('ZZ (km)')
      #plt.show()
      name   = 'prof_'+var0
      xsize  = (8,10)
      mf.savefig(fig,ax,pathout,name=name,xsize=xsize)

path   = "/Volumes/TOSHIBA/"
#var    = 'THLM' #'W_VAR' #'N' #'TA' #'SVT003' #'N' #'RCT'
vars   = ['NNV'] #['NNV','THLM','RNPM','THV','W_VAR','W2','RHO','THT','REHU','THS2','THS1','MSE','CLDFR','RCT','SVT001','SVT002','SVT003','WINDSHEAR','TA','N']
for var in vars:
  main(path,infocas,var)
