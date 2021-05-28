import numpy as np
import tools as tl

# Informations about the studies
# Ready for FIRE

def profref(var,z):
    tmp = np.ones(len(z))
    if var == 'U' or var == 'V':
      alpha = 305.; Ugeo = 6.0 
      if var == 'U':
         tmp   = tmp*Ugeo*np.cos(alpha*np.pi/180.)
      if var == 'V':
         tmp   = tmp*Ugeo*np.sin(alpha*np.pi/180.)
    elif var == 'THL' or var == 'TH':
      Tmin = 287.5; zmin = 0.595 # kms
      for tt in np.arange(len(z)):
        if z[tt] > zmin:
          tmp[tt]=299.5+0.0075*(z[tt]-zmin)*1000.
        else:
          tmp[tt]=Tmin
    elif var == 'RV' or var == 'RT':
      qmax = 9.6; zmin = 0.595 # kms
      for tt in np.arange(len(z)):
        if z[tt] > zmin:
          tmp[tt]=(6.6-0.003*(z[tt]-zmin)*1000.)#/1000.
        else:
          tmp[tt]=qmax#/1000.
    else:
       tmp = None
    return tmp


def readfire():
   path0     = '/cnrm/tropics/user/brientf/FIRE/Stephan/'
   pathfire1 = path0+'EUROCS_Florent/'
   #obs
   fileobs1 = pathfire1+'sni.rest.KG.txt'
   fileobs2 = pathfire1+'sni.diu_KG.txt'
   fileobs3 = pathfire1+'sni.detrended.top+base.KG.txt'
   nameLES  = ['lewellen','moeng','sanchLES','lock','chlond','zanten']; nbLES=len(nameLES)
   data1    = tl.readvalues(fileobs1);
   data2    = tl.readvalues(fileobs2,typ='\n',over=2);
   data3    = tl.readvalues(fileobs3);
   dataLES1D= {}; dataLES2D= {}; 
   timing   = ['36.5','23.5']
   sets     = ['D','E','F']
   for ij in nameLES:
      fileles=pathfire1+ij+'.setC_KG.txt'
      dataLES1D[ij] = tl.readvalues(fileles)
   for it in timing:
      dataLES2D[it] = {}
      for iss in sets:
        dataLES2D[it][iss] = {}
        for ij in nameLES:
          fileles = pathfire1+ij+'.set'+iss+it+'_KG.txt'
          dataLES2D[it][iss][ij] = tl.readvalues(fileles)
          
   #print 'dataLES  ',dataLES['moeng'][1]
   #print 'data2 ',fileobs2,data2,data1
   return data1,data2,data3,dataLES1D,dataLES2D
    
