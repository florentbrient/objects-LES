import netCDF4 as nc
import pylab as plt
import numpy as np
import sys
import time
from copy import deepcopy
from collections import OrderedDict

import makefigures as mf
import tools as tl
#import infocases as info
import Constants as CC


def main(data1D,data,plots,dataobjs,xy=None,zview=[],pathfig='.',nametitle='',avg=0,fluxes=0,idxzi2D=None):
    if __name__== "__main__" :
      print('ok')
      title='_'.join(data.keys())
      if len(dataobjs)>0:
        title += '_'.join(dataobjs.keys())

      # Hatch for objects
      hchcross = True
      hchzview = False

      #print(dataobjs.keys())
      if plots['cross']:
        mf.plot2D(data1D,data,dataobjs,xy=xy,pathfig=pathfig,nametitle=nametitle,avg=avg,hatchlog=hchcross,idxzi2D=idxzi2D)

      if plots['zview']:
        for zz in zview:
          mf.plot2D(data1D,data,dataobjs,zz=zview[zz],pathfig=pathfig,nametitle=nametitle+'_'+zz,avg=avg,hatchlog=hchzview)

      if plots['mean']:
        mf.plotmean(data1D,data,dataobjs=dataobjs,zz=zview,pathfig=pathfig,nametitle=nametitle,fluxes=fluxes)
      



if len(sys.argv)!=6 : 
   print("Usage : python",sys.argv[0], "object varsel ncfile")
   exit(1)

# information about the simulation
path    = sys.argv[1] #path
case    = sys.argv[2] #IHOP
sens    = sys.argv[3] #Ru0x0
prefix  = sys.argv[4] #006
vtype   = sys.argv[5] #V0301

##################################################
#    Variables of interest
#vars       = ['Reflectance','LWP','WT','THV','RNPM','THLM','DIVUV','REHU','WINDSHEAR','RCT','PRW'] #,'RNPM','RVT','THLM','RCT','THV','DIVUV','REHU','PRW','LWP','LCL','LFC','LNB']
vars       = ['SVT004','SVT005','SVT006']

#    With objects?
objectchar = 1

#    Fluxes?
fluxes     = 0 #WT by default
fluxchar   = 'WT'

#    Which plots?
plots0     = {'cross':0, 'zview':0, 'mean':1}

# Average over +/-xx grid points for cross section
avg        = 0
##################################################

### Source Path
path0   = path #'/cnrm/tropics/user/brientf/MESONH/'

### Path to store figures
pathfig = path0+'figures/';tl.mkdir(pathfig)
pathfig = pathfig+case+'/';tl.mkdir(pathfig)
pathfig = pathfig+sens+'/';tl.mkdir(pathfig)

### conditions to plot
boxch     = False
boxes,xy  = tl.cond2plot(case,prefix,boxch=boxch)

### Open file
#path    = path0

# Name
filename = 'sel_'+sens+".1."+vtype+"."+prefix
path     = path + case +'/' + sens+'/'
suffixnc = '.nc4'

ncsuff  = ['IHODC','IHOP5',]
if sens in ncsuff:
  suffixnc= '.nc'

file    = path + filename + suffixnc
DATA    = nc.Dataset(file,'r')

if objectchar:
  # Objects"
  thrs   = 2
  thch   = str(thrs).zfill(2)
  # Minimal volume for object detection vmin
  # Brient et. al 19 (GRL): 0.25 km^3
  # By default : 0.02 km^3
  
  # Select by?
  minchar = 'volume' #unit
  if minchar == 'volume':
      vmin   = 0.02
      suffixmin = '_vol'+str(vmin)
  elif minchar == 'unit':
      nbmin  = 100 #100 #1000
      suffixmin = '_nb'+str(nbmin)
  
  
  nbplus = tl.svttyp(case,sens) #1
  AddWT  = 1
  typs, objtyp = tl.def_object(thrs,nbplus=nbplus,AddWT=AddWT)
  
  if ~fluxes:
    vars += ['Frac']
else:
  typs = []; objtyp=[]

# Dimensions
var1D  = ['vertical_levels','W_E_direction','S_N_direction'] #Z,Y,X
data1D,nzyx,sizezyx = [OrderedDict() for ij in range(3)]
for ij in var1D:
  data1D[ij]  = DATA[ij][:]/1000. #km #tl.removebounds(DATA[ij][:])
  nzyx[ij]    = data1D[ij][1]-data1D[ij][0]
  sizezyx[ij] = len(data1D[ij])
  #print ij,data1D[ij]
#data1D['vertical_levels']/=1000.

nxny   = nzyx[var1D[1]]*nzyx[var1D[2]] #km^2
ALT    = data1D[var1D[0]]
dz     = [0.5*(ALT[ij+1]-ALT[ij-1]) for ij in range(1,len(ALT)-1)]
dz.insert(0,ALT[1]-ALT[0])
dz.insert(-1,ALT[-1]-ALT[-2])
nxnynz = np.array([nxny*ij for ij in dz]) # volume of each level
nxnynz = np.repeat(np.repeat(nxnynz[:, np.newaxis, np.newaxis]
                             , sizezyx[var1D[1]], axis=1), sizezyx[var1D[2]], axis=2)

print('1: ',nxnynz.shape)
print('2: ',nxnynz[:,0,0])
print('3: ',nxnynz[0,0,:])

# name figure
nametitle0='TTTT_MMMM_{var}{objch}{vtypch}_{suffix}'+'AAAA' # Cross_x1_y2_Mean_Var_updrdowndowndown_2_nb100_008


# Find cloud base, cloud middle and cloud top
subcloud,cloudbase,cloudmiddle,cloudtop,zb,zi=[0 for ij in range(6)]
zview0 = {}
try:
   rct  = DATA['RCT'][:] #tl.removebounds(DATA['RCT'])
   z    = data1D[var1D[0]]
   # Define as the first layer where RCT (ql) > epsilon (1e-6 by default)
   epsilon = 1e-6
   cloudbase,cloudmiddle,cloudtop,zb,zi = tl.cloudinfo(z,rct,epsilon)
   subcloud = int(round(cloudbase/2.))
   zview0.update(
        {'subcloud':subcloud
        ,'cloudbase':cloudbase
        ,'cloudmiddle':cloudmiddle
        ,'cloudtop':cloudtop}
             )
except:
   pass


# Find lcl
#try:
Q           = DATA['RVT'][:] #tl.removebounds(DATA['RVT'])
Q0          = Q[0,:,:]
#print(Q[:,0,:].shape,data1D['S_N_direction'])
#plt.contourf(Q[:,0,:]);plt.show()
THT         = DATA['THT'][:] #tl.removebounds(DATA['THT'])
P           = DATA['PABST'][:] #tl.removebounds(DATA['PABST'])
z           = data1D[var1D[0]][:]*1000. # m
TA          = THT*pow(100000./P,-1*CC.gammafix)
TA0         = TA[0,:,:]
idxlcl,lcl  = tl.findlcl(Q0,TA0,P,z)
idxlfc,lfc,idxlnb,lnb = tl.findlfc(idxlcl,Q,TA,z)
ss          = Q.shape
ZZ          = np.repeat(np.repeat(z[ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
offset      = 0.25
THV         = tl.createnew('THV',DATA,var1D)
#THV         = tl.removebounds(THV)
idxzi       = tl.findTHV3D(ZZ,THV,offset)

idxzi2D    = idxzi

if ~np.isnan(idxlcl).all():
  idxlcl      = int(round(np.nanmean(idxlcl)))
else:
  idxlcl      = 62
idxsublcl     = int(round(idxlcl/2.))
if ~np.isnan(idxlfc).all():
  idxlfc      = int(round(np.nanmean(idxlfc)))
if ~np.isnan(idxlnb).all():
  idxlnb      = int(round(np.nanmean(idxlnb)))
#except:
#   pass

#idxlcl,lcl  = tl.createnew('LCL',DATA)
#print type(DATA)
#DATA['LCL'] = lcl
#lcl = DATA.createVariable("LCL","f4")

zview0.update({'sublcl':idxsublcl
        ,'lcl':idxlcl
        ,'lfc':idxlfc
        ,'lnb':idxlnb})
print('zview ',zview0)


# Compute objects once
dataobjs = {}
nameobjs = []
mask     = []
for typ in typs:
  print(typ)
  nameobj   = typ+'_'+thch
  nameobjs += [nameobj] #updr_SVT001_WT_02
  try:
    dataobj   = DATA[nameobj][:] #tl.removebounds(DATA[nameobj])
    # current version
    #dataobjs[nameobj],nbr  = tl.do_delete2(dataobj,tl.do_unique(deepcopy(dataobj)),nbmin,rename=True)
    
    # new version (volume)
    tmpmask = tl.do_unique(deepcopy(dataobj))*nxnynz
    print('tmpmask ',tmpmask[0,:,0])
    dataobjs[nameobj],nbr  = tl.do_delete2(dataobj,tmpmask,vmin,rename=True)
    
    #print nbr,dataobjs[nameobj].shape,np.max(dataobjs[nameobj])
    mask.append(tl.do_unique(deepcopy(dataobjs[nameobj])))
  except:
    dataobjs[nameobj] = None
  #mask.append(tl.do_unique(dataobjs[nameobj]))
  #print mask,dataobjs[nameobj],

if len(mask)>0: # and fluxes:
   nameobj   = 'All'
   nameobjs += [nameobj]
   dataobjs[nameobj] = np.sum(mask,axis=0)

for ij in dataobjs.keys():
  print('NAME : ',ij,np.max(dataobjs[ij]))

objch=''.join(objtyp)
if objch!='':
  objch='_'+objch+'_'+thch+suffixmin


vtypch=''
if vtype!='V0301':
   vtypch='_'+vtype


for vv in vars:
  print('Variable ',vv)
  data      = {} ;
  if vv not in DATA.variables.keys():
     time1  = time.time()
     tmp    = tl.createnew(vv,DATA,var1D)
     time2  = time.time()
     print('Creating new variable %s took %0.3f ms' % ("Variable "+vv, (time2-time1)*1000.0))
     
     # Create new variables for objects
     if objectchar and  ~fluxes and vv=='Frac':
       tmp  = np.ones(DATA['THT'].shape)
  else:
     tmp    = DATA[vv][:]

  data[vv]     = tmp #tl.removebounds(tmp)
  if tmp is not None:
    data[vv] *= tl.findoffset(vv)

    plots     = deepcopy(plots0)
    zview     = deepcopy(zview0)

    if len(data[vv].shape)!=3:
      plots['cross'] = 0; plots['mean'] = 0;
      if objectchar == 0:
        zview = {'lcl':idxlcl}
    
    # clear zview
    zview = {k: zview[k] for k in zview if not np.isnan(zview[k])}

    if fluxes and len(data[vv].shape)==3:
      vv2       = fluxchar+vv
      data[vv2] = DATA[fluxchar]*tl.anomcalc(data[vv])
      #data[vv2] = tl.removebounds(DATA[fluxchar])*tl.anomcalc(data[vv])
      data.pop(vv)
      vv        = vv2

    nametitle = nametitle0.format(var=vv,objch=objch,vtypch=vtypch,suffix=prefix)
    main(data1D,data,plots,dataobjs,xy=xy,zview=zview,pathfig=pathfig,nametitle=nametitle,avg=avg,fluxes=fluxes,idxzi2D=idxzi2D)



