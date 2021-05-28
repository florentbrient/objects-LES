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

      print(dataobjs.keys())
      if plots['cross']:
        mf.plot2D(data1D,data,dataobjs,xy=xy,pathfig=pathfig,nametitle=nametitle,avg=avg,hatchlog=hchcross,idxzi2D=idxzi2D)

      [zview.pop(zz) for zz in zview.keys() if np.isnan(zview[zz])]
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
vars       = ['Reflectance'] #'Reflectance','LWP','WT','THV','RNPM','THLM','DIVUV','WINDSHEAR'] #,'RCT','PRW'] #,'RNPM','RVT','THLM','RCT','THV','DIVUV','REHU','PRW','LWP','LCL','LFC','LNB']

#    With objects?
objectchar = 0

#    Fluxes?
fluxes     = 0 #WT by default
fluxchar   = 'WT'

#    Which plots?
plots0     = {'cross':0, 'zview':1, 'mean':0}

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
filename= 'sel_'+sens+".1."+vtype+"."+prefix
path    = path + case +'/' + sens+'/'
suffixnc= '.nc4'
ncsuff  = ['IHODC','IHOP5',]
if sens in ncsuff:
  suffixnc= '.nc'
file    = path + filename + suffixnc
print(file)
DATA    = nc.Dataset(file,'r')

if objectchar:
  # Objects
  thrs   = 1
  thch   = str(thrs).zfill(2)
  nbmin  = 100 #100 #1000
  objtyp = ['updr','down','down','down']
  objnb  = ['001' ,'001' ,'003' ,'002' ]
  nbplus = tl.svttyp(case,sens) #1
  if nbplus == 1:
     objnb = [str(int(ij)+3).zfill(3) for ij in objnb]
  AddWT  = 1
  WTchar=['' for ij in range(len(objtyp))]
  if AddWT  == 1:
    WTchar = ['_WT','_WT','_WT','_WT']
  typs = [field+'_SVT'+objnb[ij]+WTchar[ij] for ij,field in enumerate(objtyp)]
  print(typs)
  if ~fluxes:
    vars += ['Frac']
else:
  typs = []; objtyp=[]

# Dimensions
var1D  = ['vertical_levels','W_E_direction','S_N_direction']
data1D = OrderedDict()
for ij in var1D:
  data1D[ij] = DATA[ij] #tl.removebounds(DATA[ij][:])
  #print ij,data1D[ij]
#data1D['vertical_levels']/=1000.


# name figure
nametitle0='TTTT_MMMM_{var}{objch}{vtypch}_{suffix}'+'AAAA' # Cross_x1_y2_Mean_Var_updrdowndowndown_2_nb100_008


# Find cloud base, cloud middle and cloud top
subcloud,cloudbase,cloudmiddle,cloudtop,zb,zi=[0 for ij in range(6)]
zview0 = {}
try:
   rct  = DATA['RCT'] #tl.removebounds(DATA['RCT'])
   z    = data1D['vertical_levels']
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
Q           = DATA['RVT'] #tl.removebounds(DATA['RVT'])
Q0          = Q[0,:,:]
print(Q[:,0,:].shape)
plt.contourf(z,data1D['S_N_direction'],Q[:,0,:]);plt.show()
THT         = DATA['THT'] #tl.removebounds(DATA['THT'])
P           = DATA['PABST'] #tl.removebounds(DATA['PABST'])
z           = data1D['vertical_levels']*1000. # m
TA          = THT*pow(100000./P,-1*CC.gammafix)
TA0         = TA[0,:,:]
idxlcl,lcl  = tl.findlcl(Q0,TA0,P,z)
idxlfc,lfc,idxlnb,lnb = tl.findlfc(idxlcl,Q,TA,z)
ss          = Q.shape
ZZ          = np.repeat(np.repeat(z[ :,np.newaxis, np.newaxis],ss[1],axis=1),ss[2],axis=2)
offset      = 0.25
THV         = tl.createnew('THV',DATA)
#THV         = tl.removebounds(THV)
idxzi       = tl.findTHV3D(ZZ,THV,offset)
#plt.contourf(idxzi);plt.show()

#plt.contourf(idxlcl);plt.show()
#plt.contourf(idxlfc);plt.show()
#plt.contourf(idxlnb);plt.show()
#pause
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
    dataobj   = DATA[nameobj] #tl.removebounds(DATA[nameobj])
    dataobjs[nameobj],nbr  = tl.do_delete2(dataobj,tl.do_unique(deepcopy(dataobj)),nbmin,rename=True)
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
  objch='_'+objch+'_'+thch+'_nb'+str(nbmin)


vtypch=''
if vtype!='V0301':
   vtypch='_'+vtype


for vv in vars:
  data      = {} ;
  if vv not in DATA.variables.keys():
     time1  = time.time()
     tmp    = tl.createnew(vv,DATA)
     time2  = time.time()
     print('Creating new variable %s took %0.3f ms' % ("Variable "+vv, (time2-time1)*1000.0))
     
     # Create new variables for objects
     if objectchar and  ~fluxes and vv=='Frac':
       tmp  = np.ones(DATA['THT'].shape)
  else:
     tmp    = DATA[vv]

  data[vv]     = tmp #tl.removebounds(tmp)
  if tmp is not None:
    data[vv] *= tl.findoffset(vv)

    plots     = deepcopy(plots0)
    zview     = deepcopy(zview0)

    print(vv,len(data[vv].shape))
    if len(data[vv].shape)!=3:
      plots['cross'] = 0; plots['mean'] = 0;
      if objectchar == 0:
        zview = {'lcl':idxlcl}

    if fluxes and len(data[vv].shape)==3:
      vv2       = fluxchar+vv
      data[vv2] = tl.removebounds(DATA[fluxchar])*tl.anomcalc(data[vv])
      data.pop(vv)
      vv        = vv2
  
    nametitle = nametitle0.format(var=vv,objch=objch,vtypch=vtypch,suffix=prefix)
    main(data1D,data,plots,dataobjs,xy=xy,zview=zview,pathfig=pathfig,nametitle=nametitle,avg=avg,fluxes=fluxes,idxzi2D=idxzi2D)



