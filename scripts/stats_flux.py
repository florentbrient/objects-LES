# -*- coding: utf-8 -*-
"""
Created on Aug 20
"""
import sys
#sys.path.append('/home/brientf/MESONH/scripts/')
from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
import tools as tl
#from itertools import chain
from copy import deepcopy
#import scipy as sp
#from scipy import ndimage as ndi
import time
from collections import OrderedDict


def format(value):
    return "%13.3f" % value
def format2(value):
    return "%13s" % value

def plottophat(alt,tot,Fi,tophat,covar,subcomp,var,filesave='output',title='No',legend=1): #tot,Fi,tophat,covar,subcomp
  fig = plt.figure()
  size= [6.0,10.0]; lw=2; fts=15
  ax  = fig.add_subplot(111)
  tl.adjust_spines(ax, ['left', 'bottom'])
  plt.plot(tot,alt,'k',label='Tot',lw=lw)
  plt.plot(Fi,alt,'r',label='Fi',lw=lw)
  plt.plot(tophat,alt,'r',linestyle='--',label='Tophat',lw=lw)
  plt.plot(covar,alt,'r',linestyle=':',label='Covar',lw=lw)
  plt.plot(subcomp,alt,'b',linestyle='-.',label='SubComp',lw=lw)
  plt.plot(tophat+subcomp,alt,'k',linestyle='-.',label='Top-hat + SubComp',lw=lw)
  #ax.set_xlim([0.,100.]) 
  #ax.set_ylim([0.,40.])
  xlabel='Fluxes'
  filesave=filesave.replace('XX',var)
  ax.set_xlabel(xlabel,fontsize=fts)
  ax.set_ylabel('Altitude (km)',fontsize=fts)
  ax.set_title(title)
  plt.legend(loc=legend, prop={'size': fts-5})
  plt.xticks(size=fts)
  plt.yticks(size=fts)
  fig.savefig(filesave + '.pdf')
  plt.close()

# call  : python stats_objects.py downdraft SVT006 L25.6 003
typ      = sys.argv[1].split(',') # downdraft
selSVT   = sys.argv[2].split(',') # SVT006
simu     = sys.argv[3] #'L25.6'
hour     = sys.argv[4] #'003'
cas      = sys.argv[5] #'FIRE'
name     = sys.argv[6] #'V0301'
try:
  cond = sys.argv[7]
except:
  cond = ''

# open file
file_dir       = "../data/"+cas+"/"+simu+"/"
filename_input = "sel_"+simu+".1."+name+"."+hour+".nc4"

# out files
dir_out        = "../data/d_stats/"+cas+"/";tl.mkdir(dir_out)
dir_out       += simu+"/";tl.mkdir(dir_out)

liste_fichiers = file_dir+filename_input
print(liste_fichiers)

# Open data
DATA = Dataset(liste_fichiers)

#################
# Dimensions (similar to plot_LES_3D.py)
var1D  = ['vertical_levels','W_E_direction','S_N_direction'] #Z,Y,X
data1D,nzyx,sizezyx = [OrderedDict() for ij in range(3)]
for ij in var1D:
  data1D[ij]  = DATA[ij][:]/1000. #km #tl.removebounds(DATA[ij][:])
  nzyx[ij]    = data1D[ij][1]-data1D[ij][0]
  sizezyx[ij] = len(data1D[ij])

nxny   = nzyx[var1D[1]]*nzyx[var1D[2]] #km^2
ALT    = data1D[var1D[0]]
dz     = [0.5*(ALT[ij+1]-ALT[ij-1]) for ij in range(1,len(ALT)-1)]
dz.insert(0,ALT[1]-ALT[0])
dz.insert(-1,ALT[-1]-ALT[-2])
nxnynz = np.array([nxny*ij for ij in dz]) # volume of each level
nxnynz = np.repeat(np.repeat(nxnynz[:, np.newaxis, np.newaxis]
                             , sizezyx[var1D[1]], axis=1), sizezyx[var1D[2]], axis=2)

#print('1: ',nxnynz.shape)
#print('2: ',nxnynz[:,0,0])
#print('3: ',nxnynz[0,0,:])
  
X    = DATA.variables[var1D[1]][:]
Y    = DATA.variables[var1D[2]][:]
nx,ny,nz = (X[1]-X[0])/1e3,(Y[1]-Y[0])/1e3,(ALT[1]-ALT[0])/1e3
#################


WT   = DATA.variables['WT'][:]
THLM = DATA.variables['THLM'][:]
RNPM = DATA.variables['RNPM'][:]
print(WT.shape)

# zmax = Altitude maximum to average fluxes !!
# Very important
zmax     = WT.shape[0]
if cas == 'FIRE':
  altmax = 620. 
  #zmax   = 62 # maximum layer
  zmax   = tl.near(ALT,altmax)

# Resize
WT       = WT[0:zmax,:,:]
THLM     = THLM[0:zmax,:,:]
RNPM     = RNPM[0:zmax,:,:]*1000.
dz       = dz[0:zmax]
nxnynz   = nxnynz[0:zmax,:,:]

# Calculate anomalies
ATHLM    = tl.anomcalc(THLM)
ARNPM    = tl.anomcalc(RNPM)
# Domain-mean fluxes
WTTHLM   = WT*ATHLM;print(np.mean(WTTHLM))
WTRNPM   = WT*ARNPM;print(np.mean(WTRNPM))
# Domain-mean average
WTM      = np.nanmean(WT,axis=(1,2))
THM      = np.nanmean(THLM,axis=(1,2))
RTM      = np.nanmean(RNPM,axis=(1,2))
# Total volume
vtot     = np.sum(nxnynz)

##### Object selection
# m value (strength of the conditional sampling)
listthrs = [0,1,2,3,4] #[0]
N1       = len(listthrs)
# Vmin (can be defined in units or volume in km^3)
minchar = 'volume' #'volume' #unit
if minchar == 'volume':
    nbmins   = [0.005,0.01,0.02,0.05,0.10,0.25]
else:
    nbmins   = [0,10,50,100,500,1000,5000,10000,20000]
N2 = len(nbmins) #,50,100,500,1000,5000,10000,20000];
#####

nbr,surf,volume,alpha,fluxTHLM,fluxRNPM,fremTHLM,fremRNPM,fluxTHLM2,fluxRNPM2,fremTHLM2,fremRNPM2,altmin,altmax       = [np.zeros((N1,N2)) for ij in range(14)]

# Test with one value
maketest = False
testch   = ''
if maketest:
  listthrs = [2] #[2]
  #nbmins   = [nbmins[4]] #[0,10,50,100,500,1000,5000,10000,20000] #[10000]
  testch   = '.test'

### Name the output file
# if condition AND/OR
jointyp  = cond
# First letter
typ_first = [x[0] for x in typ]
# info tracer
sel_first = [x[-4:] for x in selSVT]
# Important info from filename
print(filename_input.split('.'))
file_short='_'.join(np.array(filename_input.split('.'))[[0,-2]])
file_short=file_short.replace('sel_', '') 
# Name of the output file
fileout0  = dir_out
#fileout0 += 'stats.'+jointyp.join(typ)+'.'+''.join(selSVT)+'.LL.'+filename_input+'.3'+testch+'.txt'
fileout0 += 'stats.'+''.join(typ_first)+'.'+''.join(sel_first)+'.LL.'+file_short+testch+'.txt'

# Run along the number of m values (conditional sampling)
#id1=-1
for id1,ll in enumerate(listthrs):
  OBJ      = []
  nameobjs = []
  #id1     +=1
  for ik in range(len(typ)):
    nameobjs.append(typ[ik][0:4]+'_'+selSVT[ik]+'_'+str(ll).zfill(2))
    print('nameobjs ',ll,nameobjs[ik])
    OBJ.append(DATA.variables[nameobjs[ik]][:][0:zmax,:,:])

  sz        = OBJ[0].shape
  frac,fTHLM,fTHLMfr,fTHLMrem,fTHLMtot,fRNPM,fRNPMfr,fRNPMrem,fRNPMtot= [np.zeros(sz[0]) for ij in range(9)]
  WTI, THI, RTI = [np.zeros(sz[0]) for ij in range(3)]
  tophatTHLM, intraTHLM, tophatRNPM, intraRNPM = [np.zeros(sz[0]) for ij in range(4)]

  fileout = fileout0.replace('LL',str(ll)); 
  f = open(fileout, 'w')
  f.write("File name : "+liste_fichiers+"\n")
  f.write("Object types : "+''.join(selSVT)+"\n")
  f.write(jointyp.join(typ)+" threshold : "+str(ll*0.5)+"\n")
  names  = ['nbmin',' nb',' surf',' vol',' rvol',' rTHLMflux',' rRNPMflux',' altmin',' altmax\n']
  unity  = ['-','-','km2','km3','%','%','%','km','km\n'] 
  for v in names:
    f.write(format2(v))
  for v in unity:
    f.write(format2(v))
  id2     =-1
  for ij,bin in enumerate(nbmins):
    print('ij,bin : ',ij,bin)
    id2    +=1
    mask   = tl.do_unique(deepcopy(OBJ[0]))
    if minchar == 'volume':
        mask = mask*nxnynz # Not good... 2D not 3D
    time1 = time.time()
    OBJ[0],nbr[id1,id2]  = tl.do_delete2(OBJ[0],mask,bin,rename=True)
    print(nbr[id1,id2])
    time2 = time.time()
    print('%s function 1 took %0.3f s for Vmin=%7.3f' % ("delete", (time2-time1)*1.0, bin))
    # Check volume or number
    #labs = np.unique(OBJ[0])
    #print('labs ',labs)
    #for lab in labs:
    #  print(lab,np.sum(mask[OBJ[0]==lab]))
        
    if nbr[id1,id2] >-1:
     condition         = (OBJ[0]!=0)
     for ik in range(len(typ)-1):
       ikk              = ik+1
       mask             = tl.do_unique(deepcopy(OBJ[ikk]))
       if minchar == 'volume':
         mask = mask*nxnynz # Not good... 2D not 3D
       OBJ[ikk],nbrtmp  = tl.do_delete2(OBJ[ikk],mask,bin,rename=True)
       nbr[id1,id2]    += nbrtmp
       if jointyp  == 'and':
         condition       = (condition & (OBJ[ikk]!=0))
       elif jointyp  == 'or':
         condition       = (condition | (OBJ[ikk]!=0))

     surf[id1,id2]     = np.sum(condition)*nxny #km2
     volume[id1,id2]   = np.sum(nxnynz*condition)#*nxnynz #km3
     alpha[id1,id2]    = volume[id1,id2]/vtot
     
     # ratio of fluxes (first try)
     #fluxTHLM[id1,id2] = 100.*(alpha[id1,id2])*np.nanmean(WTTHLM[condition])/np.mean(WTTHLM)
     #fluxRNPM[id1,id2] = 100.*(alpha[id1,id2])*np.nanmean(WTRNPM[condition])/np.mean(WTRNPM)
     # ratio of fluxes (second try -> good)
     for ik in range(sz[0]):
       cd       = condition[ik,:,:]
       tmp0     = WTTHLM[ik,:,:]; tmp1 = WTRNPM[ik,:,:]
       WT0      = WT[ik,:,:]
       TH0      = THLM[ik,:,:];    RT0 = RNPM[ik,:,:]

       frac[ik] = np.sum(cd)/float(sz[1]*sz[2])
       #print ik,cd.shape,np.sum(cd),float(sz[1]*sz[2]),frac[ik]

       # CS flux contribution
       #if ik==30:
       #  print tmp0[:]
       #  print tmp0[cd]
       #  print len(tmp0.ravel()),len(tmp0[cd]),frac[ik]
       #  print np.nansum(tmp0[cd])/len(tmp0[cd])
       #  print np.nanmean(tmp0[cd])
       #  #pause
       fTHLMfr[ik]  = frac[ik]*np.nanmean(tmp0[cd])
       fTHLM[ik]    = np.nanmean(tmp0[cd])
       # remaining part
       fTHLMrem[ik] = (1.0-frac[ik])*np.nanmean(tmp0[~cd])
       fTHLMtot[ik] = np.mean(tmp0)
       # top-hat approx
       WTI[ik]        = np.nanmean(WT0[cd])
       THI[ik]        = np.nanmean(TH0[cd])
       tophatTHLM[ik] = (WTI[ik]-WTM[ik])*(THI[ik]-THM[ik])
       # internal variability
       intraTHLM[ik]  = np.nanmean((WT0[cd]-WTI[ik])*(TH0[cd]-THI[ik]))
       #print 'z, decomp, top-hat, intra ',ik, fTHLMfr[ik], tophatTHLM[ik], intraTHLM[ik]
       
       

       # CS flux contribution
       fRNPMfr[ik]  = frac[ik]*np.nanmean(tmp1[cd])
       fRNPM[ik]    = np.nanmean(tmp1[cd])
       # remaining part
       fRNPMrem[ik] = (1.0-frac[ik])*np.nanmean(tmp1[~cd])
       fRNPMtot[ik] = np.mean(tmp1)
       # top-hat approx
       RTI[ik]        = np.nanmean(RT0[cd])
       tophatRNPM[ik] = (WTI[ik]-WTM[ik])*(RTI[ik]-RTM[ik])
       # internal variability
       intraRNPM[ik]  = np.nanmean((RT0[cd]-RTI[ik])*(WT0[cd]-WTI[ik]))
       #print 'z, decomp, top-hat, intra ',ik, fTHLMfr[ik], tophatTHLM[ik], intraTHLM[ik]

     #plt.plot(WTI);plt.show()
     #plt.plot(WTM);plt.show()
     #plt.plot(THI);plt.show()
     alts = ALT[:zmax]
     #plt.plot(frac,alts);plt.show()

     # Flux totaux, Fi, tophat, intra-inter, sub compensatoire
     #filesave='decomposition_flux_'+var+'_'+typ
     filesave  = dir_out+'decomposition_flux_'+'XX'+'_'+jointyp.join(typ)+'_'+''.join(selSVT)
     filesave += '_'+simu+'_'+hour+'_'+str(ll)+'_'+str(bin)
     if name != 'V0301':
        filesave += '_'+name
     title     = 'Flux '+cas+' '+simu+' '+hour \
               +' '+jointyp.join(typ)+' (m='+str(float(ll)/2.)+',Vmin='+str(bin)+')'

     alphai = frac
     tot=fTHLMtot; Fi=alphai*fTHLM; tophat=alphai*tophatTHLM; 
     covar=alphai*intraTHLM; subcomp=(alphai*alphai/(1.0-alphai))*tophatTHLM
     plottophat(alts,tot,Fi,tophat,covar,subcomp,'THLM',filesave,title,legend=3)

     tot=fRNPMtot; Fi=alphai*fRNPM; tophat=alphai*tophatRNPM; 
     covar=alphai*intraRNPM; subcomp=(alphai*alphai/(1.0-alphai))*tophatRNPM
     plottophat(alts,tot,Fi,tophat,covar,subcomp,'RNPM',filesave,title,legend=1) #4


     # Vertical averaging
     fluxTHLM[id1,id2] = 100.*np.nansum(abs(fTHLMfr*dz))/np.nansum(abs(fTHLMtot*dz))
     fluxRNPM[id1,id2] = 100.*np.nansum(abs(fRNPMfr*dz))/np.nansum(abs(fRNPMtot*dz))
     print(id1,id2,np.nansum(abs(fTHLMfr)),np.nansum(abs(fTHLMtot)),np.nansum(abs(fTHLMfr*dz))/np.nansum(abs(fTHLMtot*dz)))
     print(id1,id2,np.nansum(abs(fRNPMfr)),np.nansum(abs(fRNPMtot)),np.nansum(abs(fRNPMfr*dz))/np.nansum(abs(fRNPMtot*dz)))

     fremTHLM[id1,id2] = 100.*np.nansum(abs(fTHLMrem*dz))/np.nansum(abs(fTHLMtot*dz))
     fremRNPM[id1,id2] = 100.*np.nansum(abs(fRNPMrem*dz))/np.nansum(abs(fRNPMtot*dz))


#     fluxTHLM2[id1,id2] = 100.*np.nansum(abs(fTHLMfr*dz))/np.nansum(abs(fTHLMfr*dz)+abs(fTHLMrem*dz))
#     fluxRNPM2[id1,id2] = 100.*np.nansum(abs(fRNPMfr*dz))/np.nansum(abs(fRNPMfr*dz)+abs(fRNPMrem*dz))
#     print id1,id2,np.nansum(abs(fTHLMfr)),np.nansum(abs(fTHLMfr)+abs(fTHLMrem)),np.nansum(abs(fTHLMfr*dz))/np.nansum(abs(fTHLMfr*dz)+abs(fTHLMrem*dz))
#     print id1,id2,np.nansum(abs(fRNPMfr)),np.nansum(abs(fRNPMfr)+abs(fRNPMrem)),np.nansum(abs(fRNPMfr*dz))/np.nansum(abs(fRNPMfr*dz)+abs(fRNPMrem*dz))
#     fremTHLM2[id1,id2] = 100.*np.nansum(abs(fTHLMrem*dz))/np.nansum(abs(fTHLMfr*dz)+abs(fTHLMrem*dz))
#     fremRNPM2[id1,id2] = 100.*np.nansum(abs(fRNPMrem*dz))/np.nansum(abs(fRNPMfr*dz)+abs(fRNPMrem*dz))

     #print fluxTHLM[id1,id2],fluxRNPM[id1,id2],alpha[id1,id2]
     
     # size
     if len(typ)==1:
       listObjs          = np.unique(OBJ[0][OBJ[0]!=0]);#print listObjs
       min0,max0         = [],[]
       for il in listObjs: #range(max(listObjs)):
          cd = np.where(OBJ[0].reshape((sz[0],sz[1]*sz[2]))==il)
          # Error? ? si nz change with height???
          min0.append(np.min(cd[0])*nz*1e3); max0.append(np.max(cd[0])*nz*1e3) 
          #print il,cd,altmin,altmax
       altmin[id1,id2]   = np.mean(min0)
       altmax[id1,id2]   = np.mean(max0)

    tab               = [bin,nbr[id1,id2],surf[id1,id2],volume[id1,id2],100.*alpha[id1,id2],fluxTHLM[id1,id2],fluxRNPM[id1,id2],altmin[id1,id2],altmax[id1,id2]]
    print(tab)
    for v in tab:
      f.write(str(format(v)))
    f.write("\n")
  f.close()




