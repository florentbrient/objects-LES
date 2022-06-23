# -*- coding: utf-8 -*-

import sys
sys.path.append('/cnrm/tropics/user/brientf/MESONH/scripts/')
from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
import tools2 as tl
import Constants as CC
from itertools import chain
from copy import deepcopy
import scipy as sp
from scipy import ndimage as ndi
import time
from skimage import measure
import netCDF4 as ncdf
from collections import OrderedDict
import makefigures as mf
# open infocas
from infoplot_home import *

def format(value):
    return "%9.3f" % value

def findvars1D(var):
    name1D = ['Bowen','H0','Wstar','Tstar','Ra']
    if var in name1D:
      #if var == 'E0':
      #  var  = ['E0','MEAN_TH','MEAN_RV','MEAN_PRE']
      #if var == 'Q0':
      #  var  = ['Q0','MEAN_TH','MEAN_RV','MEAN_PRE']
      if var == 'Bowen':
        var  = ['E0','Q0']
      if var == 'H0':
        var  = ['E0','Q0','MEAN_TH','MEAN_RV','MEAN_PRE']
      if var == 'Wstar' or var == 'Tstar' or var == 'Ra':
        var  = ['E0','Q0','MEAN_TH','MEAN_RV','MEAN_THV','MEAN_PRE']
    else:
      var    = [var]
    return var

def findsvt2(cas,subcloud=False):
#svt = findsvt2(scat,cas,subcloud=subcloud)
    svt       = 'SVT003'
    namenocld = {'IHOP','IHODC','AYOT03S','AYOT05S','AYOT05W','AYOT24S'\
                ,'FIREnoc','FIREopen','FIRE','BOMEXnoc','ARMCUnoc'}
    if (subcloud and cas not in namenocld):
       svt       = svt.replace('003','002')
    return svt

def findsvt(name):
    svtchar = {}
    svtchar['rel_down_hum']={'styp':'SVT003_WT','sname':'downdraft','index':'rRNPMflux'}
    svtchar['rel_down_th'] ={'styp':'SVT003_WT','sname':'downdraft','index':'rTHLMflux'}
    svtchar['rel_up_hum']  ={'styp':'SVT001_WT','sname':'updraft','index':'rRNPMflux'}
    svtchar['rel_up_th']   ={'styp':'SVT001_WT','sname':'updraft','index':'rTHLMflux'}
    svtchar['pos_down_hum']={'styp':'SVT003_WT','sname':'downdraft','index':'pRNPMfl'}
    svtchar['neg_down_hum']={'styp':'SVT003_WT','sname':'downdraft','index':'nRNPMfl'}
    svtchar['pos_down_th'] ={'styp':'SVT003_WT','sname':'downdraft','index':'pTHLMfl'}
    svtchar['neg_down_th'] ={'styp':'SVT003_WT','sname':'downdraft','index':'nTHLMfl'}
    svtchar['pos_up_hum']  ={'styp':'SVT001_WT','sname':'updraft','index':'pRNPMfl'}
    svtchar['neg_up_hum']  ={'styp':'SVT001_WT','sname':'updraft','index':'nRNPMfl'}
    svtchar['pos_up_th']   ={'styp':'SVT001_WT','sname':'updraft','index':'pTHLMfl'}
    svtchar['neg_up_th']   ={'styp':'SVT001_WT','sname':'updraft','index':'nTHLMfl'}
    svtchar['WTm_up']      ={'styp':'SVT001_WT','sname':'updraft','index':'WTmax'}
    svtchar['WTm_down']    ={'styp':'SVT003_WT','sname':'downdraft','index':'WTmax'}
    svtchar['WTm_down_1']  ={'styp':'SVT001_WT','sname':'downdraft','index':'WTmax'}
    svtchar['vol_up']      ={'styp':'SVT001_WT','sname':'updraft','index':'vol'}
    svtchar['vol_down']    ={'styp':'SVT003_WT','sname':'downdraft','index':'vol'}
    svtchar['MF_down']     ={'styp':'SVT003_WT','sname':'downdraft','index':'MassFlux'}
    svtchar['MF_up']       ={'styp':'SVT001_WT','sname':'updraft','index':'MassFlux'}
    return svtchar[name]


def updatename(scat,name,subcloud=False):
    #upname = {}
    #upnames= ['rel_down_hum','rel_down_th','rel_up_hum','rel_up_th',\
    #          'WTm_down','WTm_up',\
    #          'pos_down_hum','neg_down_hum','pos_down_th','neg_down_th']
    #if name in upnames: #for key in upnames:
    #  upname={'ARMCU','RICO','BOMEX'}
    namenocld={'IHOP','IHODC','AYOT24S','AYOT03S','AYOT05S','AYOT05W','FIREnoc',\
               'FIRE','FIREopen','BOMEXnoc','ARMCUnoc'}
    #upname['rel_down_hum']={'ARMCU','RICO','BOMEX'}#,'ASTEX','FIRE'}
    #upname['rel_down_th'] ={'ARMCU','RICO','BOMEX'}#,'ASTEX','FIRE'}
    #upname['WTm_down']    ={'ARMCU','RICO','BOMEX'}#,'ASTEX','FIRE'}
    #print cas,upname,name
    subprefix = ''
    svt       = findsvt(scat)['styp']
#    if cas in upname:
#      if 'down' in name:
#        svt = svt.replace('003','002');subprefix='.scld'
    print subcloud,name,namenocld,(not subcloud),scat
    if (subcloud and name not in namenocld): #or \
       #(not subcloud  and 'down' in scat):
       svt       = svt.replace('003','002');
       subprefix = '.scld'
      #if 'up' in name:
      #  subprefix='.scld'
    return svt,subprefix

def opennetcdf(path,info,name,hour,scat,var1D,subcloud):
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
    Data,filenc,cas,simu  = tl.opennc(path,info,name,hour,nc=nc,vtype=vtype,rmsimu=True)
    # Open 1D netcdf 
    in1 = 0; in2 = 0; Data1D=[]
    if var1D:
     while in1 == 0:
      try:
        in2 += 1;
        vtype2   = vtype
        if vtype=='V0301':
          vtype2 = 'V03'+str(in2).zfill(2)
        if not(vtype2[0:2]=='V0'):
          in1=1
        print path,info,name,vtype2
        tmp,tmp1,tmp2,tmp3 = tl.opennc(path,info,name,'000',nc=nc,vtype=vtype2,rmsimu=True)
        #print in1,vtype2,tmp
        Data1D=np.append(Data1D,tmp)
      except:
        in1 = 1
        pass
    print scat,var1D

    if not var1D:
#      svt,subprefix = updatename(scat,cas,subcloud=subcloud)
      svt = findsvt2(name,subcloud=subcloud)
      if 'svt' in info.keys():
         svt  = tl.replacename(svt)
    else:
      svt = None

    # Create index
    print 'svt ',svt
    data = createindex(Data,Data1D,scat,svt=svt,hour=hour,var1D=var1D,name1=name1)
    return data

def opentxt(fileout0,info,name,hour,scat,mm,vtype,subcloud):
    typ     = findsvt(scat)['sname']
    simu    = info['sens']
    cas     = name
    if 'cas' in info.keys():
       cas  = info['cas']
    print name
    svt,subprefix = updatename(scat,name,subcloud=subcloud)   
    if 'svt' in info.keys():
       svt  = tl.replacename(svt)   
    if 'vtyp' in info.keys():
       vtype= info['vtyp']

    fileout = fileout0.format(cas=cas,typ=typ,svt=svt,mm=mm,simu=simu,vtype=vtype).replace('SSSS',subprefix)
    fileout = fileout.replace('TTTT',hour+'.nc4') 
    if 'nc' in info.keys():
       fileout=fileout.replace('nc4',info['nc'])
    print fileout
    tab,nametab,Nmin,Nval,ijmin  = tl.openfilestat(fileout)
    idx     = nametab.index(findsvt(scat)['index'])
    data    = np.zeros((Nmin,Nval))
    for ij in range(Nmin):
      for ik in range(Nval):
        data[ij,ik]=float(tab[ijmin+ij][ik])
    idxvmin   = tl.near(data[:,0],info['vmin'])
    #print fileout,idxvmin,data[idxvmin,0]
    tmp     = data[idxvmin,idx]
    return tmp

def find_diff(ZZ,data,idx=None):
    #dz     = 3
    idz    = 0.1 #km
    dz     = ZZ[1]-ZZ[0]
    idzp   = int(round(idz/dz)) #4 or more 
    idxall = np.arange(idx-idzp,idx+idzp+1)
    idx0   = np.arange(0,2)
    #datax  = (np.mean(data[idxall])-np.mean(data[idx0])) #\
             #/(np.mean(ZZ[idxall])-np.mean(ZZ[idx0]))
    idx1   = np.arange(idx+1,idx+idzp+1)
    idx2   = np.arange(idx-idzp-1,idx-1)
    datax  = np.mean(data[idx1])-np.mean(data[idx2])
    print idx,ZZ[idx],datax
    #plt.plot(data,ZZ);plt.show()
    return datax

def gradient(data,ZZ,nn):
    #print int(nn/2)
    #plt.plot(data,ZZ);plt.show()
    data  = tl.smooth(data,nn)
    dTHdz = np.gradient(data)/np.gradient(ZZ)
    dTHdz[:int(nn/2)+1]=np.nan
    dTHdz[len(dTHdz)-int(nn/2):]=np.nan
    idx   = np.argmax((dTHdz>0))
    #plt.plot(dTHdz,ZZ);plt.plot(dTHdz[idx],ZZ[idx],'ro');plt.show()
    return dTHdz,idx

def THVtop(data,ZZ,offset):
   THVint = sp.integrate.cumtrapz(data)/np.arange(1,len(data))
   offset = 0.25
   DT     = data[:-1]-(THVint+offset)
   idx    = np.argmax(DT>0)
   return idx

def create1D(Data1D,ij,xaxis,vars1D,name1='--',idx=None):
    datax = dict()
    for vv in vars1D:
      name  = vv+name1+'PROC1'
      print vv,name,np.squeeze(Data1D[ij][name][:]).shape
      datax[vv]=np.squeeze(Data1D[ij][name][:])*tl.findoffset(vv) #CC.RLVTT
      
    if xaxis=='Bowen': #
      datax[xaxis] = datax['Q0']/datax['E0']
    if xaxis=='E0' or xaxis=='Q0': #
      THT   = datax['MEAN_TH'][:,0]
      P     = datax['MEAN_PRE'][:,0]
      T     = tl.tht2temp(THT,P)
      RHO   = tl.createrho(T,P)
      datax[xaxis]*= RHO
    if xaxis=='H0' or xaxis=='Ra' or 'star' in xaxis: #
      THT   = datax['MEAN_TH'][:,0]
      P     = datax['MEAN_PRE'][:,0]
      T     = tl.tht2temp(THT,P)
      RHO   = tl.createrho(T,P)

      q0  = datax['MEAN_RV'][:,0]
      t0  = datax['MEAN_TH'][:,0]
      datax[xaxis] =(1.0+0.61*q0)*datax['Q0']\
                   +0.61*t0*datax['E0']*CC.RCP/CC.RLVTT
      datax[xaxis]*=RHO
      #print xaxis,datax[xaxis],datax['E0'],datax['Q0']
      #plt.plot(datax[xaxis]);plt.show()

      if xaxis=='Wstar' or xaxis=='Tstar' or xaxis=='Ra': #
        # g*zi*H0/theta_v(0-zi)
        THV = datax['MEAN_THV']
        ZZ  = tl.removebounds(Data1D[ij]['ZHAT'][:]) #datax['ZHAT'][:]
        ss  = THV.shape
        #ZZ  = np.repeat(ZZ[ :,np.newaxis],ss[1],axis=1)
        nn  = 5; offset = 0.25
        if xaxis=='Tstar'or xaxis=='Ra':
          H0 = deepcopy(datax[xaxis])
        for ij in range(ss[0]):
          #dTHdz,idx = gradient(THV[ij,:],ZZ,nn)
          #if idx is None:
          idx       = THVtop(THV[ij,:],ZZ,offset)
          #print idx,datax[xaxis][ij]
          #if ij == int(ss[0]/2.):
          #  plt.plot(THV[ij,:],ZZ);plt.plot(THV[ij,idx],ZZ[idx],'ro');plt.show()
          #print idx,ZZ[idx],datax[xaxis][ij],np.nanmean(THV[ij,0:idx])
          datax[xaxis][ij]=CC.RG*ZZ[idx]*datax[xaxis][ij]/np.nanmean(THV[ij,0:idx])
          datax[xaxis][ij]=pow(datax[xaxis][ij]/CC.RCP,1./3.)
          if xaxis=='Ra': #Ra = Beta*G*Qs*zi^3/(tau*nu)/Wstar
            tmp              = (CC.RG*pow(ZZ[idx],3.)*(H0[ij]/CC.RCP))/np.nanmean(THV[ij,0:idx])
            datax[xaxis][ij] = tmp/datax[xaxis][ij]
            nu               = 15e-6 #m2/s
            mu               = 21.4e-6 #m2/s Viscosity
            datax[xaxis][ij]/= nu*mu # Replace by rho (density?)
          
        if xaxis=='Tstar': #Qs/Wstar
          print xaxis,H0[-1],datax[xaxis][-1]
          datax[xaxis]    = (H0/CC.RCP)/datax[xaxis]
    return datax[xaxis]

def createindex(Data,Data1D,xaxis,svt=None,hour=None,var1D=False,name1='--'):
    # Index for PBL height (top or sub-cloud)
    newidx=None
    if svt is not None:
      newidx = tl.findmax(np.squeeze(Data[svt[:6]]))
    print newidx

    ZZ     = Data['ZHAT'][:]/1000.
    THT    = np.mean(np.squeeze(Data['THT']),axis=(1,2,))
    P      = np.mean(np.squeeze(Data['PABST']),axis=(1,2,))
    #RVT    = np.mean(np.squeeze(Data['RVT']),axis=(1,2,))
    T      = tl.tht2temp(THT,P)
    #RHO    = createrho(T,P)

    print 'Xaxis is : ',xaxis
    
    if xaxis == 'THLMdiff':
      THLM   = np.mean(np.squeeze(Data['THLM']),axis=(1,2,))
      datax = find_diff(ZZ,THLM,idx=newidx)
    if xaxis == 'THS1diff':
      THS1  = np.mean(tl.createnew('THS1',Data),axis=(1,2,))
      datax = find_diff(ZZ,THS1,idx=newidx)
    if xaxis == 'QTdiff':
      RNPM   = np.mean(np.squeeze(Data['RNPM']),axis=(1,2,))
      datax = find_diff(ZZ,RNPM,idx=newidx)
    if xaxis == 'THTdiff':
      datax = find_diff(ZZ,THT,idx=newidx)
    if xaxis == 'MSEdiff':
      datax = find_diff(ZZ,Data,idx=newidx)
    if xaxis == 'THVdiff' or xaxis == 'we' or xaxis == 'DZe':
      data  = tl.createnew('THV',Data)
      THV   = np.mean(np.squeeze(data),axis=(1,2,))
      print newidx,ZZ[newidx],THV.shape
      datax = find_diff(ZZ,THV,idx=newidx)
      if xaxis == 'we':
        WT    = np.squeeze(Data['WT'])
        data2 = np.nanmean(tl.anomcalc(WT)[newidx,:,:]*tl.anomcalc(data)[newidx,:,:])
        datax = data2/datax
      if xaxis == 'DZe':
        WT    = np.squeeze(Data['WT'])
        data2 = np.nanmean(tl.anomcalc(WT)*tl.anomcalc(data),axis=(1,2))
        idx0  = np.where(data2<0)[0][0]
        print 'idx0  ',idx0,ZZ[idx0]
        idx1  = np.where(data2[idx0:]>0)[0][0]+idx0
        print 'idx1  ',idx1,ZZ[idx1]
        datax = ZZ[idx1]-ZZ[idx0]
      # essai
      #nn    = 5
      #print int(nn/2)
      #THV   = tl.smooth(THV,nn)
      #dTHdz = np.gradient(THV)/np.gradient(ZZ)
      #dTHdz[:int(nn/2)-1]=np.nan
      #dTHdz[len(dTHdz)-int(nn/2):]=np.nan
      #idx   = np.argmax((dTHdz>0))
      #plt.plot(dTHdz,ZZ);plt.plot(dTHdz[idx],ZZ[idx],'ro');plt.plot(dTHdz[newidx],ZZ[newidx],'bo');plt.show()
    if xaxis == 'maxSV': #
      #print newidx,svt[:6],np.nanmean(np.squeeze(Data[svt[:6]]),axis=(1,2,))
      datax = ZZ[newidx]
    if xaxis == 'Kappa': #
      print ZZ,THLM,newidx
      THLM   = np.mean(np.squeeze(Data['THLM']),axis=(1,2,))
      dataT = find_diff(ZZ,THLM,idx=newidx)
      RNPM   = np.mean(np.squeeze(Data['RNPM']),axis=(1,2,))
      dataQ = find_diff(ZZ,RNPM,idx=newidx)
      epsilon = 10**-6
      if abs(dataQ)>epsilon:
        datax = 1.0+dataT/((CC.RLVTT/CC.RCP)*(dataQ))
      else:
        datax = np.nan
      print datax
    if xaxis == 'NNmax' or xaxis == 'NNmin':
      data = tl.createnew('NN',Data)
      if data is not None:
        data = np.mean(np.squeeze(data),axis=(1,2,))
      if xaxis == 'NNmax':
        maxposition = np.argmax(data)
      if xaxis == 'NNmin':
        maxposition = np.argmin(data)
      print 'max position ',maxposition,ZZ[maxposition],data[maxposition]
      datax = data[maxposition]
    if 'skew' in xaxis:
      typ   = xaxis.split('_')[-1]
      if newidx is None:
        idx =  Data[typ].shape[0]
      else:
        idx = newidx
      data  = tl.anomcalc(Data[typ][:idx,:,:]).flatten()
      datax = sp.stats.skew(data)
      

    if var1D:#xaxis == 'fluxLH':
      print len(Data1D)
      N1D   = len(Data1D)
      Y1D   = True
      datax = np.nan
      t0    = 0
      vars1D=findvars1D(xaxis)
      name  = vars1D[0]+name1+'PROC1'
      print vars1D,xaxis
      for ij in range(N1D):
        if Y1D:
          nametime = name.replace('PROC1','DATIM')
          time  = Data1D[ij][nametime][:][:,-1] # hour
          DT    = (time[1]-time[0])
          print ij,time[-1]/60./60.
          time  = t0+(np.arange(DT,(DT)*(time.shape[0]+1),DT))/60./60. #hours
          idxtime = np.where(time==float(hour))[0]
          print ij,N1D,time[-1],hour,idxtime
          if idxtime:
            datax = create1D(Data1D,ij,xaxis,vars1D,name1=name1,idx=newidx)
            datax = datax[idxtime]
            print idxtime,time[idxtime],datax
            Y1D   = False
          else:
            t0    = time[-1]
      #datax =  np.nanmean((Data['WT'][1,:,:])*(Data['RNPM'][1,:,:]-RNPM[1]))
    if xaxis == 'gradRH' or xaxis == 'gradQT' or xaxis == 'gradTHT'  or xaxis == 'gradTHLM':
      idxpbl,toppbl,grad = tl.findpbltop(xaxis,Data,idx=newidx)
      datax = toppbl
    #if xaxis == 'gradTHT':
    #  idxpbl,toppbl,grad = tl.findpbltop(xaxis,Data,idx=newidx)
    #  datax = grad[newidx] #toppbl
    if xaxis == 'maxTHV':
      idxpbl,toppbl,grad = tl.findpbltop(xaxis,Data,idx=newidx)
      datax = toppbl
    if xaxis == 'PBLhRH':
      idxpbl,toppbl,datax= tl.findpbltop(xaxis,Data)
    return datax

def findmax(name):
    out    = (None,None)
    maxmin = dict()
    maxmin = {'WTm_up':(0.,2.7),'WTm_down':(-2.7,0.),'MF_up':(0.,0.3),'MF_down':(-0.1,0.),\
              'vol_up':(0.,100.),'vol_down':(0.,50.),}
    if name in maxmin.keys():
      out = maxmin[name]
    return out


def main(path,infocas,scatter,vtype='V0301',mm='2',testch='',name1D=[],vartxt=[],subcloud=False,colorplot=None):
    if __name__== "__main__" :
      pathout   = path+'d_unif/'
      tl.mkdir(pathout)
      pathncdf  =  '/Volumes/TOSHIBA/'
      #file0     = pathncdf+"{cas}/{simu}.1.{vtype}.TTTT.nc4"
      # ex: stats/d_FIRE/stats.downdraft.SVT006.2.Ls2x0.1.V0301.024.nc4.3.txt
      fileout0  = pathncdf+'/stats/d_{cas}/'
      fileout0 += 'stats.{typ}.{svt}.{mm}.{simu}.1.{vtype}.TTTT.3'+testch+'SSSS'+'.txt'

      fig = plt.figure()
      ax  = fig.add_subplot(111)
      cm  = plt.cm.get_cmap('Oranges')
      print scatter
      names     = infocas.keys()
      #datax     = np.zeros(len(names))
      #datay     = np.zeros(len(names))
      #print names
       
      xx,yy,colorall  = [],[],[]
      for nn,name in enumerate(names):
        info  = infocas[name]
        print name
        cas   = name
        if 'cas' in info.keys():
         cas  = info['cas']
        for hh,hour in enumerate(info['hour']):
          data = []
          for ij,scat in enumerate(scatter):
            var1D   = False
            if scat[0] in name1D:
              var1D = True
            if cas == 'AYOTTE':
              scat = [s.replace('THLM','THT') for s in scat]
              scat = [s.replace('THV','THT') for s in scat]

            print scat,vartxt
            if scat[0] in vartxt:
              tmp = []
              for stmp in scat: 
                tmp  = np.append(opentxt(fileout0,info,name,hour,stmp,mm,vtype,subcloud),tmp)
              print 'LENGTH TMP TXT: ',len(tmp)
              if len(tmp)>1:
                tmp = np.abs(tmp[0]/tmp[1])
              data = np.append(data,tmp)
            else:
              # Open netcdf
              tmp = []
              for stmp in scat: 
                tmp  = np.append(opennetcdf(pathncdf,info,name,hour,stmp,var1D,subcloud),tmp)
              print 'LENGTH TMP NC: ',len(tmp)
              if len(tmp)>1:
                tmp = np.abs(tmp[0]/tmp[1])
              data = np.append(data,tmp)

          if colorplot is not None:
            if colorplot in vartxt:
              color  = opentxt(fileout0,info,name,hour,colorplot,mm,vtype,subcloud)
            else:
              # Open netcdf
              var1D   = False
              if colorplot in name1D:
                var1D = True 
              color  = opennetcdf(pathncdf,info,name,hour,colorplot,var1D,subcloud)
            colorall += [color[0]]
            color    = None
            print 'color : ',colorall
          else:
            color = info['color']

          #print name
          xx += [data[0]]
          yy += [data[1]]
          print scatter[0],xx
          print scatter[1],yy

          marker = 'o'
          if 'marker' in info.keys():
             marker=info['marker']
          
          # Figure
          ax.plot(data[0],data[1],marker=marker,markersize=5, markerfacecolor='none',color=info['color'],label=name)
          ax.text(data[0],data[1],hour,color=info['color'],fontsize=10)

      ii = 0
      for nn,name in enumerate(names):
        info  = infocas[name]
        for hh,hour in enumerate(info['hour']):
          print name,hour,xx[ii],yy[ii]
          ii += 1
          
      # stats
      idx     = (~np.isnan(xx) & ~np.isnan(yy))
      print idx[:],type(xx[:])
      xx      = np.array(xx)[idx]
      yy      = np.array(yy)[idx]
      corr    = np.corrcoef(xx,yy)[0,1]
      #### Confidence interval slope
      p       = np.polyfit(xx,yy, 1) 
      #print 'p ',p
      y_model = tl.equation(p, xx)
      print corr,p

      
      slope='y={aa}*x+{bb} (r={rr})'.format(aa="%4.2f" % p[0],bb="%4.2f" % p[1],rr="%4.2f" % corr)
      ax.plot(xx,y_model,lw=1,color='k',label=slope)
      if colorplot is not None:
        z = colorall[:] #-max(colorall)
        sc = ax.scatter(xx,yy,c=z, s=100, cmap=cm) # vmin=min(colorall), vmax=max(colorall)
        cbar   = plt.colorbar(sc)

      xmax = findmax(scatter[0])
      ymax = findmax(scatter[1])
      plt.xlim(xmax);plt.ylim(ymax)
      xlim1, xlim2 = plt.xlim()
      ylim1, ylim2 = plt.ylim()
 
      #### plot 1:1 line
      nb   = 100.
      x11  = np.arange(xlim1, xlim2,(xlim2-xlim1)/nb)
      y11  = [np.sign(p[0])*ij for ij in x11]
      ax.plot(x11,y11,lw=1,color='k',ls=':')

      # remove duplicate in legend
      handles, labels = plt.gca().get_legend_handles_labels()
      by_label = OrderedDict(zip(labels, handles))
      plt.legend(by_label.values(), by_label.keys())

      xscat=scatter[0];yscat=scatter[1]
      print xscat
      if len(xscat)>1:
        xscat=['Ratio'+''.join(scatter[0])]
      if len(yscat)>1:
        yscat=['Ratio'+''.join(scatter[1])]
      plt.xlabel(xscat[0])
      plt.ylabel(yscat[0])
      #plt.show()
      subname = ''
      if subcloud:
        subname = '_scld'
      name   = xscat[0]+'_'+yscat[0]+subname
      if colorplot is not None:
        name+= '_'+colorplot
      xsize  = (10,10)
      mf.savefig(fig,ax,pathout,name=name,xsize=xsize)
     
path     = "/Users/florentbrient/Dropbox/GitHub/scripts_LES/" #"/cnrm/tropics/user/brientf/MESONH/"
#relab   = True # by default

scat1   = ['DZe'] #'we'] #'Wstar' #'H0' #'THVdiff' #'Bowen' 'Q0' #'Kappa' #'fluxLH' #'THS1diff' #'THLMdiff'  #'PBLhRH' #'MSEdiff' #'THTdiff' 'gradTHT' 'gradRH' 'gradTHLM'
scat2   = ['we'] #skew_WT #['NNmax','NNmin'] #'vol_up' #'vol_up' 'WTm_down' #'rel_down_hum' #'neg_down_th' #'WTm_down' #'rel_down_th'
scatter = [scat1,scat2] #x-axis,y-axis

vartxt  = ['WTm_down','WTm_up','rel_up_hum','rel_up_th','rel_down_hum','rel_down_th',\
           'neg_down_th','pos_down_th','neg_down_hum','pos_down_hum',\
           'neg_up_th'  ,'pos_up_th'  ,'neg_up_hum'  ,'pos_up_hum',\
           'pos_up_th','vol_up','vol_down','MF_up','MF_down',\
           'WTm_down_1']
name1D  = ['E0','Q0','H0','Bowen','Wstar','Tstar','Ra','W\*','BL_H'] #'Kappa',

subcloud= True

colorplot=None #'WTm_up' #'H0'

main(path,infocas,scatter,name1D=name1D,vartxt=vartxt,subcloud=subcloud,colorplot=colorplot)



