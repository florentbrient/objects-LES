#import netCDF4 as nc
import pylab as plt
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
#import matplotlib.colors as colors
import numpy as np 
#import math as mp
from copy import deepcopy
#from collections import OrderedDict
import time
import tools as tl
#matplotlib.use('Agg')


def coloryourself(start,end,nb):
   color1=np.linspace(start[0],end[0],nb)
   color2=np.linspace(start[1],end[1],nb)
   color3=np.linspace(start[2],end[2],nb)
   colors=np.vstack((color1,color2,color3))
   return colors

def savefig(fig,ax,pathfig,title='',fts=15,xsize=(4,5),zmax=None,dpi=1000,quality=95,makepdf=False):
    plt.xticks(size=fts)
    plt.yticks(size=fts)
    plt.title(title)
    #pathfig = pathfig+'d_'+name+'/';tl.mkdir(pathfig)
    fig.set_size_inches(xsize)
    print(pathfig,title)
    fig.savefig(pathfig+title+'.png')#,dpi='figure')
    if makepdf:
      fig.savefig(pathfig+title+'.pdf')#,quality=quality)#, dpi=dpi)
    plt.close()

def opendata(data):
    print(data.keys())
    var           =  list(data.keys())[0]
    datav         = {}
    datav['Mean'] = data[var]

    # Anomaly of data1
    ip = -1
    if len(data[var].shape)==3:
       ip = 0
    datav['Anom'] = tl.anomcalc(datav['Mean'],ip=ip)
    return var,datav

#####

# plotmean()
def plotmean(data1D,data,dataobjs={},zz=None,pathfig='./',nametitle=None,fluxes=0):
    # plot standard deviation in the figure
    plotstd = False
    
    # Plot domain mean
    zyx  = [data1D[ij] for ij in data1D.keys()]

    styles  = ['-','--','-.',':']
    zzplot  = {'cloudbase':'--','cloudtop':'--'
             ,'lcl':'-','lfc':'-','lnb':'-'}
    zzcolor = {'cloudbase':'k','cloudtop':'k'
             ,'lcl':'gray','lfc':'gray','lnb':'gray'}
    lw      = 2
    lw2     = 0.5
    fts     = 12.0

    meanch    = ['Mean','Anom']
    if len(dataobjs)==0:
      meanch.remove('Anom')

    case      = pathfig.split('/')[-3]
    print('case : ',case)
    #print(data)
    var,datav = opendata(data)
    nz,ny,nx  = datav[meanch[0]].shape

    # Adjust objects:
    nbobjmax = 0 # 0 = do not plot each object
    #nbobjmax = 20 # plot each object if nbobj < nbobjmax
    #varmax   = ['WT','THV']

    # plot error bar every Di layers
    Di       = 4

    frac = {}; slctobjs={}; stdobjs ={}; numobj={}; nbobj={}
    time1  = time.time()
    if dataobjs is not None:
      for iobj,obj in enumerate(dataobjs.keys()):
        if dataobjs[obj] is not None:
          #print dataobjs[obj].shape
          zm           = deepcopy(dataobjs[obj])
          numobj[obj]  = np.unique(zm[zm>0]); #print 'OBJ',numobj[obj],len(numobj[obj]),np.max(zm)
          nbobj[obj]   = len(numobj[obj])
          slctobjs[obj]= {}; stdobjs[obj]= {}; 
          zm[zm>0]     = 1
          zm           = np.ma.masked_where(zm==0,zm)
          frac[obj]    = np.nansum(zm,axis=(1,2,))/float(ny*nx)
          for ij in meanch:
            #print obj,ij,np.nansum(zm),ny,nx
            tmp               = np.ma.masked_where(zm==0,datav[ij])
            slctobjs[obj][ij] = np.nanmean(tmp,axis=(1,2,)) #  datav[ij][zm],axis=(1,2,))
            stdobjs[obj][ij]  = np.nanstd(tmp,axis=(1,2,)) 
    time2  = time.time()
    print('Loop 1 took %0.3f ms' % ((time2-time1)*1000.0)) # 48 sec
     
     
    datav2       = {}
    for ij in meanch:
      datav2[ij] = np.nanmean(np.mean(datav[ij],axis=2),axis=1)

    # Name figures
    nameplot  = 'prof'
    nameTTT   = nameplot
    nametitle = nametitle.replace('TTTT',nameTTT)

    nametitle = nametitle.replace('AAAA','')
    
    pathfig += 'd_'+nameplot+'/'
    tl.mkdir(pathfig)

    for ij,field in enumerate(meanch):
      # start figure
      fig   = plt.figure()
      ax    = fig.add_subplot(111)
      title =  nametitle.replace('MMMM',field)

      # infos for figures
      infofig   = tl.infosfigures(case,var,mtyp=field)
#      levels  = infofig['levels'] # use it?
      zminmax = infofig['zminmax']
      
      if not(var=='Frac'):
        ax.plot(datav2[field],zyx[0],color='k',linestyle=styles[1],linewidth=lw*1.5)

      # plot objects
      if dataobjs is not None:
        for iobj,obj in enumerate(dataobjs.keys()):
          if dataobjs[obj] is not None:
            color    = tl.findhatch(obj,typ='color')
            tmpobj   = slctobjs[obj][field]
            #if field == 'Anom':  # Only for Anom?
            if fluxes or var=='Frac':
              tmpobj = tmpobj*frac[obj]
              
            if plotstd:
              #print 'every line : ',nbobj[obj],nbobjmax
              # Plot every line
              
              tmp0 = np.zeros((nbobj[obj],nz))
              time1 = time.time()
              for iobs,fobs in enumerate(numobj[obj]):
                ####tmp0[iobs] = np.nanmean(datav[field][dataobjs[obj]==fobs],axis=(1,2,))
                ####print iobs,fobs
                tmp          = np.ma.masked_where(dataobjs[obj]!=fobs,datav[field])
                tmp0[iobs,:] = np.nanmean(tmp,axis=(1,2,)) #  datav[ij][zm],axis=(1,2,))
                if nbobj[obj]<=nbobjmax and field=='Anom':
                  ax.plot(tmp0[iobs,:],zyx[0],color=(0.6,0.6,0.6),linestyle=styles[3],linewidth=lw/4.0)
              time2 = time.time()
              print('Loop 2 for %0.1f took %0.3f ms' % (len(numobj[obj]),(time2-time1)*1000.0))
              #260 sec for 53 objects
              stdobj = np.nanstd(tmp0,axis=0)
               
            else: 
              stdobj = stdobjs[obj][field]*0.
              
            # plot only some interquartile
            ax.errorbar(tmpobj[::Di],zyx[0][::Di],xerr=stdobj[::Di],fmt='none',linestyle=styles[0],ecolor=color,color=color,lw=lw/2.0)

            ax.plot(tmpobj,zyx[0],color=color,linestyle=styles[0],linewidth=lw*1.5)



      # plot z lines (cloud base, cloud top, lcl, lfc)
      if zz is not None:
        for iz,line in enumerate(zzplot.keys()):
          if line in zz.keys():
            ax.axhline(y=zyx[0][zz[line]], color=zzcolor[line], linewidth = 0.5, linestyle=zzplot[line])
     
      print('zminmax ',zminmax)
      if zminmax is not None:
         ax.set_ylim(zminmax)

      axes = plt.axis();# print axes
      if (np.sign(axes[0])!=np.sign(np.nanmax(axes[1]))):
        ax.axvline(x=0, color='k', linewidth=1)
      # adjust figure
      tl.adjust_spines(ax, ['left', 'bottom'])
      ax.get_yaxis().set_tick_params(direction='out')
      ax.get_xaxis().set_tick_params(direction='out')
      # Save figures
      savefig(fig,ax,pathfig,title=title,fts=fts)


# cross
def plot2D(data1D,data,dataobjs,xy=None,zz=None,pathfig='./',nametitle=None,avg=0,hatchlog=False,idxzi2D=None):
    # Plot cross section
    #      z-view

    if xy is not None:
      nameplot='cross'
      slct    = (0,xy)
      suffplot='x'+str(xy[0,0])+str(xy[0,1])+'_y'+str(xy[1,0])+str(xy[1,1])
      xsize   = (12,6)
    elif zz is not None:
      nameplot='zview'
      slct    = (1,zz)
      suffplot='z'+str(zz)
      xsize   = (12,9)
    else:
      print('Problem with plotcross')
      
    var,datav = opendata(data)

    case      = pathfig.split('/')[-3]
    print('case : ',case)

    zyx       = [data1D[ij] for ij in data1D.keys()]

    lw,lw2    = 2.,0.5
    fts       = 15
    meanch    = ['Mean','Anom']

    # Specific type of plot (cross, z-view)

    # Adjust objects:
    zm = {}; slctobjs={}
    if len(dataobjs)>0: #dataobjs is not None:
      # Remove 'All'
      objects = deepcopy(dataobjs)
      print(objects.keys())
      objects.pop('All')

      for iobj,obj in enumerate(objects.keys()):
        if objects[obj] is not None:
          #print objects[obj].shape
          slctobjs[obj],axisobj = tl.selectdata(zyx,objects[obj],slct=slct,avg=0)
          slctobjs[obj][slctobjs[obj]>0] = 1
          zm[obj]               = np.ma.masked_where(slctobjs[obj]==0,slctobjs[obj])
          

    # Name figures
    nameTTT   = nameplot+'_'+suffplot
    nametitle = nametitle.replace('TTTT',nameTTT)

    avgch = ''
    if avg != 0:
     avgch = '_avg'+str(avg)
    nametitle = nametitle.replace('AAAA',avgch)
    
    # Save figures
    pathfig += 'd_'+nameplot+'/'
    tl.mkdir(pathfig)

    for ij,field in enumerate(meanch):
      title =  nametitle.replace('MMMM',field)
      
      # infos for figures
      infofig   = tl.infosfigures(case,var,mtyp=field)
      cmaps     = infofig['cmap']
      levels,zminmax = None,None
      levels  = infofig['levels']
      if xy is not None:
        zminmax = infofig['zminmax']
      #print('levels: ',xy,levels)

      # Data to plot
      datap,axis = tl.selectdata(zyx,datav[field],slct=slct,avg=avg)
      # LCL altitude
      if idxzi2D is not None:
        idxzi     = np.repeat(idxzi2D[np.newaxis, :,:],len(axis[0]),axis=0)
        idxzi,tmp = tl.selectdata(zyx,idxzi,slct=slct,avg=avg)
        #print idxlcl[0,:]
        #print axis[0],axis[0][62]
        ZI        = [axis[0][int(idx)] for idx in idxzi[0,:]]
        print(ZI)
        
        
      #print('axis : ',axis[1])
      
      # Start plot 
      fig    = plt.figure(figsize=xsize)
      ax     = fig.add_subplot(111)

      #img    = plt.imread("water_worldview.jpg")
#      img    = plt.imread("dark-blue-ocean.jpg")
#      ax.imshow(img, extent=[0, max(axis[0]), 0, max(axis[1])]) #, aspect='auto')

      x0,y0 = axis[1],axis[0]
      if xy is not None:
          if xy[0,0] != xy[0,1]:
              x0 = [ij+xy[0,0] for ij in x0]
          elif xy[1,0] != xy[1,1]:
              x0 = [ij+xy[1,0] for ij in x0]
      xx, yy = np.meshgrid(x0,y0)

      cmap = cmaps['Mean']
      norm = None
      nmin,nmax = np.nanmin(datap),np.nanmax(datap)
      if levels is not None:
        nmin,nmax = np.nanmin(levels),np.nanmax(levels)
      if (np.sign(nmin)!=np.sign(nmax)) and nmin!=0:
        cmap  = cmaps['Anom']
        norm = mpl.colors.Normalize(vmin=nmin, vmax=abs(nmin))

      CS     = ax.contourf(xx,yy,datap,levels=levels,cmap=cmap,norm=norm) #,extend='both')

      if len(dataobjs)>0: #objects is not None:
        for iobj,obj in enumerate(objects.keys()):
          if objects[obj] is not None:
            if hatchlog:
              hatch = tl.findhatch(obj,typ='hatch')
              CS2   = ax.pcolor(xx,yy,zm[obj],hatch=hatch,alpha=0.)
            else:
              color = tl.findhatch(obj,typ='color')
              #print 'color ',color,obj,slctobjs[obj].any(),np.nanmax(slctobjs[obj])
              if slctobjs[obj].any():
                CS2   = ax.contour(xx,yy,slctobjs[obj],colors=color,linewidths=lw/4.)
                #CS2   = ax.contour(xx,yy,zm[obj],colors=color,linewidths=lw/4)
      else:
        if idxzi2D is not None:
          #CS3 = plt.plot(axis[1],ZI,color='k',lw=lw2,linestyle='--')
          CS3 = plt.plot(x0,ZI,color='k',lw=lw2,linestyle='--')
        elif (np.sign(nmin)!=np.sign(nmax)):
          CS3 = plt.contour(xx,yy,datap,levels=[0],colors='k',lw=lw2,linestyles='--')
        

      cbar   = plt.colorbar(CS,)
      cbar.ax.tick_params(labelsize=fts)

      if zminmax is not None:
         ax.set_ylim(zminmax)
      savefig(fig,ax,pathfig,title=title,fts=fts,xsize=xsize)



def plotstats(datax,datay,data2=None,data3=None,logx=False,filesave=None,xlab=None,ylab=None,unitx=None,xaxis=None,yaxis=None,size=[10.0,5.0],extrax=None,extray=None,nbmin=None,mins=None,mins1=None,mins2=None):
# subroutine to plot the evolution of object characteristics with threshold defining objects
# related to the routine plot_stats.py
    colors  = ['b','r','purple','g','k']
    markers = ['o','x','s','d']
    lw      = 2
    fts     = 15
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    sz  = datax.shape; print(sz)
    ms  = 70
    if nbmin is not None:
      colors  = colorBltoRed(sz[0])
      ms      = [np.sqrt(n)/1.5 for n in nbmin]
      print nbmin
      idx     = np.where(nbmin==10000)[0]
      msextra = ms[idx]
    else:
      idx     = np.arange(0,sz[1]) 
      msextra = ms
    if logx:
      datax[datax==0.]=np.nan
      datay[datay==0.]=np.nan



    for ij in range(sz[0]):
      #if colors[0]=='r':
      color = colors[-1] #[ij]
      #else:
      #  color = colors[:,ij]
      #if sz[0]==1:
      #  color = [0.5,0,0.5]
        #print colors
      #l1, =ax.plot(datax[ij,:],datay[ij,:],color=color,linewidth=lw)
      #if logx:
      #  l1, =ax.semilogx(datax[ij,:],datay[ij,:],color=color,linewidth=lw)
      print datay.shape
      ax.scatter(datax[ij,:],datay[ij,:],marker=markers[0],s=ms,color=color)
      ax.plot(datax[ij,:],datay[ij,:],color=color,linewidth=lw)

      if extrax is not None and extray is not None:
        for ik in range(len(extrax)):
          ax.scatter(extrax[ik][ij,idx],extray[ik][ij,idx],marker=markers[1],s=msextra,color=colors[ik])
        #ax.scatter(extrax[0][ij,idx],extray[0][ij,idx],marker=markers[1],s=msextra,color=color)
        #ax.scatter(extrax[1][ij,idx],extray[1][ij,idx],marker=markers[2],s=msextra,color=color)
          if len(idx) != 1 :
            ax.plot(extrax[ik][ij,:],extray[ik][ij,:],color=colors[ik],linewidth=lw,linestyle='--')
          #ax.plot(extrax[0][ij,:],extray[0][ij,:],color=color,linewidth=lw,linestyle='--')
          #ax.plot(extrax[1][ij,:],extray[1][ij,:],color=color,linewidth=lw,linestyle='--')
        if mins1 is not None:
          ax.fill_between(datax[ij,:],mins1[ik][0,:],mins1[ik][1,:],color='grey',alpha=0.3)
        #if mins2 is not None:
        #  ax.fill_between(datax[ij,:],mins2[0,:],mins2[1,:],color='grey',alpha=0.3)

      if mins is not None:
        print mins
        ax.fill_between(datax[ij,:],mins[0,:],mins[1,:],color='grey',alpha=0.3)

        #l1, =ax.plot(extrax[0][ij,:],extray[0][ij,:],color=color,marker=markers[1],linewidth=lw)
        #if logx:
        #  l1, =ax.semilogx(extrax[0][ij,:],extray[0][ij,:],color=color,marker=markers[1],linewidth=lw)

      ax.axis('tight')
      if logx:
        ax.set_xscale('log')

      if xaxis is not None:
        if 'time' in xaxis.keys():
          xx = xaxis[xlab]
          barloc = float(xx[-1])-(xx[1]-xx[0])
          #print barloc+0.3*(xx[1]-xx[0]) # used to be 24
          Delta = 0.2 #1.0/float(len(extrax))
          ax.errorbar(barloc+Delta*(xx[1]-xx[0]),np.nanmean(datay,axis=1),yerr=np.std(datay,axis=1),marker=markers[0],markersize=np.sqrt(ms),color='k')
          if extrax is not None and extray is not None:
            for ik in range(len(extrax)):
              #print 'ik: ',ik,barloc+(ik*Delta)*(xx[1]-xx[0]),colors[ik]
              ax.errorbar(barloc+((ik+2)*Delta)*(xx[1]-xx[0]),np.nanmean(extray[ik],axis=1),yerr=np.nanstd(extray[ik],axis=1),marker=markers[1],markersize=np.sqrt(ms),color=colors[ik])
              ax.scatter(extrax[0][ij,idx],extray[0][ij,idx]+extray[1][ij,idx],marker=markers[0],s=msextra,edgecolors='k',color='none')
           
            #ax.errorbar(barloc+0.3*(xx[1]-xx[0]),np.mean(extray[0],axis=1),yerr=np.std(extray[0],axis=1),marker=markers[1],markersize=np.sqrt(ms),color='k')
            #ax.errorbar(barloc+0.7*(xx[1]-xx[0]),np.mean(extray[1],axis=1),yerr=np.std(extray[1],axis=1),marker=markers[2],markersize=np.sqrt(ms),color='k')
            #print filesave,np.mean(datay,axis=1),np.mean(extray[0],axis=1),np.mean(extray[1],axis=1)
            #print filesave,np.std(datay,axis=1),np.std(extray[0],axis=1),np.std(extray[1],axis=1)
            # plot extra circle for downdraft + updraft
            ax.scatter(extrax[0][ij,idx],extray[0][ij,idx]+extray[1][ij,idx],marker=markers[0],s=msextra,edgecolors='k',color='none')

        
    if yaxis is not None:
      if ylab in yaxis.keys():
        #ax.set_ylim([yaxis[0][ylab],yaxis[1][ylab]])
        plt.ylim(ymin=min(yaxis[ylab]),ymax=max(yaxis[ylab]))
        ax.axhline(y=0,   color='k', linewidth=1, linestyle='--')
        #plt.ylim()
        #ax.axhline(y=100, color='k', linewidth=1, linestyle='--')
        if xaxis is not None:
          if xlab in xaxis.keys() and max(yaxis[ylab])==max(xaxis[xlab]):
            plt.plot([min(yaxis[ylab]),max(yaxis[ylab])],[min(yaxis[ylab]),max(yaxis[ylab])],'k--')

    if xaxis is not None: # and not logx:
      if xlab in xaxis.keys():
        if not logx:
          plt.xticks(xaxis[xlab],size=fts)
        else:
          plt.xticks(size=fts)
        ax.set_xlim([min(xaxis[xlab]),max(xaxis[xlab])])
       


    ax.axvline(x=0, color='k', linewidth=1)
    tl.adjust_spines(ax, ['left', 'bottom'])
    ax.get_yaxis().set_tick_params(direction='out')
    ax.get_xaxis().set_tick_params(direction='out')
    ax.set_xlabel(xlab+'['+unitx+']',fontsize=fts)
    ax.set_ylabel(ylab,fontsize=fts)
    plt.xticks(size=fts)
    plt.yticks(size=fts)
    fig.set_size_inches(size[0], size[1])
    fig.savefig(filesave + '.png')
    fig.savefig(filesave + '.pdf')
    plt.close()

