from collections import OrderedDict

infocas  = OrderedDict()
infocas['FIRE']    = {'sens':'Ls2x0','hour':('006','012','018','024',),'vmin':800 ,'color':'orange','dxdy':50,'svt':True,'cas':'FIRE'}
infocas['AYOT24S']  = {'sens':'24Sx0','hour':('006',),'vmin':3200,'color':'gray','dxdy':25,'svt':True,'cas':'AYOTTE'}
infocas['IHOP']    = {'sens':'Ru0x0','hour':('004','006','008',),'vmin':3200,'color':'g','dxdy':25}
infocas['ARMCU']   = {'sens':'Ru0x0','hour':('006','008','010',),'vmin':3200,'color':'k','dxdy':25}
infocas['BOMEX']   = {'sens':'Ru0x0','hour':('008',),'vmin':3200,'color':'r','dxdy':25}
infocas['RICO']    = {'sens':'Ru0x0','hour':('008','012',),'vmin':3200,'color':'b','dxdy':25}
infocas['ASTEX']   = {'sens':'BIG07','hour':('006','010','014','018',),'vmin':50,'color':'c','dxdy':160,'cas':'ASTEX'}
infocas['FIREnoc'] = {'sens':'LNOCF','hour':('006','012',),'vmin':800 ,'color':'y','marker':'s','dxdy':50,'cas':'FIRE'}
infocas['IHODC']   = {'sens':'IHODC','hour':('006',),'vmin':3200,'color':'darkgreen','marker':'D','dxdy':25\
                     ,'svt':True,'cas':'IHOP','vtyp':'RK4DI','nc':'nc','name1':'___'}
#infocas['FIREopen']= {'sens':'L2RA2','hour':('006','012','018','024',),'vmin':800 ,'color':'gold','marker':'s','dxdy':50,'svt':True,'cas':'FIRE'}
infocas['BOMEXnoc']= {'sens':'Ru0NC','hour':('008','012',),'vmin':3200,'color':'red','marker':'s','dxdy':25,'cas':'BOMEX'}
infocas['ARMCUnoc']= {'sens':'Ru0NC','hour':('006','008','010',),'vmin':3200,'color':'k','marker':'s','dxdy':25,'cas':'ARMCU'}
infocas['AYOT03S']  = {'sens':'03Su0','hour':('006',),'vmin':3200,'color':'gray','dxdy':25,'svt':True,'marker':'s','cas':'AYOTTE','short':True}
infocas['AYOT05S']  = {'sens':'05Su0','hour':('006',),'vmin':3200,'color':'gray','dxdy':25,'svt':True,'marker':'D','cas':'AYOTTE','short':True}
infocas['AYOT05W']  = {'sens':'05Wu0','hour':('006',),'vmin':3200,'color':'gray','dxdy':25,'svt':True,'marker':'x','cas':'AYOTTE','short':True}
