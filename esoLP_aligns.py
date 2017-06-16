#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 17:36:20 2016

@author: mmarsset
"""

#! /usr/bin/env python


import imp
try:
    imp.find_module('astropy')
    astropyFound = True
except ImportError:
    astropyFound = False
if astropyFound:
    from astropy.io import fits as pyf
else:
    import pyfits as pyf
try:
    imp.find_module('mp_ephem')
    orbfitFound = True
except ImportError:
    orbfitFound=False
if orbfitFound:
    from mp_ephem.ephem import EphemerisReader
    from mp_ephem.bk_orbit import BKOrbit 
else:
    import mpcTools


import numpy as num
import pylab as pyl
import stsci.numdisplay as numdisplay
from astLib import astWCS as WCS,astImages as Im
from trippy.scamp import scamp, makeParFiles
from matplotlib.patches import Circle
import sys,glob,os
from trippy import bgFinder
import optparse
from numpy import ma
import string

def trimCatWithMask(catalog,mask,trimWidth=.5,trimRound=.3):#trimWidth=1.,trimRound=0.5):

    newCatalog={}
    keys=catalog.keys()
    good=[]

    good=num.where((catalog['XWIN_IMAGE']>=10) & (catalog['XWIN_IMAGE']<=1020) & 
        (catalog['YWIN_IMAGE']>=10) & (catalog['YWIN_IMAGE']<=1020) & (catalog['AWIN_IMAGE']>trimWidth) & 
        (catalog['BWIN_IMAGE']/catalog['AWIN_IMAGE']>trimRound) & (catalog['FLUX_AUTO']>0))
    
    print "\nAccepting %s good sources.\n"%(len(good[0]))
    for ii in keys:
        newCatalog[ii]=catalog[ii][good]
    return newCatalog


class Cursor:
    def __init__(self, ax,drawax):
        self.ax = ax
        self.drawax=drawax
        self.lx = drawax.axhline(color='r',linestyle='--')  # the horiz line
        self.ly = drawax.axvline(color='r',linestyle='--')  # the vert line

        self.lxl = ax.axhline(color='r',linestyle='--')  # the horiz line
        self.lyl = ax.axvline(color='r',linestyle='--')  # the vert line

        # text location in axes coords
        self.txt = ax.text( 0.7, 0.9, '', transform=ax.transAxes)

    def mouse_move(self, event):
        if not event.inaxes: return

        x, y = event.xdata, event.ydata
        # update the line positions
        self.lx.set_ydata(y )
        self.ly.set_xdata(x )
        self.lxl.set_ydata(y )
        self.lyl.set_xdata(x )

        self.txt.set_text( 'x=%1.2f, y=%1.2f'%(x,y) )
        pyl.draw()


def JD2Greg(inp,fracyear=False,MJD=False):

    x=str(inp)
    try:
        jd=float(x)
    except ValueError:
        pass
    if MJD: jd+=2400000.5
    jd=jd+0.5
    Z=int(jd)
    F=jd-Z
    alpha=int((Z-1867216.25)/36524.25)
    A=Z + 1 + alpha - int(alpha/4)

    B = A + 1524
    C = int( (B-122.1)/365.25)
    D = int( 365.25*C )
    E = int( (B-D)/30.6001 )

    dd = B - D - int(30.6001*E) + F

    if E<13.5:
        mm=E-1

    if E>13.5:
        mm=E-13

    if mm>2.5:
        yyyy=C-4716

    if mm<2.5:
        yyyy=C-4715

    months=["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "Dec\
ember"]
    months=["01","02","03","04","05","06","07","08","09","10","11","12"]
    daylist=[31,28,31,30,31,30,31,31,30,31,30,31]
    daylist2=[31,29,31,30,31,30,31,31,30,31,30,31]

    h=int((dd-int(dd))*24)
    min=int((((dd-int(dd))*24)-h)*60)
    sec=86400*(dd-int(dd))-h*3600-min*60

    # Now calculate the fractional year. Do we have a leap year?                                                               
    if (yyyy%4 != 0):
        days=daylist2
    elif (yyyy%400 == 0):
        days=daylist2
    elif (yyyy%100 == 0):
        days=daylist
    else:
        days=daylist2

    if fracyear:
        daysum=0
        for y in range(mm-1):
            daysum=daysum+days[y]
        daysum=daysum+dd-1

        if days[1]==29:
            fracyear=yyyy+daysum/366
        else:
            fracyear=yyyy+daysum/365
        return fracyear
    else:
        sec=float(int(sec*10.0)/10.0)
        return (str(yyyy)+'/'+months[mm-1]+'/'+str(int(dd)),string.zfill(h,2)+":"+string.zfill(min,2)+":"+string.zfill(sec,4))



def getShift(event):
    global catalog1,catalog2,sp1,sp2,pathList1,patchList2,selected1,selected2,starx1,stary1,starx2,stary2

    (x,y)=event.xdata,event.ydata
    #print x,y
    #print selected1,selected2
    if event.inaxes==sp1 and selected1: 
        patchList1[len(patchList1)-1].remove()
        selected1=False
        patchList2[len(patchList2)-1].remove()
        selected2=False
    
    d=((catalog1['XWIN_IMAGE']-x)**2+(catalog1['YWIN_IMAGE']-y)**2)**0.5
    arg=num.argmin(d)
    starx,stary=catalog1['XWIN_IMAGE'][arg],catalog1['YWIN_IMAGE'][arg]
    circle=Circle((starx,stary),15,facecolor="none",edgecolor='red',linestyle='solid',linewidth=2, alpha=0.75,zorder=10)

    sp1.add_patch(circle)
    patchList1.append(circle)
    selected1=True
    starx1,stary1=starx,stary

    d=((catalog2['XWIN_IMAGE']-x)**2+(catalog2['YWIN_IMAGE']-y)**2)**0.5
    arg=num.argmin(d)
    starx,stary=catalog2['XWIN_IMAGE'][arg],catalog2['YWIN_IMAGE'][arg]
    circle=Circle((starx,stary),15,facecolor="none",edgecolor='red',linestyle='solid',linewidth=2, alpha=0.75,zorder=10)


    sp2.add_patch(circle)
    patchList2.append(circle)
    starx2,stary2=starx,stary

    pyl.draw()


def computeOffsets(event):
    global catalog1,catalog2,selected1,selected2,starx1,stary1,starx2,stary2,opt,ox,oy,eox,eoy,maxShift
    
    if event.key==' ':
        if not selected1:
            print 'Need to select the common source in BOTH frames.'
            return

        dx=starx2-starx1
        dy=stary2-stary1
        #print dx,dy

        X,Y=catalog2['XWIN_IMAGE']*1.,catalog2['YWIN_IMAGE']*1.
        X-=dx
        Y-=dy
        dxs=[]
        dys=[]
        for ii in range(len(catalog1['XWIN_IMAGE'])):
            x,y=catalog1['XWIN_IMAGE'][ii],catalog1['YWIN_IMAGE'][ii]
            d=((X-x)**2+(Y-y)**2)**0.5
            arg=num.argmin(d)
            if d[arg]<opt.maxDist:
                #print len(dxs)+1,x,y,catalog2['XWIN_IMAGE'][arg],catalog2['YWIN_IMAGE'][arg]
                dxs.append(catalog2['XWIN_IMAGE'][arg]-x)
                dys.append(catalog2['YWIN_IMAGE'][arg]-y)
        dxs=num.array(dxs)
        dys=num.array(dys)
        ox,eox=num.mean(dxs),num.std(dxs)*len(dxs)**-0.5
        oy,eoy=num.mean(dys),num.std(dys)*len(dys)**-0.5
        #print dxs
        if len(dxs)<1.:
            print 'Not enough stars to get a reliable shift.'
            ox,oy=0.0,0.0
            eox,eoy=-1.,-1.
        print 'Offsets %s,%s with errors %s %s'%(ox,oy,eox,eoy)
        pyl.close()




    elif event.key=='q':
        pyl.close()
        sys.exit()
    elif event.key=='d':
        print 'Accepting default offsets of 0.0,0.0'
        ox=0.0
        eox=-1.
        oy=0.0
        eoy=-1.
        pyl.close()

    
    

def getScaledImage(sec,A=501,B=501,nsamp=200,contrast=0.5,lims=[]):
    if lims==[]:
        lims=[0,A,0,B]
    im=num.zeros([A,B]).astype(float)
    (z1,z2)=numdisplay.zscale.zscale(sec,nsamples=nsamp,contrast=contrast)
    norm=Im.normalise(sec,[z1,z2])
    im[lims[0]:lims[1],lims[2]:lims[3]]=norm
    return im



if __name__=="__main__":

    usage="> esoLP_aligns.py [-d]"
    parser=optparse.OptionParser(usage=usage)
    parser.add_option('-d', dest='defaultShifts', default=False,
                      action='store_true',
                      help = 'Accept default shifts contained in image headers? DEFAULT=%default')
    parser.add_option('-t',
                      dest='target', action='store', default='O13BL3TA',
                      type='str',
                      help='The target name (and directory). DEFAULT=%default')
    parser.add_option('--maxDist',
                      dest='maxDist', type='float', default=5.,
                      action='store',
                      help = 'Maximum distance between user provided shift and location of target to accept as good for measuring the shift of an image. DEFAULT=%default')
    parser.add_option('--maskName',
                      default='bpMap.fits', action='store', dest='maskName',
                      type='str',
                      help='Bad pixel mask name. DEFAULT=%default.')
    parser.add_option('--sidereal',
                      default=False,dest='sidereal',action='store_true',
                      help='Process a sidereal target? DEFAULT=%default')
    parser.add_option('--siderealOut',
                      action='store',default='shifted.fits', type='str',
                      dest='siderealfn',
                      help='Sidereal output image name. DEFAULT=%default')
    parser.add_option('--nonsiderealOut',
                      action='store',default='nshifted.fits', type='str',
                      dest='nsiderealfn',
                      help='non-Sidereal output image name. DEFAULT=%default')
    parser.add_option('-m',
                      action='store_true',dest='manualShift', default=False,
                      help='Do a manual alignment rather than take the auto generated version.')
    parser.add_option('--skipSex',
                      action='store_true',dest='skipSex',default=False,
                      help='Skip running sextractor as you presumably have already done so.')
    parser.add_option('--maxShift',
                      dest='maxShift',action='store',type='float',
                      default=8.,
                      help='Maximum distance between nominal wcs-based shift and actual shift. DEFAULT=%default.')
    parser.add_option('--sexThreshold',
                      dest='sexThreshold', action='store', type='float',
                      default=4.0,
                      help="Sextractor SNR threshold when identifying bright sources to shift on. DEFAULT=%default.")
    parser.add_option('--sexArea',
                      dest='sexArea', action='store', type='float',
                      default=1.0,
                      help="Sextractor AREA threshold when identifying bright sources to shift on. DEFAULT=%default")
    parser.add_option('--imageList',
                      dest='imageList',default=None,action='store',
                      help='Provide a list file name to load images from a list rather than by target name.')
     
    (opt,args)=parser.parse_args()
    maxShift=100.#opt.maxShift

    with pyf.open(opt.target+'/CALIB/'+opt.maskName) as dhan:
        bpMask=dhan[0].data
        #bpMask=bpMask*0.+1.
    
    os.chdir(opt.target+'/REDUCED')
    

    writeParFiles=True
    if writeParFiles:
        os.system('rm michascience.sex')
        makeParFiles.writeSex('michascience.sex',
                              minArea=opt.sexArea,
                              threshold=opt.sexThreshold,
                              zpt=27.0,
                              aperture=15.,
                              catalogType='FITS_LDAC',
                              saturate=15000.)
        makeParFiles.writeConv()
        makeParFiles.writeParam(numAps=1)


        
    if opt.imageList<>None:
        science=[]
        with open(opt.imageList) as han:
            data=han.readlines()
        for i in range(len(data)):
            science.append(data[i].split()[0])
    else:
        science=glob.glob('rE*fits')

        

    wheaders=[]
    headers=[]
    datas=[]
    catalogs=[]
    for i in range(len(science)):
        print 'Loading data for',science[i]
        with pyf.open(science[i]) as dhan:
            headers.append(dhan[0].header)
            datas.append(dhan[0].data)
            data = dhan[0].data
            
    
    data = data*0+1
    data[0:1030,0:1030] = bpMask
    bpMask = data

    
    #get target relative positions in each image (assuming PA=0 for now)
    if not opt.sidereal:
        wcs=WCS.WCS(headers[i],mode='pyfits')

        obs = EphemerisReader().read(filename='../../../mpcs/%s.mpc'%(opt.target))
        orb = BKOrbit(obs)
            
        MJD=headers[0]['MJD-OBS']
        orb.predict(MJD+2400000.5)
        ra1,dec1=orb.coordinate.ra.degree,orb.coordinate.dec.degree
        (x1,y1)=wcs.wcs2pix(ra1,dec1)
            
        MJD=headers[0]['MJD-OBS']+1./24.
        orb.predict(MJD+2400000.5)
        ra2,dec2=orb.coordinate.ra.degree,orb.coordinate.dec.degree
        (x2,y2)=wcs.wcs2pix(ra2,dec2)
        
        dx=(x2-x1)*24.
        dy=(y2-y1)*24.

        rx,ry=[],[]
        for i in range(len(headers)):
            dMJD=headers[i]['MJD-OBS']-headers[0]['MJD-OBS']
            rx.append(dx*dMJD)
            ry.append(dy*dMJD)
        rx=num.array(rx)
        ry=num.array(ry)
        #print rx
    else:
        wcs=WCS.WCS(headers[i],mode='pyfits')
        (ra1,dec1)=wcs.pix2wcs(0,0)
        rx=num.zeros(len(science)).astype('float32')
        ry=num.zeros(len(science)).astype('float32')

    useWCS=True
    if useWCS:
        jx=[]
        jy=[]
        for i in range(len(headers)):
            wcs=WCS.WCS(headers[i],mode='pyfits')
            (x,y)=wcs.wcs2pix(ra1,dec1)
            jx.append(x)
            jy.append(y)
        DXs=-(num.array(jx)-jx[0])
        DYs=-(num.array(jy)-jy[0])
    else:

        DXs=[]
        DYs=[]
        for j in range(len(datas)):
            DX=-(headers[j]['POFFSET']-headers[0]['POFFSET'])/0.1171
            DY=(headers[j]['QOFFSET']-headers[0]['QOFFSET'])/0.1171
            DXs.append(DX)
            DYs.append(DY)
        DXs=num.array(DXs)
        DYs=num.array(DYs)

    origDXs=DXs*1.
    origDYs=DYs*1.
    #rx+=DXs
    #ry+=DYs
    #sys.exit()
    DXs-=num.min(DXs)
    DYs-=num.min(DYs)
    maxX=num.max(DXs)
    maxY=num.max(DYs)


    (A,B)=datas[0].shape
    w=num.where(bpMask==0) #the zeros are the bad pixels!
    bigDatas=num.zeros((len(datas),A+(maxY),B+(maxX))).astype('float32')

    for i in range(len(DXs)):
        bgf=bgFinder.bgFinder(datas[i])
        mode=bgf.smartBackground()
        #mode=jjkmode.stats(datas[i])[0]
        datas[i][w]=mode
        bigDatas[i]+=mode
        bigDatas[i,DYs[i]:DYs[i]+A, DXs[i]:DXs[i]+B]=datas[i]
        
        print 'Source extracting',science[i]
        if not opt.skipSex:
            scamp.runSex('michascience.sex','%s'%(science[i]),options={'CATALOG_NAME':'s%s.cat'%(science[i].split('.fits')[0])},verbose=True)
            pyf.writeto('s'+science[i],bigDatas[i],clobber=True,header=headers[i])
        
        
        catalog=scamp.getCatalog('s%s'%(science[i].split('.fits')[0]+'.cat'),paramFile='wes.param')
        catalog=trimCatWithMask(catalog.copy(),bpMask)


        for j in range(len(catalog['XWIN_IMAGE'])):
            catalog['XWIN_IMAGE'][j]+=DXs[i]
            catalog['YWIN_IMAGE'][j]+=DYs[i]
        catalogs.append(catalog.copy())
        
        

    OX,OY,eOX,eOY,dFlux=[0.],[0.],[0.],[0.],[1.]
    (junk,A,B)=bigDatas.shape
    
    for j in range(1,len(science)):

        if opt.defaultShifts:
            OX.append(0.0)
            OY.append(0.0)
            eOX.append(-1.)
            eOY.append(-1.)
            continue

        fig=pyl.figure(figsize=(30,15))
        fig.subplots_adjust(wspace=0)
        ox,eox=-32768,-32768
        oy,eoy=-32768,-32768
        
        selected1=False
        catalog1=catalogs[0].copy()


        selected2=False
        catalog2=catalogs[j].copy()
        
        if not opt.manualShift:
            #do an autosort on the good stars
            dist=[]
            for ii in range(len(catalog1['XWIN_IMAGE'])):
                d=((catalog1['XWIN_IMAGE'][ii]-catalog2['XWIN_IMAGE'])**2+(catalog1['YWIN_IMAGE'][ii]-catalog2['YWIN_IMAGE'])**2)**0.5
                arg=num.argmin(d)
                if d[arg]<maxShift:
                    #print catalog2['FLUX_AUTO'][arg]/catalog1['FLUX_AUTO'][ii],'&&'
                    dist.append([catalog2['XWIN_IMAGE'][arg]-catalog1['XWIN_IMAGE'][ii],
                                 catalog2['YWIN_IMAGE'][arg]-catalog1['YWIN_IMAGE'][ii],
                                 d[arg],
                                 arg,
                                 catalog2['FLUX_AUTO'][arg]/catalog1['FLUX_AUTO'][ii]])
                    
            dist=num.array(dist)
            print j, science[j]

            ox=num.median(dist[:,0])
            oy=num.median(dist[:,1])
            eox=num.std(dist[:,0])
            eoy=num.std(dist[:,1])
            dflux=num.mean(dist[:,4])

            print 'Offsets %1.4f,%1.4f with errors %1.3f %1.3f from %s sources with a mean flux variation of %1.3f.'%(
                ox,oy,eox,eoy,len(dist[:,0]),dflux)
            OX.append(ox)
            OY.append(oy)
            eOX.append(eox)
            eOY.append(eOY)
            dFlux.append(dflux)
            continue


        scaled1=getScaledImage(bigDatas[0],A=A,B=B)
        scaled2=getScaledImage(bigDatas[j],A=A,B=B)

        sp1=fig.add_subplot(121)
        sp1.set_aspect('equal','datalim')
        pyl.title(science[0])
        axIm1=pyl.imshow(scaled1*1.,cmap='gray',origin='lower')
        
        
        sp2=fig.add_subplot(122)
        sp2.set_aspect('equal','datalim')
        pyl.title(science[j])
        axIm=pyl.imshow(scaled2*1.,cmap='gray',origin='lower')
        
        
        patchList1=[]
        for i in range(len(catalog1['XWIN_IMAGE'])):
            x=catalog1['XWIN_IMAGE'][i]
            y=catalog1['YWIN_IMAGE'][i]
            circle=Circle((x,y),20,facecolor="none",edgecolor='blue',linestyle='dashed',linewidth=2, alpha=0.75,zorder=10)
            sp1.add_patch(circle)
            patchList1.append(circle)
        patchList2=[]
        for i in range(len(catalog2['XWIN_IMAGE'])):
            x=catalog2['XWIN_IMAGE'][i]
            y=catalog2['YWIN_IMAGE'][i]
            circle=Circle((x,y),20,facecolor="none",edgecolor='blue',linestyle='dashed',linewidth=2, alpha=0.75,zorder=10)
            sp2.add_patch(circle)
            patchList2.append(circle)
            
        pyl.connect('button_press_event',getShift)
        pyl.connect('key_press_event',computeOffsets)

        cursor = Cursor(sp1,sp2)
        pyl.connect('motion_notify_event', cursor.mouse_move)
        pyl.show()

        OX.append(ox)
        OY.append(oy)
        eOX.append(eox)
        eOY.append(eOY)


    dFlux=num.array(dFlux)
    dFlux/=num.mean(dFlux)

    #now repeat the offsets but with the new ideal shifts for the stationary target
    #print OX
    #print OY
    DXs=origDXs*1.
    DYs=origDYs*1.
    DXs=num.array(DXs)-num.array(OX)
    DYs=num.array(DYs)-num.array(OY)
    DXs-=num.min(DXs)
    DYs-=num.min(DYs)
    DXs=(DXs+0.5).astype(int)
    DYs=(DYs+0.5).astype(int)
    maxX=num.max(DXs)
    maxY=num.max(DYs)
 
    outHeader=headers[0]
    airMass=0.
    airMassStart=0.
    airMassEnd=0.
    mjdObs=0.
    if astropyFound:
        outHeader.set('CRPIX1',outHeader['CRPIX1']+DXs[0])
        outHeader.set('CRPIX2',outHeader['CRPIX2']+DXs[0])
        outHeader.set('NSTACKED',int(len(DXs)))
        outHeader.set('EXPTIME',headers[0]['EXPTIME'])
        outHeader.set('HIERARCH ESO DET OUT1 GAIN',headers[0]['HIERARCH ESO DET OUT1 GAIN'])
        outHeader.set('HIERARCH ESO DET OUT1 RON',headers[0]['HIERARCH ESO DET OUT1 RON'])
        for jj in range(len(headers)):
            airMassStart += (headers[jj]['HIERARCH ESO TEL AIRM START'])
            airMassEnd += (headers[jj]['HIERARCH ESO TEL AIRM END'])
            airMass = (airMassStart+airMassEnd)/2.
            mjdObs += (headers[jj]['MJD-OBS'])
        outHeader.set('AIRMASS',airMass/len(headers))#headers[len(headers)/2]['HIERARCH ESO TEL AIRM START'])
        outHeader.set('MJD-OBS',mjdObs/len(headers)+outHeader['EXPTIME']/24./3600./2.)#(headers[len(headers)/2]['MJD-OBS']+outHeader['EXPTIME']/24./3600./2.)
    else:
        outHeader.update('CRPIX1',outHeader['CRPIX1']+DXs[0])
        outHeader.update('CRPIX2',outHeader['CRPIX2']+DXs[0])
        outHeader.update('NSTACKED',int(len(DXs)))
        outHeader.update('EXPTIME',headers[0]['EXPTIME'])
        outHeader.update('HIERARCH ESO DET OUT1 GAIN',headers[0]['HIERARCH ESO DET OUT1 GAIN'])
        outHeader.update('HIERARCH ESO DET OUT1 RON',headers[0]['HIERARCH ESO DET OUT1 RON'])
        for jj in range(len(headers)):
            airMassStart += (headers[jj]['HIERARCH ESO TEL AIRM START'])
            airMassEnd += (headers[jj]['HIERARCH ESO TEL AIRM END'])
            airMass = (airMassStart+airMassEnd)/2.
            mjdObs += (headers[jj]['MJD-OBS'])
        outHeader.set('AIRMASS',airMass/len(headers))#headers[len(headers)/2]['HIERARCH ESO TEL AIRM START'])
        outHeader.set('MJD-OBS',mjdObs/len(headers)+outHeader['EXPTIME']/24./3600./2.)#(headers[len(headers)/2]['MJD-OBS']+outHeader['EXPTIME']/24./3600./2.)

    (A,B)=datas[0].shape

    bigDatas=num.zeros((len(datas),A+(maxY),B+(maxX))).astype('float32')
    masks=num.ones((len(datas),A+(maxY),B+(maxX))).astype('float32')
    cosmicMasks=num.zeros((len(datas),A,B)).astype('float32')
    modes=[]
    for i in range(len(DXs)):
        bgf=bgFinder.bgFinder(datas[i])
        mode=bgf.smartBackground()
        #mode=jjkmode.stats(datas[i])[0]

        #generate masks for cosmic ray strikes which saturate (produce negative pixels in high BG images)
        neg=num.where(datas[i]<0)
        for j in range(len(neg[0])):
            x=neg[1][j]
            y=neg[0][j]
            cosmicMasks[i:min(len(datas)-1,i+3),max(0,y-1):min(A,y+2),max(0,x-1):min(B,x+2)]=1 #do 3 images to get rid of persistence of cosmic rays

        bigDatas[i,DYs[i]:DYs[i]+A, DXs[i]:DXs[i]+B]=datas[i]-mode
        masks[i,DYs[i]:DYs[i]+A, DXs[i]:DXs[i]+B]=0.0
        masks[i,DYs[i]:DYs[i]+A, DXs[i]:DXs[i]+B]+=1.-bpMask #in the numpy array, good pixels are 0, in bpMask they are 1
        masks[i,DYs[i]:DYs[i]+A, DXs[i]:DXs[i]+B]+=cosmicMasks[i]

        neg = num.where( bigDatas[i]<-4*mode**0.5)
        for j in range(len(neg[0])):
            x=neg[1][j]
            y=neg[0][j]
            masks[i,y,x] = 1.0

        modes.append(mode)
        #pyf.writeto('blfn'+science[i],bigDatas[i],clobber=True)
        #mask out the entire image if it appears too faint or too bright!
        if not opt.defaultShifts:
            if dFlux[i]<0.5 or dFlux[i]>2.0:
                print 'Trimming sidereal due to flux'
                masks[i,:,:]=1.0

    w=num.where(masks>1.0)
    masks[w]=1.0


    bigDatas=ma.array(bigDatas,mask=masks)
    
 
    if bigDatas.shape[0]>9:
        shifted=ma.mean(bigDatas,axis=0).filled(0.0)+num.mean(modes)
    else:
        shifted=ma.median(bigDatas,axis=0).filled(0.0)+num.median(modes)
    pyf.writeto(opt.siderealfn,shifted,clobber=True,header=outHeader)
    

    if opt.sidereal: sys.exit()

    #now repeat the offsets but with the new ideal shifts for the moving target
    DXs=origDXs*1.
    DYs=origDYs*1.
    DXs=num.array(DXs)-num.array(OX)
    DYs=num.array(DYs)-num.array(OY)
    DXs-=rx
    DYs-=ry
    DXs-=num.min(DXs)
    DYs-=num.min(DYs)
    DXs=(DXs+0.5).astype(int)
    DYs=(DYs+0.5).astype(int)
    maxX=num.max(DXs)
    maxY=num.max(DYs)

    (A,B)=datas[0].shape
    bigDatas=num.zeros((len(datas),A+(maxY),B+(maxX))).astype('float32')
    masks=num.ones((len(datas),A+(maxY),B+(maxX))).astype('float32')
    cosmicMasks=num.zeros((len(datas),A,B)).astype('float32')
    modes=[]
    for i in range(len(DXs)):
        bgf=bgFinder.bgFinder(datas[i])
        mode=bgf.smartBackground()
        #mode=jjkmode.stats(datas[i])[0]
        #print mode

        #generate masks for cosmic ray strikes which saturate (produce negative pixels)
        neg = num.where((datas[i] < 0) | (datas[i] < mode - 3 * mode ** 0.5))
        for j in range(len(neg[0])):
            x=neg[1][j]
            y=neg[0][j]
            cosmicMasks[i:min(len(datas)-1,i+3),max(0,y-1):min(A,y+2),max(0,x-1):min(B,x+2)]=1 #do 3 images to get rid of persistence of cosmic rays

        bigDatas[i,DYs[i]:DYs[i]+A, DXs[i]:DXs[i]+B]=datas[i]-mode
        masks[i,DYs[i]:DYs[i]+A, DXs[i]:DXs[i]+B]=0.0
        masks[i,DYs[i]:DYs[i]+A, DXs[i]:DXs[i]+B]+=1.-bpMask
        masks[i,DYs[i]:DYs[i]+A, DXs[i]:DXs[i]+B]-=cosmicMasks[i]
        #pyf.writeto('bfln'+science[i],bigDatas[i],clobber=True)

        neg = num.where( bigDatas[i]<-4*mode**0.5)
        for j in range(len(neg[0])):
            x=neg[1][j]
            y=neg[0][j]
            masks[i,y,x] = 1.0

        modes.append(mode)

    w=num.where(masks>1)
    masks[w]=1

    bigDatas=ma.array(bigDatas,mask=masks)

  
    
    if bigDatas.shape[0]>9:
        shifted=ma.mean(bigDatas,axis=0).filled(0.0)+num.mean(modes)
    else:
        shifted=ma.median(bigDatas,axis=0).filled(0.0)+num.median(modes)
    pyf.writeto(opt.nsiderealfn,shifted,clobber=True,header=outHeader)
    

