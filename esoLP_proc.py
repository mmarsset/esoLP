#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 13:41:38 2016

@author: mmarsset
"""

#! /usr/bin/env python

import numpy as num
from numpy import ma
from astropy.io import fits
from scipy.ndimage import median_filter
import imp

#import pdb; pdb.set_trace()
try:
    imp.find_module('astropy')
    astropyFound = True
except ImportError:
    astropyFound = False
if astropyFound:
    from astropy.io import fits as pyf
else:
    import pyfits as pyf
    
    
def doProcStar(star,filt):
    """
    process the data.
    """
    
    print 'Processing', star
    
    hdu=fits.open('CALIB/mBias.fits')
    mBias=hdu[0].data
    hdu.close()
    
    hdu=fits.open('CALIB/m'+filt+'Flat.fits')
    mFlat=hdu[0].data
    hdu.close()

    #hdu=fits.open('CALIB/bpMap.fits')
    #bpMap=hdu[0].datan
    #hdu.close()
    
    hdu=fits.open(star)
    head=hdu[0].header
    data=hdu[0].data
    hdu.close()
        
    data-=mBias
    data=num.array(data*1.)
    data/=mFlat

    #badPixels=num.where(bpMap==0)
    #data[badPixels]=num.nan

    #cosmic ray rejection?
    if not opt.skipCosmics:
        data = removeCosmicRays(data,head['HIERARCH ESO DET OUT1 GAIN'],\
                                head['HIERARCH ESO DET OUT1 RON'],opt.sigclip,\
                                opt.sigfrac,opt.objlim,opt.niter)
    star=star[star.find('/')+1:]
    data=data.astype('int16')
    pyf.writeto('REDUCED/r'+star,data,clobber=True,header=head)

    
    
def removeCosmicRays(data,gain,readnoise,sigclip,sigfrac,objlim,niter):
    c = cosmics.cosmicsimage(data, gain=gain, readnoise=readnoise, sigclip=sigclip, sigfrac=sigfrac, objlim=objlim)
    c.run(maxiter = niter)
    data=c.cleanarray
    
    return (data)

    

def genMasterBias(bias):
    """
    generate the master bias.
    """
    
    print 'Generating master bias'
    datas=[]
    for ii in range(len(bias)):
        hdu=fits.open(bias[ii])
        data=hdu[0].data
        hdu.close()
        datas.append(data*1.)   
    
    datas=num.array(datas)
    mBias=num.median(datas,axis=0)
    #std=num.std(datas,axis=0)
    
    hdu.writeto('CALIB/mBias.fits', clobber=True)

    return (mBias)



def genMasterFlat(flats,bpMap,filt):
    """
    generate the master flat.
    """
    
    print 'Generating master flat'
    datas=[]
    masks=[]

    hdu=fits.open('CALIB/mBias.fits')
    mBias=hdu[0].data
    hdu.close()

    bpMap[:190,:]=0
    bpMap[:,:190]=0    
    bpMap[840:,:]=0
    bpMap[:,840:]=0
    bpMap[314:714,314:714]=0
    goodPix=num.where(bpMap==1)
    #pyf.writeto('CALIB/test.fits',bpMap,clobber=True)
    # badPix=num.where(bpMap==0)
    

    for ii in range(len(flats)):
        hdu=fits.open(flats[ii])
        data=hdu[0].data
        (A,B)=data.shape
        hdu.close()
        data-=mBias
        data.astype(float)        
        data=data/num.median(data[goodPix])
        
        scamp.runSex('REDUCED/michascience.sex','%s'%(flats[ii]),options={'CATALOG_NAME':'%s.cat'%(flats[ii].split('.fits')[0])}, verbose=True)
        catalog=scamp.getCatalog('%s'%(flats[ii].split('.fits')[0]+'.cat'),paramFile='REDUCED/wes.param')
        
        mask=data*0.0
        for jj in range(len(catalog['FLUX_RADIUS'])):
            width=catalog['FLUX_RADIUS'][jj]
            x=int(catalog['XWIN_IMAGE'][jj])
            y=int(catalog['YWIN_IMAGE'][jj])
            cond=width>0 and width <10 and x>20 and x<(A-20) and y>20 and y<(B-20) \
                    and not (x>320 and y>440 and x<650 and y<510)
            if cond: mask[max(0,y-6.*width):min(A,y+6.*width),max(0,x-6.*width):min(B,x+6.*width)]=1.0

            # mask[badPix]=1.0
    
        hdu = fits.PrimaryHDU(mask)
        mName='CALIB/m'+flats[ii][flats[ii].find('/')+1:]
        pyf.writeto(mName,mask,clobber=True)

        masks.append(mask)
        datas.append(data*1.)

   # masks=num.array(masks)
    datas=num.array(datas)
   # for kk in range(len(datas)):
   #     datas[kk,:,:]/=num.median(datas[kk,:,:])
   
    # VERIFIER AVEC WES QUEL EST LE PROBLEME AVEC LES MASQUES
    maskedData=ma.array(datas,mask=masks)
    hdu = fits.PrimaryHDU(maskedData)
    
    for kk in range(len(flats)):
        
        # START HACKING for the g FILTER:
        if filt=='g': masks[kk][313:355,579:634]=0
        if filt=='g': masks[kk][541:572,293:326]=0
        # END HACKING

        mName='CALIB/md'+flats[kk][flats[kk].find('/')+1:]
        maskedSlice=maskedData[kk,:,:]
        maskedSlice=num.array(maskedSlice)
        pyf.writeto(mName,maskedSlice,clobber=True)

        mmName='CALIB/mask_md'+flats[kk][flats[kk].find('/')+1:]
        pyf.writeto(mmName,masks[kk],clobber=True)
        
    mFlat=ma.median(maskedData,axis=0).filled(-32768)
    mFlat/=num.median(mFlat)
    
    
   # blurredFlat=median_filter(mFlat, size=2)
   # difference=mFlat-blurredFlat
   # threshold=10*num.std(difference)
   # bad=num.where(difference>threshold)
   # mFlat[bad]=blurredFlat[bad]
    
   # bad=num.where(mFlat==0)
   # mFlat[bad]=num.nan
    
    hdu = fits.PrimaryHDU(mFlat)
    hdu.writeto('CALIB/m'+filt+'Flat.fits', clobber=True)
    
    return (mFlat)

    
 
def genbpMap(flats,flatEXPTIMES,longStdCut=3.5,shortStdCut=3.):
    """
    generate a bad pixel map.
    """
    
    print 'Generating bad pixel map from flats'
    datas=[]
    longFlats=[]
    shortFlats=[]

    hdu=fits.open('CALIB/mBias.fits')
    mBias=hdu[0].data
    hdu.close()
    
    for ii in range(len(flats)):
        hdu=fits.open(flats[ii])
        data=hdu[0].data
        hdu.close()
        datas.append(data)
        if flatEXPTIMES[ii]>7:
            longFlats.append(data)
        if flatEXPTIMES[ii]<1.5:
            shortFlats.append(data)

    # DEMANDER A PEDRO SI C'EST UNE BONNE IDEE D'UTILISER LE BIAIS POUR EVITER LES ETOILES
    # PRESENTES SUR LES FLATS DE COURTE DUREE
    medLong=num.median(longFlats,axis=0)
    medShort=mBias#num.median(shortFlats,axis=0)
    
    bpMap=medLong*0.0+1.
    
    w=num.where((num.abs(medLong-num.median(medLong))>longStdCut*num.median(medLong)**0.5)|(medLong<=0.0))        
    bpMap[w]=0.0

    w=num.where((num.abs(medShort-num.median(medShort))>shortStdCut*num.median(medShort)**0.5))
    bpMap[w]=0.0

    hdu = fits.PrimaryHDU(bpMap)
    hdu.writeto('CALIB/bpMap.fits', clobber=True)
    
    return (bpMap)
    
    
if __name__=='__main__':
    from trippy.scamp import scamp
    import glob,sys,os,optparse,cosmics

    usage="""
    >esoLP_proc.py target
    Provide the target to create all preprocessed NIRI images.
    """
    parser=optparse.OptionParser(usage=usage)
    parser.add_option('-l',
                      action='store',dest='flatLoc',type='str',
                      default='/Users/fraserw/git/ColOSSOS/NIRI',
                      help='Location of the bad pixel mask file, flat_bpm.fits, provided by Gemini. This is in your github directory! DEFAULT=%default.')
    parser.add_option('-d',
                      action='store',dest='date',type='str',
                      default='2015-01-17',
                      help='Night to process. DEFAULT=%default.')
    parser.add_option('--sigclip',
                      action='store',default=4.0,
                      dest='sigclip', type="float",
                      help='Detection limit for cosmic rays (sigma). DEFAULT=%default')
    parser.add_option('--sigfrac',
                      action='store',default=0.5,
                      dest='sigfrac', type="float",
                      help='Fractional detection limit for neighbouring pixels. DEFAULT=%default')
    parser.add_option('--objlim',
                      action='store',default=4.0,
                      dest='objlim', type="float",
                      help='Contrast limit between cosmic rays and underlying object. DEFAULT=%default')
    parser.add_option('--niter',
                      action='store',default=4, 
                      dest='niter', type="int",
                      help='Number of iteration for cosmic ray rejection. DEFAULT=%default')
    parser.add_option('--skipCosmics',action='store_true',dest='skipCosmics',default=False,
                      help='Skip the cosmic ray rejection step. DEFAULT=%default.')


    (opt,args)=parser.parse_args()


    if len(args)<>1:
        parser.print_help()
        sys.exit()

    target=args[0]

    os.chdir(target)
    

    calFiles=glob.glob('CALIB/E*fits')
    calFiles.sort()
    bias=[]
    flats=[]
    flatEXPTIMES=[]
    flatFILTS = []
    for i in range(len(calFiles)):
        han=pyf.open(calFiles[i])
        header=han[0].header
        han.close()
        if 'BIAS' in header['OBJECT']:
            print calFiles[i],header['OBJECT']
            bias.append(calFiles[i])
        if 'SKY,FLAT' in header['OBJECT']:
            print calFiles[i],header['OBJECT']
            flats.append(calFiles[i])
            flatEXPTIMES.append(float(header['EXPTIME']))
            flatFILTS.append(header['ESO INS FILT1 NAME'])
    
    rawFiles=glob.glob('RAW/E*fits')
    raw=[]
    for i in range(len(rawFiles)):
        han=pyf.open(rawFiles[i])
        header=han[0].header
        han.close()
        if target in header['OBJECT']:
            print rawFiles[i],header['OBJECT']
            raw.append(rawFiles[i])            
    
    mBias=genMasterBias(bias)
    bpMap=genbpMap(flats,flatEXPTIMES,longStdCut=6.,shortStdCut=6.)
    
    gFlats = []
    rFlats = []
    for ii in range(len(flats)):
        if 'g#782' in flatFILTS[ii]:
            gFlats.append(flats[ii])
        if 'r#784' in flatFILTS[ii]:
            rFlats.append(flats[ii])
    

    mgFlat=genMasterFlat(gFlats,bpMap,'g')
    mrFlat=genMasterFlat(rFlats,bpMap,'r')
    
    for ii in range(len(rawFiles)):
        han=pyf.open(rawFiles[ii])
        header=han[0].header
        filt=header['ESO INS FILT1 NAME'][0]
        doProcStar(rawFiles[ii],filt)