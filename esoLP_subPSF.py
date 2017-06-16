#! /usr/bin/env pythonw

import optparse,sys,os,pickle,glob

from trippy import pill,psf,psfStarChooser, MCMCfit
from trippy.scamp import scamp, makeParFiles

import numpy as num
import scipy as sci
import notoriousTEX as ntt


import pylab as pyl
from stsci import numdisplay
from astLib import astWCS as WCS,astImages as Im

import Ellipse as Ell
from matplotlib.patches import Circle

from weightedMeanSTD import weightedMeanSTD

import os.path

from utils import JD2Greg
import sqlcl

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
    imp.find_module('orbfit_pyephem')
    orbfitFound = True
except ImportError:
    orbfitFound=False
if orbfitFound:
    from mp_ephem.ephem import EphemerisReader
    from mp_ephem.bk_orbit import BKOrbit 
else:
    import mpcTools
import ps1_sdss
import urllib2
from astropy.visualization import interval


def extent(r1,r2,n):
    lr1=num.log10(r1)
    lr2=num.log10(r2)

    return 10.0**(num.linspace(lr1,lr2,n))




def catalogTrim(catalog,snrCut,aperInd=False,trimDistance=35):
    global opt
    X=num.array(catalog['XWIN_IMAGE'])
    Y=num.array(catalog['YWIN_IMAGE'])
    AI=num.array(catalog['AWIN_IMAGE'])
    BI=num.array(catalog['BWIN_IMAGE'])
    F=num.array(catalog['FLUX_AUTO'])
    FE=num.array(catalog['FLUXERR_AUTO'])
    ABI=AI/BI
    if aperInd:
        AUTO_APER=num.array(catalog['MAG_AUTO'])-num.array(catalog['MAG_APER'])[aperInd]
    else:
        AUTO_APER=num.array(catalog['MAG_AUTO'])-num.array(catalog['MAG_APER'])

    w=num.where((BI>1)&(ABI<opt.maxABr)&(F/FE>snrCut)) #throw away the obvious ones
    medBI=num.median(BI[w])
    stdBI=num.std(BI[w])
    medABI=num.median(ABI[w])
    stdABI=num.std(ABI[w])
    medAA=num.median(AUTO_APER[w])
    stdAA=num.std(AUTO_APER[w])

    w=num.where((BI>medBI-stdBI)&(BI<medBI+stdBI)&(medABI<ABI+stdABI)&(AUTO_APER>medAA-stdAA)&(num.abs(AUTO_APER-medAA)<stdAA)&(F/FE>snrCut)) #only round, good sources
    #W=num.where( (((X>850)&(X<2025)) | ((X>2100)&(X<4110)) | ((X>4190)&(X<5260))) & (Y>45)&(Y<4550) )
    W=num.where( (X>45)&(X<1000) & (Y>45)&(Y<1000) )

    fw=num.intersect1d(w[0],W[0])

    """
    #now trim on distance
    gfw=[]
    for i in fw:
        d=((catalog['XWIN_IMAGE'][i]-catalog['XWIN_IMAGE'])**2+(catalog['YWIN_IMAGE'][i]-catalog['YWIN_IMAGE'])**2)**0.5
        d=num.sort(d)
        if d[1]<trimDistance:
            gfw.append(i)
    gfw=num.array(gfw)
    """
    newcatalog={}
    for i in catalog:
        newcatalog[i]=catalog[i][fw]

    #print num.min(newcatalog['YWIN_IMAGE'])

    return newcatalog


def statCut(x):
    fullSTD=num.std(x)
    fullMean=num.mean(x)
    args=num.argsort(num.abs(x-fullMean))
    X=x[args]
    for ii in range(1,int(len(x)*0.5)+1):
        std=num.std(X[:len(x)-ii])
        rat=fullSTD/std

        testRat=test(len(x)-ii,ii)
        print ii,fullSTD,std,rat,testRat
        if testRat>rat: break
    return ii-1


def test(sampSize=28,cut=1):
    a=sci.randn(sampSize)
    delt=num.abs(a-num.mean(a))
    args=num.argsort(delt)
    b=a[args][:len(a)-cut]
    return (num.std(a)/num.std(b))

def sigTrim(y,x,e,g,r,z,cut=2.):
    #print num.mean(x),num.std(x),'**'
    w=num.where(num.abs(x-num.mean(x))<cut*num.std(x))
    return (y[w],x[w],e[w], g[w], r[w], z[w])

def getSDSS(wcs):
    (minRA,maxRA,minDEC,maxDEC)=wcs.getImageMinMaxWCSCoords()
    query = querySDSS(minRA-0.1, maxRA+0.1, minDEC-0.1, maxDEC+0.1)
    xy = wcs.wcs2pix(query[:,0], query[:,1])
    out = []
    for ii in range(len(query)):
        out.append([xy[ii][0],xy[ii][1],
                    query[ii][2], query[ii][3], query[ii][4], query[ii][5], query[ii][6],
                    query[ii][0], query[ii][1],
                    query[ii][7], query[ii][8], query[ii][9], query[ii][10], query[ii][11]])
    return num.array(out)
    
def querySDSS(minRA, maxRA, minDec, maxDec):
    query = 'http://skyserver.sdss.org/dr12/SkyserverWS/SearchTools/SqlSearch?cmd=select ra,dec,u,g,r,i,z,err_u,err_g,err_r,err_i,err_z from star where ra%3E{} and ra%3C{} and dec%3E{} and dec%3C{} and err_r%3C0.25 and err_g%3C0.25&format=csv'.format(
        minRA, maxRA, minDec, maxDec).replace(' ', '%20').replace('>', '%3E').replace('<', '%EC')
    han = urllib2.urlopen(query)
    data = han.read().split('\n')
    han.close()

    sdss = []
    for ii in range(2, len(data)):
        s = data[ii].split(',')
        if len(s) < 12: break
        (ra, dec, u, g, r, i, z, du, dg, dr, di, dz) = (
        float(s[0]), float(s[1]), float(s[2]), float(s[3]), float(s[4]), float(s[5]), float(s[6]), float(s[7]),
        float(s[8]), float(s[9]), float(s[10]), float(s[11]))
        sdss.append([ra, dec, u, g, r, i, z, du, dg, dr, di, dz])
    return num.array(sdss)






def getClosestSource(catalog,cx,cy):
    #hack for slight speed improvement
    w=num.where((num.abs(catalog['XWIN_IMAGE']-cx)<40)&(num.abs(catalog['YWIN_IMAGE']-cy)<40))
    if len(w[0])>0:
        arg=num.argmin((catalog['XWIN_IMAGE'][w]-cx)**2+(catalog['YWIN_IMAGE'][w]-cy)**2)
        return ([catalog['XWIN_IMAGE'][w][arg],
                 catalog['YWIN_IMAGE'][w][arg],
                 catalog['MAG_AUTO'][w][arg],
                 catalog['FLUX_AUTO'][w][arg],
                 catalog['FLUX_AUTO'][w][arg]/catalog['FLUXERR_AUTO'][w][arg],
                 catalog['MAG_APER'][w][arg],
                 catalog['ALPHA_SKY'][w][arg],catalog['DELTA_SKY'][w][arg],
                 catalog['AWIN_IMAGE'][w][arg],catalog['BWIN_IMAGE'][w][arg],catalog['THETA_IMAGE'][w][arg]],
                 catalog['FLAGS'][w][arg])
    else:
        arg=num.argmin((catalog['XWIN_IMAGE']-cx)**2+(catalog['YWIN_IMAGE']-cy)**2)
        return ([catalog['XWIN_IMAGE'][arg],
                 catalog['YWIN_IMAGE'][arg],
                 catalog['MAG_AUTO'][arg],
                 catalog['FLUX_AUTO'][arg],
                 catalog['FLUX_AUTO'][arg]/catalog['FLUXERR_AUTO'][arg],
                 catalog['MAG_APER'][arg],
                 catalog['ALPHA_SKY'][arg],catalog['DELTA_SKY'][arg],
                 catalog['AWIN_IMAGE'][arg],catalog['BWIN_IMAGE'][arg],catalog['THETA_IMAGE'][arg]],
                 catalog['FLAGS'][arg])


def getScaledImage(sec,A=501,B=501,nsamp=200,contrast=0.5,lims=[]):
    if lims==[]:
        lims=[0,A,0,B]
    im=num.zeros([A,B]).astype(float)
    (z1,z2)=numdisplay.zscale.zscale(sec,nsamples=nsamp,contrast=contrast)
    norm=Im.normalise(sec,[z1,z2])
    #(a,b)=sec.shape
    #im[:a,:b]=norm
    im[lims[0]:lims[1],lims[2]:lims[3]]=norm
    return im



def KillCal(event):
    global line1Artist,mcsp,sdssZmagsTrim,dzpTrim,edzpTrim,sdssRAsTrim,sdssDecsTrim,mm,ss


    if event.button==1:
        px,py=event.xdata,event.ydata

        dd=((px-sdssZmagsTrim)**2+(dzpTrim-py)**2)**0.5
        arg=num.argmin(dd)
        kills[arg]=1

        w=num.where(kills==1)
        mcsp.errorbar(sdssZmagsTrim[w],dzpTrim[w],yerr=edzpTrim[w],linestyle='none',marker='o',linewidth=2,color='r')

        line1Artist[0].remove()

        w=num.where(kills==0)
        (mm,ss)=weightedMeanSTD(dzpTrim[w],edzpTrim[w]**-2)
        print 'New Median %.3f with standard uncertainty %.3f.'%(mm,ss)
        line1Artist=mcsp.plot([num.min(sdssZmagsTrim),num.max(sdssZmagsTrim)],[mm,mm],'k--')

        pyl.draw()
    elif event.button==3:
        px,py=event.xdata,event.ydata

        dd=((px-sdssZmagsTrim)**2+(dzpTrim-py)**2)**0.5
        arg=num.argmin(dd)
        kills[arg]=0

        w=num.where(kills==0)
        mcsp.errorbar(sdssZmagsTrim[w],dzpTrim[w],yerr=edzpTrim[w],linestyle='none',marker='o',linewidth=2,color='b')

        line1Artist[0].remove()

        w=num.where(kills==0)
        (mm,ss)=weightedMeanSTD(dzpTrim[w],edzpTrim[w]**-2)
        print 'New Median %.3f with standard uncertainty %.3f.'%(mm,ss)
        line1Artist=mcsp.plot([num.min(sdssZmagsTrim),num.max(sdssZmagsTrim)],[mm,mm],'k--')

        pyl.draw()


def KBOget(event):
    global secSize,cx,cy,selected,sdssSelected,sdssShifted,sp1,patchList,catalog,kron,sdss,redraw,dx,dy
    if event.button==1:
        rcx=event.xdata
        rcy=event.ydata
        if rcx==None or rcy==None: return

        cx=secSize['B']/2+(rcx-secSize['wx'])
        cy=secSize['A']/2+(secSize['wx']+1-rcy)


        #get the closest source
        (cat,flag)=getClosestSource(catalog,cx,cy)

        #clear all remaining patches and rest list
        for ii in range(len(patchList)):
            patchList[ii].remove()
        patchList=[]



        posx=cat[0]-secSize['B']/2.+secSize['wx']
        posy=-(cat[1]-secSize['A']/2.-secSize['wy']-1)
        theta=-cat[len(cat)-1]
        eA=cat[len(cat)-3]*kron
        eB=cat[len(cat)-2]*kron


        #print posx,posy,eA,eB,theta
        ellipse=Ell.TransformedEllipse((posx,posy), eA,eB, angle=theta, facecolor="none", edgecolor="red",zorder=10, linewidth=2, alpha=0.75, fix_x=1.)
        sp1.add_patch(ellipse)
        patchList.append(ellipse)

        circle=Circle((posx,posy),20,facecolor="none",edgecolor='red',linestyle='dashed',linewidth=2, alpha=0.75,zorder=10)
        sp1.add_patch(circle)
        patchList.append(circle)


        #draw in sdss
        for ii in range(len(sdss)):
            posx=sdss[ii][0]-secSize['B']/2.+secSize['wx']
            posy=-(sdss[ii][1]-secSize['A']/2.-secSize['wy']-1)
            patchList.append(Circle((posx,posy),5,facecolor="none",edgecolor='blue',linestyle='solid',linewidth=2, alpha=0.75,zorder=10))
            sp1.add_patch(patchList[len(patchList)-1])

        pyl.title('KBO Selected. Now right click on a catalog circle')
        pyl.draw()


        selected=True
        #if filts[i]<>'z_G0304':

        return

    if event.button==3:
        if not selected:
            print 'Select the KBO first please.'
            print 'Then right click on a well isolated, obviously linked sdss source (blue circle).'
            return

        rcx=event.xdata
        rcy=event.ydata

        Cx=secSize['B']/2+(rcx-secSize['wx'])
        Cy=secSize['A']/2+(secSize['wx']+1-rcy)
        if not sdssSelected:
            arg=num.argmin(((sdss[:,0]-Cx)**2+(sdss[:,1]-Cy)**2)**0.5)

            posx=sdss[arg][0]-secSize['B']/2.+secSize['wx']
            posy=-(sdss[arg][1]-secSize['A']/2.-secSize['wy']-1)
            patchList.append(Circle((posx,posy),5,facecolor="none",edgecolor='magenta',linestyle='solid',linewidth=2, alpha=0.75,zorder=10))
            sp1.add_patch(patchList[len(patchList)-1])
            pyl.title('Now right click on the corresponding star')
            pyl.draw()
            sdssSelected=(sdss[arg][0],sdss[arg][1])
        elif sdssSelected:
            pyl.title('Now middle click to process the positions.')
            (cat,flag)=getClosestSource(catalog,Cx,Cy)
            posx=cat[0]-secSize['B']/2.+secSize['wx']
            posy=-(cat[1]-secSize['A']/2.-secSize['wy']-1)
            theta=-cat[len(cat)-1]
            eA=cat[len(cat)-3]*kron
            eB=cat[len(cat)-2]*kron


            ellipse=Ell.TransformedEllipse((posx,posy), eA,eB, angle=theta, facecolor="none", edgecolor="magenta",zorder=10, linewidth=2, alpha=0.75, fix_x=1.)
            sp1.add_patch(ellipse)
            patchList.append(ellipse)

            pyl.draw()
            sdssShifted=(cat[0],cat[1])
        return
    if event.button==2:
        if (sdssSelected and sdssShifted) or len(sdss)==0:
            newPatchList=patchList[:2]
            for ii in range(len(patchList)):
                patchList[ii].remove()
            patchList=newPatchList[:]
            if len(sdss)>0:
                dx=+(-sdssShifted[0]+sdssSelected[0])
                dy=+(-sdssShifted[1]+sdssSelected[1])
                for ii in range(len(sdss)):
                    posx=sdss[ii][0]-secSize['B']/2.+secSize['wx']-dx
                    posy=-(sdss[ii][1]-secSize['A']/2.-secSize['wy']-1)+dy
                    patchList.append(Circle((posx,posy),5,facecolor="none",edgecolor='blue',linestyle='solid',linewidth=2, alpha=0.75,zorder=10))
                for ii in range(len(patchList)):
                    sp1.add_patch(patchList[ii])

            pyl.draw()
        return



def screenControl(event):
    global selected
    if event.key in ['q','Q'] and selected:
        pyl.close()
    elif not selected:
        print 'Click on the KBO!'
    else:
        print 'q, and Q  are all I understand'




description="""Description:
Photometry processing of reduced images. If the database does not contain an entry for the object, one will be created.
"""

usage="""
> esoLP_lc.py target [-o -d]

Main option is the object name. This must be provided
Use -o if you'd like to overwrite file names.
Use -g to get files but don't wipe directory.
Use -d if you'd like to reenter the database entry.
"""

parser=optparse.OptionParser(usage=usage,description=description)
parser.add_option('-o','--overwriteFiles',
                  dest='overwrite',default=False,action='store_true',
                  help="Over write the files in the object working directory. DEFAULT=%default")
#parser.add_option('-d','--overwriteDatabase',
#                  dest='database',default=False,action='store_true',
#                  help="Overwrite the database entry. DEFAULT=%default.")
parser.add_option('-g','--getFiles',
                  dest='getFiles',default=False,action='store_true',
                  help='(re)Download the image files from vospace. DEFAULT=%default')
parser.add_option('-s','--dontRunSex',
                  dest='runSex', action='store_false',default=True,
                  help="Don't run Sextractor on the target files. DEFAULT=%default")
parser.add_option('-t','--sexThreshold',
                  dest='threshold', action='store', type='float',
                  default=2,
                  help='The Sextreactor SNR threshold used for extraction. DEFAULT=%default')
parser.add_option('-a','--minArea',
                  dest='minArea',type='float',action='store',
                  default=2,
                  help='The Sextractor minium area used in the initial extraction. DEFAULT=%default')
parser.add_option('-k','--killSexFiles',
                  dest='killSexFiles',action='store_true',default=False,
                  help="Erase the sextactor files. DEFAULT=%default")
parser.add_option('--moffatWidth',
                  dest='moffatWidth', action='store', default=40,
                  type='int',
                  help='Width of the PSF fitting box when determining the moffat profile of the sources. DEFAULT=%default')
parser.add_option('--moffatSNR',
                  dest='moffatsnr',action='store',default=30.,
                  type='float',
                  help='the SNR of the sources from which to calculate the psf moffat profile. DEFAULT=%default')
parser.add_option('--PSFRadius',
                  dest='psfRadius',action='store',default=6,
                  type='float',
                  help='The PSF radius in FWHM to be generated from the stars. DEFAULT=%default.')
parser.add_option('--aperMulti',
                  dest='aperMulti',action='store',default=2.,#default=1.3
                  type='float',
                  help='The round aperture radius (times the FWHM) to be used in the photometry. DEFAULT=%default')
parser.add_option('--lineAperMulti',
                  dest='lineAperMulti',action='store',default=2.,#default=1.2
                  type='float',
                  help='The round aperture radius (times the FWHM) to be used in the photometry. DEFAULT=%default')
parser.add_option('--skyRadiusMulti',
                  action='store',type='float',dest='skyRadius',
                  default=4.,
                  help='The sky radius multps1_iplier. Sky value taken x*FWHM radius from the object. DEFAULT=%default.')
parser.add_option('--mode',
                  dest='bgMedian',action='store_true',default=False,
                  help='Use this to use the median for background estimation. Default is smart method.')
parser.add_option('-i',
                  dest='fileNum',action='store',default=0,
                  type='int',
                  help='Which file in the list to perform photometry on. Starts at 0. DEFAULT=%default')
parser.add_option('--maxABr',
                  dest='maxABr',action='store',default=1.5,
                  type='float',
                  help='The maximum A/B ratio when initially trimming away sources that arent stars. Increase this if there arent enough stars for a good PSF model. DEFAULT=%default')
parser.add_option('--displayWidth',
                  dest='displayWidth', action='store', type='int',
                  default=500,
                  help='The number of pixels (+- centre) to display when identifying the KBO. DEFAULT=%default.')
parser.add_option('--loadPSF',action='store_true',dest='loadPSF',default=False,
                  help='Use the saved PSF and aperture corrections instead of cteating from scratch. DEFAULT=%default')
parser.add_option('--bgCutFactor',action='store',dest='bgCutFactor',
                  default=3.5,type='float',
                  help='The sigma factor with which to cut background pixels when determining the backougrnd level. '
                       'DEAFULT=%default')
parser.add_option('--doPSFPhotometry',
                  action='store_true',dest='doPSFPhotometry',default=False,
                  help='Do PSF photmetry. DEFAULT=%default')
parser.add_option('--figHeight',
                  action='store',type='int',dest='figHeight',
                  default=15,
                  help='Figure height in inches. DEFAULT=%default')
parser.add_option('--dontUseSuperFlats',
                  action='store_false',dest='useSuperFlats',default=True,
                  help='Dont use the super flats (from Meg Schwamb) to process the images before photometry. DEFAULT is true.')
parser.add_option('--useGeminiSystem',
                  action='store_true', dest='useGeminiSystem', default=False,
                  help='Convert SDSS/PS1 magnitudes to the Gemini filter system? DEFAULT=%default')
parser.add_option('--PS1',action='store_true',dest='PS1',default=True,
                  help='Use PS1 Photometry instead of SDSS.')
parser.add_option('--list',default='/data/ColOSSOS/2014B.list',
                  action='store',type='str',dest='list',
                  help='The semester image list with full path. DEFAULT=%default')
parser.add_option('--imageList',default=None,
                  action='store',type='str',dest='imageList',
                  help='List of images where to perform PSF subtraction. DEFAULT=%default')
parser.add_option('--manualZeropoint',
                  dest='manualZeropoint',type='float',action='store',
                  default=-32768.,
                  help='Use a manual zeropoint, best for those frames outside the SDSS. DEFAULT=-32768.')
parser.add_option('--preliminaryResultsFile',
                  action='store',type='string',
                  dest='preliminaryResultsFile',default='/Users/fraserw/git/ColOSSOS/Preliminary.Results',
                  help='File from which to load the old parameters. ONLY USED WHEN RECALCULATING THE ZEROPOINT. DEFAULT=%default')
parser.add_option('--redoZeropoints',action='store_true',default=False,dest='redoZeropoint',
                  help='Redo zeropoint only? DEFAULT=%default')
parser.add_option('--masterPath', default = '/data/ColOSSOS/',
                  action = 'store', type = 'string', dest = 'masterPath',
                  help = 'The master inside which the super flats and mpc file diretories will be found. DEFAULT=%default')

(opt,args)=parser.parse_args()

if len(args)<>1:
    parser.print_help()
    sys.exit()

target=args[0].replace(' ','')


if opt.masterPath[-1]<>'/': opt.masterPath+='/'




#image display size
secSize={'wy':opt.displayWidth,'wx':opt.displayWidth}

#other globals
cx,cy=-1.,-1.
selected=False
sdssSelected=False
sdssShifted=False
redraw=False
dx,dy=0.,0.#-32768.,-32768.

#sex parameters
aperture=10#range(55,4,-3)
kron=2.
min_rad=3.5


bgmode='smart'
if opt.bgMedian:
    bgmode='median'
    print 'Using the median for background estimation.\n'



files=[]
filts=[]


os.chdir(target+'/REDUCED')


if opt.imageList<>None:
    with open(opt.imageList) as han:
        data=han.readlines()
    for i in range(len(data)):
        files.append(data[i].split()[0])
else:
    files = glob.glob('r*.fits')
    wfiles = []
    for i in range(len(files)):
        if 'psf' not in files[i] and 'PSFsub' not in files[i]: 
            wfiles.append(files[i])
    files=wfiles


for ii in range(len(files)):
    with pyf.open(files[ii]) as han:
            header=han[0].header
    filts.append(header['HIERARCH ESO INS FILT1 NAME'])

names=[target]

ntt.plotSettings(10,opt.figHeight*2./3.)
out=[]
zeros=[]


wfiles=[]
wfilts=[]
for i in range(len(files)):
    if filts[i]=='r#784': 
        wfiles.append(files[i])
        wfilts.append(filts[i])
files=wfiles
filts=wfilts


if 1:
    for i in range(len(files)):

        # load the previous entry string if we are redoing the zeropoint
        redoEntry=None
        if opt.redoZeropoint:
            with open(opt.preliminaryResultsFile) as rhan:
                ddata=rhan.readlines()
            for k in range(len(ddata)):
                if files[i] in ddata[k]:
                    redoEntry=ddata[k].split()
            for k in range(2,len(redoEntry)):
                redoEntry[k]=float(redoEntry[k])
        #end previous entry loading


        entry=[files[i],filts[i]]
        if filts[i]=='g#782': nominalZPT=28.24
        elif filts[i]=='r#784': nominalZPT=28.34
        else: nominalZPT=0.0

        if opt.manualZeropoint<>-32768.: nominalZPT=opt.manualZeropoint

        sdssZmags=[]
        sdssZerrs=[]
        dzp=[]
        edzp=[]
        sdssgs=[]
        sdssrs=[]
        sdsszs=[]
        mm=-1.
        ss=-1.


        cx=-1
        cy=-1
        dx=-32768.
        dy=-32768.
        selected=False
        sdssSelected=False
        sdssShifted=False
        redraw=False

        
        fname = files[i].split('.fits')[0]
        catName=fname+'.cat'
        print
        if opt.runSex: scamp.runSex('michascience.sex', fname+'.fits',options={'CATALOG_NAME':catName},verbose=True)
        #       scamp.runSex('science.sex','rn%s'%(science[ii]),options={'CATALOG_NAME':'rn%s.cat'%(science[ii].split('.')[0])})        
        catalog=scamp.getCatalog(catName,paramFile='def.param')
        

        trimmedCatalog=catalogTrim(catalog, opt.moffatsnr)


        

        with pyf.open(files[i]) as han:
            data=han[0].data
            header=han[0].header

        wcs=WCS.WCS(header,mode='pyfits')
        
        airMassStart = (header['HIERARCH ESO TEL AIRM START'])
        airMassEnd = (header['HIERARCH ESO TEL AIRM END'])
        AIRMASS = (airMassStart+airMassEnd)/2.


        MJD=header['MJD-OBS']+float(header['EXPTIME'])/3600./24.
        EXPTIME=float(header['EXPTIME'])

        print 'Image MJD {}'.format(MJD)

        mpcFile=None
        
        for j in range(len(names)):
            if os.path.isfile('../../../mpcs/{}.mpc'.format(names[j])):
                mpcFile='../../../mpcs/{}.mpc'.format(names[j])
                break
        if mpcFile==None:
            raise Exception("These aren't the mpc droids you are looking for.")

        if orbfitFound:
            obs = EphemerisReader().read(filename='../../../mpcs/%s.mpc'%(target))
            orb = BKOrbit(obs)
            orb.predict(MJD+2400000.5)
            ra1,dec1=orb.coordinate.ra.degree,orb.coordinate.dec.degree
            orb.predict(MJD+2400000.5+(1./24.))
            ra2,dec2=orb.coordinate.ra.degree,orb.coordinate.dec.degree
        else:
            DATE=JD2Greg(MJD+2400000.5)[0]+str(MJD-int(MJD))[1:7]
            (ra1,dec1,er,ed,ang)= mpcTools.runPredict(opt.masterPath+'../../../mpcs/%s.abg'%(target),DATE)
            DATE=JD2Greg(MJD+2400000.5)[0]+str(MJD-int(MJD)+1./24.)[1:7]#+'.541666'
            (ra2,dec2,er,ed,ang)= mpcTools.runPredict(opt.masterPath+'../../../mpcs/%s.abg'%(target),DATE)

        rate=((ra2-ra1)**2+(dec2-dec1)**2)**0.5*3600.

        (x1,y1)=wcs.wcs2pix(ra1,dec1)
        (x2,y2)=wcs.wcs2pix(ra2,dec2)
        angle=num.arctan2(y2-y1,x2-x1)*180./num.pi
        if angle>90: angle-=180.
        elif angle<-90:angle+=180.


        entry.append(MJD)
        entry.append(AIRMASS)
        entry.append(rate)
        entry.append(angle)

        print '\nObject %s moving at %s "/hr at angle %s degrees counter clockwise from horizontal.'%(files[i],rate,angle)
        print 'Object is near ',x1,x2
        print 'RA, Dec',ra1,dec1
        print

        phot=pill.pillPhot(data,repFact=10)


        starChooser=psfStarChooser.starChooser(data,trimmedCatalog['XWIN_IMAGE'],trimmedCatalog['YWIN_IMAGE'],trimmedCatalog['FLUX_AUTO'],trimmedCatalog['FLUXERR_AUTO'])
        (goodFits,goodMeds,goodSTDs)=starChooser(opt.moffatWidth,opt.moffatsnr,xWidth=51,yWidth=51)


        aa=int(2*opt.psfRadius*goodMeds[0])+1
        if aa%2==0: aa+=1
        goodPSF=psf.modelPSF(num.arange(aa),
                             num.arange(aa),
                             alpha=goodMeds[2],beta=goodMeds[3],repFact=10)
        fwhm=goodPSF.FWHM(fromMoffatProfile=True)
        print '\nProfile has FWHM of %s pix or %s arcsec.'%(fwhm,fwhm*0.0728)
        print '\n'

        print 'Generating the lookup table...'
        goodPSF.genLookupTable(data,goodFits[:,4],goodFits[:,5])
        goodPSF.genPSF()

        #print goodFits[:,4]
        #print goodFits[:,5]
        #print goodMeds[2],goodMeds[3],goodMeds[0]
        #sys.exit()
        print '   Done.\n'


        sourceAperCorrs=[]

        print 'Computing aperture corrections of each good star...'
        for j in range(len(goodFits)):


            #using the lineAperRadius to get the uncertainty on the **line** aperture correction
            #this is not a bug!
            phot.computeRoundAperCorrFromSource(goodFits[j][4],goodFits[j][5],
                                                num.array([opt.aperMulti*fwhm,4.*fwhm]),
                                                skyRadius=4.5*fwhm,width=7*fwhm,mode=bgmode,
                                                displayAperture=False,display=False)
            sourceAperCorrs.append(phot.roundAperCorr(opt.aperMulti*fwhm))
        #print sourceAperCorrs

        roundAperCorrError=num.nanstd(sourceAperCorrs)*len(sourceAperCorrs)**-0.5

        print '   done.\n'

        print
        print '\nComputing round aperture corrections'

        goodPSF.computeRoundAperCorrFromPSF(extent(min(0.8,opt.aperMulti*0.99)*fwhm,4*fwhm,10),display=False,
                                                      displayAperture=False,
                                                      useLookupTable=True)
        roundAperCorr=goodPSF.roundAperCorr(opt.aperMulti*fwhm)
        if roundAperCorrError==0.0: roundAperCorrError=0.01
        print "Using round aperture of %.4f pixel radius with aperture correction of %.4f mags and uncertainty %.4f mags."%(fwhm*opt.aperMulti,roundAperCorr,roundAperCorrError)

        print '\nComputing linear aperture corrections'
        goodPSF.line(rate,angle,EXPTIME/3600.,pixScale=0.0728,useLookupTable=True)
        goodPSF.computeLineAperCorrFromTSF(extent(min(0.8,opt.lineAperMulti*0.99)*fwhm,4*fwhm,10),
                                                     l=(EXPTIME/3600.)*rate/.0728,a=angle,display=False,displayAperture=False)
        #goodPSF.computeLineAperCorrFromMoffat(num.linspace(0.8*opt.lineAperMulti*fwhm,4*fwhm,7),l=(EXPTIME/3600.)*rate/.0728,a=angle,display=False,displayAperture=False)
        lineAperCorr=goodPSF.lineAperCorr(opt.lineAperMulti*fwhm)
        lineAperCorrError=-2.5*num.log10(1.-(10.0**(roundAperCorrError/2.5)-1)*10.0**((lineAperCorr-roundAperCorr)/2.5))


        print 'Using line aperture correction of %.4f magnitudes and uncertainty %.4f.\n'%(lineAperCorr,lineAperCorrError)

        """
        fig=pyl.figure('Aperture Corrections')
        x=num.linspace(0.5*opt.lineAperMulti,4*opt.aperMulti*fwhm,20)
        pyl.plot(x,goodPSF.roundAperCorr(x),'k-',label='Round')
        pyl.plot(x,goodPSF.lineAperCorr(x),'r-',label='Line')
        pyl.xlabel('Radius')
        pyl.ylabel('Aperture Correction (mags)')
        pyl.legend()
        pyl.show()
        pyl.close()
        """

        goodPSF.psfStore(fname+'.psf.fits')


        with open(fname+'.psfStars.pickle','w+') as psfHan:
            pickle.dump([goodFits,sourceAperCorrs,opt.aperMulti*fwhm],psfHan)



        entry.append(fwhm)
        entry.append(fwhm*opt.aperMulti)
        entry.append(fwhm*opt.lineAperMulti)
        entry.append(roundAperCorr)
        entry.append(lineAperCorr)
        entry.append(lineAperCorrError)

        #retrieve the SDSS magnitudes
        if opt.PS1:
            print '  Using PS1 photometry.'
            ps1=ps1_sdss.sdssPS1(catPath='/Users/mmarsset/data/PSCatalogs/')
            ps1.loadCat(target+'.tsv')
            psc=ps1.convert()

            sdss=[]
            for jj in range(len(psc)):
                (x,y)=wcs.wcs2pix(psc[jj][0],psc[jj][1])
                gmr_sdss=psc[jj][3]-psc[jj][4]
                #if psc[jj,ind]>0.03 and len(sdss)>20:break
                sdss.append([x,y,psc[jj][2],
                            psc[jj][3]-0.012, #hacks to get the PS1 catalog onto the SDSS catalog ***MUST BE NEGATIVE
                             psc[jj][4]-0.048,
                             psc[jj][5],
                             psc[jj][6],
                             psc[jj][0],
                             psc[jj][1],
                             psc[jj][7],
                             psc[jj][8],
                             psc[jj][9],
                             psc[jj][10],
                             psc[jj][11]])
            sdss=num.array(sdss)
        else:
            sdss=getSDSS(wcs)

        if opt.useGeminiSystem:
            for iii in range(len(sdss)):
                gmr=sdss[iii][3]-sdss[iii][4]
                sdss[iii][3]-=0.146*gmr
                sdss[iii][4]-=0.051*gmr


        fig1=pyl.figure(1)
        sp1=fig1.add_subplot(1,1,1)
        sp1.set_aspect('equal','datalim')
        sp1.autoscale_view(True,True,True)
        patchList=[]



        (A,B)=data.shape
        secSize['A']=A
        secSize['B']=B

        #scale the secsize to ensure that at least two SDSS sources are in the frame
        numSDSSPlotted=0
        while numSDSSPlotted<3:
            for k in range(len(sdss)):
                if abs(sdss[k][0] - secSize['B'] / 2. )<secSize['wx'] and abs(sdss[k][1] - secSize['A'] / 2.)<secSize['wy']:
                    numSDSSPlotted+=1
            if numSDSSPlotted<3:
                secSize['wx'] *= 1.2
                secSize['wy'] *= 1.2
                numSDSSPlotted = 0

        sec=data[A/2-secSize['wy']:A/2+secSize['wy']+1,B/2-secSize['wx']:B/2+secSize['wx']+1][::-1,:]

        print sec.shape
        print secSize
        scaled=getScaledImage(sec,A=2*secSize['wy']+1,B=2*secSize['wx']+1)
        axIm=pyl.imshow(scaled*1.,cmap='gray')

        px=secSize['wx']+(x1-B/2)
        py=secSize['wy']-(y1-A/2)

        ellipse=Ell.TransformedEllipse((px,py), 40.,40., angle=0.0, facecolor="none", edgecolor="yellow",zorder=10, linewidth=3, alpha=0.35, fix_x=1.)
        sp1.add_patch(ellipse)


        pyl.title('Mark the KBO (left click). Yellow is apporximately where the KBO should be.')
        pyl.connect('button_press_event',KBOget)
        pyl.connect('key_press_event',screenControl)

        pyl.show()
        
        (cat,flag)=getClosestSource(catalog,cx,cy)
        cx=cat[0]
        cy=cat[1]

        psfHdu=pyf.open(fname+'.psf.fits')
        psfData=psfHdu[0].data
        psfHdu.close()
        
        bg = phot.bg
        trimData=data[cy-150.:cy+150.,cx-150.:cx+150.]-bg
        
        fitter=MCMCfit.MCMCfitter(goodPSF,trimData)
        fitter.fitWithModelPSF(150+cx-int(cx)-1,150+cy-int(cy)-1, m_in=1000.,fitWidth=10, nWalkers=20, nBurn=20, nStep=20, bg=phot.bg, useLinePSF=True, verbose=True,useErrorMap=False)
        
        
        (fitPars,fitRange)=fitter.fitResults(0.67)
        print fitPars
        print fitRange
        
        
        modelImage=goodPSF.plant(fitPars[0],fitPars[1],fitPars[2],trimData,addNoise=False,useLinePSF=True,returnModel=True)
        pyl.imshow(modelImage)
        pyl.show()
        
        removedData=goodPSF.remove(fitPars[0],fitPars[1],fitPars[2],trimData,useLinePSF=True)
        
        (z1,z2)=numdisplay.zscale.zscale(removedData)
        normer=interval.ManualInterval(z1,z2)
        
        pyl.imshow(normer(trimData))
        pyl.show()
        
        pyl.imshow(normer(removedData))
        pyl.show()
    
        data[cy-150.:cy+150.,cx-150.:cx+150.]=removedData+bg
        

        pyf.writeto('t'+fname+'.fits', trimData, clobber=True)
        pyf.writeto('t'+fname+'_PSFsub.fits', removedData, clobber=True)
        
        pyf.writeto(fname+'_PSFsub.fits', data, clobber=True, header=header)