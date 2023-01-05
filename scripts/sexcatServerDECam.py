#! /usr/bin/env python
#
# $Rev:: 207                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2015-08-03 11:20:13 -0700 (Mon, 03 Aug 2015)     $:  
#
# output postage stamps from a sextractor catalog to a given directory
#
import numpy as np
import argparse
import os
from astropy.io import ascii
from astropy.io import fits as pyfits
import string
from collections import OrderedDict
from donutlib.donututil import clipPostageStamp
from donutlib.decamutil import decaminfo


parser = argparse.ArgumentParser(prog='sexcatServerDECam')
parser.add_argument("-c", "--catPrefix", dest="catPrefix",
                  help="input catalog filename prefix")
parser.add_argument("-f", "--fitsFile", dest="fitsFile",
                  help="input fits file")
parser.add_argument("-od", "--outputDirectory", dest="outputDirectory",
                  help="output directory")
parser.add_argument("-o", "--outputPrefix", dest="outputPrefix",
                  help="output filename prefix")

parser.add_argument("-e", "--extName", dest="extName",
                  help="extension Name")


parser.add_argument("-i", "--fileNum",
                  dest="fileNum",type=int,default=None,
                  help="fits file number")
parser.add_argument("-cut", "--cutString",
                  dest="cutString",default="",
                  help="python code for cuts")
parser.add_argument("-nP", "--nPixels",
                  dest="nPixels",type=int,
                  default=64,
                  help="nPixels by nPixels")


# define clearRegions:
#clearRegion = {"FS1":[1081,2104],"FS2":[1081,2104],"FS3":[1081,2104],"FS4":[1081,2104],"FN1":[57,1080],"FN2":[57,1080],"FN3":[57,1080],"FN4":[57,1080]}

# add Science sensors to clear region - all good there
dinfo = decaminfo()
infoDict = dinfo.info()
#for extName in infoDict:
#    if not infoDict[extName]["FAflag"]:
#        clearRegion[extName] = [57,2104]

# method to add data to catalog
def addData(mysexcat,extName):
    """ add some useful columns 
    """
    mysexcat['AVE_IMAGE'] = (mysexcat['A_IMAGE'] + mysexcat['B_IMAGE'])/2.

    rdecam_list = []
    for irow in range(mysexcat.as_array().shape[0]):
        # calculate whether this is in the good region (use only 1/2 of the chips for now)
        iX = mysexcat["X_IMAGE"][irow] 
        iY = mysexcat["Y_IMAGE"][irow]

        # calculate rdecam
        xdecam,ydecam = dinfo.getPosition(extName,iX,iY)
        rdecam_list.append(np.sqrt(xdecam*xdecam + ydecam*ydecam))

    mysexcat['RDECAM'] = np.array(rdecam_list)

# method to check that cats are passed
def checkCuts(mysexcat,irow,extName,cutString):

    passCuts = False

    # unpack variables into scalars
    flux_auto = mysexcat["FLUX_AUTO"][irow]
    a_image = mysexcat["A_IMAGE"][irow]
    b_image = mysexcat["B_IMAGE"][irow]
    ellipticity = mysexcat["ELLIPTICITY"][irow]
    ave_image = (a_image+b_image)/2.
    flags = mysexcat["FLAGS"][irow]

    # calculate whether this is in the good region (use only 1/2 of the chips for now)
    iX = mysexcat["X_IMAGE"][irow] 
    iY = mysexcat["Y_IMAGE"][irow]

    #goodRegion = False
    #if iX>=clearRegion[extName][0] and iX<=clearRegion[extName][1]:
    #    goodRegion = True

    # calculate rdecam
    xdecam,ydecam = dinfo.getPosition(options.extName,iX,iY)
    rdecam = np.sqrt(xdecam*xdecam + ydecam*ydecam)
    
    # exec the Cuts
    if cutString != "":
        passCuts = eval(cutString) 
    else:
        passCuts = True
    return passCuts

    

# unpack options
options = parser.parse_args()

# nhalfPixel
nhalfPixels = options.nPixels/2

# check that outputDirectory exists, if not make it
if not os.path.exists(options.outputDirectory):
    os.mkdir(options.outputDirectory)

# open fits file
fitsFileName = options.fitsFile
hdu = pyfits.open(fitsFileName)
inputHeader = hdu[0].header

# get some words from original header
headerOutWords = ["G-SEEING","TELFOCUS","FILTER","DATE-OBS","TIME-OBS","OBSID","OBSTYPE","DTACQNAM","FILENAME","RECNO","DIMMSEE","TELRA","TELDEC","ZD","HA","LST","UPTRTEMP","LWTRTEMP","EXPTIME","RA","DEC","MSURTEMP","MAIRTEMP","MJD-OBS"]
originalHeader = OrderedDict()

# if TELFOCUS is present parse it
if list(inputHeader.keys()).count("TELFOCUS")>0:
    try:
        telfocusstr = inputHeader["TELFOCUS"]
        words = telfocusstr.split(",")

        dx = float(words[0])
        dy = float(words[1])
        dz = float(words[2])
        tx = float(words[3])
        ty = float(words[4])
        tz = float(words[5])
    except:
        dx = -999999.9
        dy = -999999.9
        dz = -999999.9
        tx = -999999.9
        ty = -999999.9
        tz = -999999.9
    # add to originalHeader
    originalHeader["DX"] = dx
    originalHeader["DY"] = dy
    originalHeader["DZ"] = dz
    originalHeader["TX"] = tx
    originalHeader["TY"] = ty
    originalHeader["TZ"] = tz
    
# if DTACQNAM
if list(inputHeader.keys()).count("DTACQNAM")>0:
    try:
        ifile = int(inputHeader["DTACQNAM"][-16:-8])
        originalHeader["IFILE"] = ifile
    except:
        print("sexcatServerDECam: problem unpacking DTACQNAM")

# if FILENAME
if list(inputHeader.keys()).count("FILENAME")>0:
    try:
        ifile = int(inputHeader["FILENAME"][6:14])
        originalHeader["IFILE"] = ifile
    except:
        print("sexcatServerDECam: problem unpacking FILENAME")

# loop over words
for key in headerOutWords:
    if list(inputHeader.keys()).count(key)>0:
        originalHeader[key] = inputHeader[key]

    
# reset counter
nStamp = 0

# open cat file name
catFileName = "%s.txt" % (options.catPrefix)
print(catFileName)
mysexcat = ascii.read(catFileName)

# get the data array
dataArr = hdu[options.extName].data

# Loop over Donuts
for irow in range(mysexcat.as_array().shape[0]):

    # apply cuts as desired, based on cutSet
    passCuts = checkCuts(mysexcat,irow,options.extName,options.cutString)
    if passCuts :

        nStamp = nStamp + 1
        xval = int(mysexcat["X_IMAGE"][irow])
        yval = int(mysexcat["Y_IMAGE"][irow])
        print("Found another at ",options.extName,xval,yval)

        # get the data
        stampArr = clipPostageStamp(xval,yval,dataArr,options.nPixels)

        # output to fits file in outputDirectory
        outputPostfix = ".%s.%04d.stamp.fits" % (options.extName,nStamp)
        fullOutputName = os.path.join(options.outputDirectory,options.outputPrefix+outputPostfix)

        # dump postage stamp to new file
        stamphdu = pyfits.PrimaryHDU(stampArr)
        stamphdulist = pyfits.HDUList([stamphdu])

        # make output header
        stamphdr = stamphdu.header

        # fill header: start with words from original file's header
        for key in originalHeader:
            stamphdr[key] = originalHeader[key]

        # more header variables from sextractor
        extraheader = OrderedDict()
        extraheader["EXTNAME"] = options.extName
        extraheader["IX"] = xval - 56  # output is 1 to 2048
        extraheader["IY"] = yval - 51  # output is 1 to 2048
        extraheader["X_IMAGE"] = xval
        extraheader["Y_IMAGE"] = yval
        extraheader["ISTAMP"] = nStamp
        
        xdecam,ydecam = dinfo.getPosition(options.extName,extraheader["IX"],extraheader["IY"])
        extraheader["XDECAM"] = xdecam
        extraheader["YDECAM"] = ydecam

        extraheader["SEX_FLUX"] = mysexcat["FLUX_AUTO"][irow]
        extraheader["SEX_FLGS"] = mysexcat["FLAGS"][irow]
        extraheader["SEX_ELPT"] = mysexcat["ELLIPTICITY"][irow]
        extraheader["SEX_AVEI"] = (mysexcat["A_IMAGE"][irow] + mysexcat["B_IMAGE"][irow])/2.
            
        # fill header words of the postage stamp
        for key in extraheader:
            stamphdr[key] = extraheader[key]                    
                       
        # write out postage stamp
        stamphdulist.writeto(fullOutputName,overwrite=True)
