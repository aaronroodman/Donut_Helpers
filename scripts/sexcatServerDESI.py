#! /usr/bin/env python
#
# output postage stamps from a sextractor catalog to a given directory
#
import numpy
import argparse
import os
import asciidata
import pyfits
import string
from collections import OrderedDict
from donutlib.donututil import clipPostageStamp
#from donutlib.decamutil import decaminfo


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

    # calculate rdesi
    xdesi = (iX-1024.)*0.257
    ydesi = (iY-516.)*0.257
    rdesi = numpy.sqrt(xdesi*xdesi + ydesi*ydesi)
    
    # exec the Cuts
    passstring = "passCuts = %s" % (cutString)
    if cutString != "":
        exec(passstring)
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
headerOutWords = ["EXPREQ","OBSNUM","MJD-OBS"]
originalHeader = OrderedDict()

# loop over words
for key in headerOutWords:
    if list(inputHeader.keys()).count(key)>0:
        originalHeader[key] = inputHeader[key]
    
# reset counter
nStamp = 0

# open cat file name
catFileName = "%s.txt" % (options.catPrefix)
print(catFileName)
mysexcat = asciidata.open(catFileName)

# get the data array
#dataArr = hdu[options.extName].data
dataArr = hdu[0].data

# Loop over Donuts
for irow in range(mysexcat.nrows):

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
        
        xdesi = (xval-1024.)*0.257
        ydesi = (yval-516.)*0.257
        extraheader["XDESI"] = xdesi
        extraheader["YDESI"] = ydesi

        extraheader["SEX_FLUX"] = mysexcat["FLUX_AUTO"][irow]
        extraheader["SEX_FLGS"] = mysexcat["FLAGS"][irow]
        extraheader["SEX_ELPT"] = mysexcat["ELLIPTICITY"][irow]
        extraheader["SEX_AVEI"] = (mysexcat["A_IMAGE"][irow] + mysexcat["B_IMAGE"][irow])/2.
            
        # fill header words of the postage stamp
        for key in extraheader:
            stamphdr[key] = extraheader[key]                    
                       
        # write out postage stamp
        stamphdulist.writeto(fullOutputName,clobber=True)
