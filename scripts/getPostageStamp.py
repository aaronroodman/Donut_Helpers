#! /usr/bin/env python
#
# User interface to get a postage stamp from an image
#
import numpy
from astropy.io import fits as pyfits
import subprocess
import os
import argparse
import time
import pyds9 as ds9
import pdb
import string
from collections import OrderedDict

parser = argparse.ArgumentParser(prog='getPostageStamp')
parser.add_argument("-i", "--inputFile",
                  dest="inputFile",
                  default=None,
                  help="input file name")
parser.add_argument("-e", "--extName",type=str,
                  dest="extName",
                  default=None,
                  help="extension name")
parser.add_argument("-es", "--extString",type=str,
                  dest="extString",
                  default=None,
                  help="extension name")
parser.add_argument("-o", "--outputPrefix", dest="outputPrefix",
                    default=None,
                    help="output file name prefix")
parser.add_argument("-n", "--npixels", dest="npixels",type=int,
                    default=64,
                    help="npixels")


# collect the options
options = parser.parse_args()

# set extension name
if options.extString == None:
    options.extString = options.extName

# find screen size
screen = subprocess.getoutput("xdpyinfo | grep dimensions | awk '{print $2;}'")
screenW,screenH= screen.split('x')
screenWidth = int(screenW)
screenHeight = int(screenH)


# open ds9 window with input image
ds9window = ds9.DS9(target="Main",start=True,wait=20)  
ds9window.set('height %d' % (int(screenHeight-368)))
ds9window.set('width %d' % (int(screenHeight-368)))

# load entire image
hdu = pyfits.open(options.inputFile)

# get desired extension in ds9
if options.extName != None:
    command = "file %s[%s]" % (options.inputFile,options.extName)
else:
    command = "file %s" % (options.inputFile)
ds9window.set(command)
ds9window.set("zoom to fit")
ds9window.set("scale scope local")
ds9window.set("scale zscale")

# get the data
hdu.info()
if options.extName != None:
    dataArr = hdu[options.extName].data
else:
    dataArr = hdu[0].data
inputHeader = hdu[0].header

# stamp Window
try:
    stampWin = ds9.DS9(target="Stamp",start=True,wait=30)
except:
    stampWin = ds9.DS9(target="Stamp",start=False,wait=30)

stampWin.set("view panner no")
stampWin.set("view magnifier no")
stampWin.set("width 360")
stampWin.set("height 360")
stampWin.set('height %d' % (int((screenHeight-368)/1.5)))
stampWin.set('width %d' % (int(screenHeight-368)))

# nhalfPixels
nhalfPixels = int(options.npixels/2)

# output file name
if options.outputPrefix==None:

    if options.inputFile[-8:] == '.fits.fz':
        options.outputPrefix = options.inputFile[:-8]
    else:
        options.outputPrefix = options.inputFile[:-8]


# while loop
keepgoing = True
print("getPostageStamp:  instructions: in extension window: move mouse to donut and press d to Select Donut ")
print("                                in same window press y to Write Out Donut, n to Discard,")
print("                                press q to Quit")

count = 0

while keepgoing:

    reply = ds9window.get("imexam key coordinate image")

    words = reply.split()

    if words[0] == "d":
        xval = int(float(words[1]))
        yval = int(float(words[2]))
        print((xval,yval,words))

        # now display this stamp in a separate window
        stampArr = dataArr[yval-nhalfPixels:yval+nhalfPixels,xval-nhalfPixels:xval+nhalfPixels].copy()
        istampArr = numpy.array(stampArr,dtype=int)
        stampWin.set_np2arr(istampArr)
        stampWin.set("zoom to fit")
        stampWin.set("scale minmax")

        reply = ds9window.get("imexam key data")
        words = reply.split()
        if words[0]=="y":
            ok = True
        else:
            ok = False

        # write to a file
        if ok:

            # make file name
            count = count + 1
##  old naming            outputPostfix = ".%s.stamp_x%dy%d.fits" % (options.extString,xval,yval)
            outputPostfix = ".%s.%d.stamp.fits" % (options.extString,count)
            outputName = options.outputPrefix + outputPostfix

            # dump postage stamp to new file
            stamphdu = pyfits.PrimaryHDU(stampArr)
            stamphdulist = pyfits.HDUList([stamphdu])

            # header words from input file's main header for the stamp header
            headerOutWords = ["G-SEEING","TELFOCUS","FILTER","DATE-OBS","TIME-OBS","OBSID","DTACQNAM","RECNO","DIMMSEE","TELRA","TELDEC","ZD","HA","LST","UPTRTEMP","LWTRTEMP","EXPTIME"]

            # make output header
            stamphdr = stamphdu.header

            # loop over words
            for key in headerOutWords:
                if key in inputHeader:
                    stamphdr[key] = inputHeader[key]

            # if TELFOCUS is present parse it
            if "TELFOCUS" in inputHeader:
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
                    print("getPostageStamp:  error in TELFOCUS")
                    dx = 0.
                    dy = 0.
                    dz = 0.
                    tx = 0.
                    ty = 0.
                    tz = 0.
        
            # create extraheader
            outputHeader = OrderedDict()
            outputHeader["EXTNAME"] = options.extString
            outputHeader["IX"] = xval
            outputHeader["IY"] = yval
            if "TELFOCUS" in inputHeader:
                outputHeader["DX"] = dx
                outputHeader["DY"] = dy
                outputHeader["DZ"] = dz
                outputHeader["TX"] = tx
                outputHeader["TY"] = ty
                outputHeader["TZ"] = tz

            # if DTACQNAM
            if "DTACQNAM" in inputHeader:
                ifile = int(inputHeader["DTACQNAM"][-16:-8])
                outputHeader["IFILE"] = ifile

            
            # fill some header words of the postage stamp
            for key in outputHeader:
                stamphdr[key] = outputHeader[key]                    
                       
            # write out postage stamp
            stamphdulist.writeto(outputName,clobber=True)


            
    
    elif words[0] == "q":        
        keepgoing = False

# done

ds9window.set("quit")
stampWin.set("quit")
