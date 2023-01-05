#! /usr/bin/env python
# 
# $Rev:: 190                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2014-09-03 10:58:53 -0700 (Wed, 03 Sep 2014)     $: 
#
# convert sextractor ASCII_HEAD file to a DS9 region file
# assume we are working on an extracted single extention file

import numpy
import argparse
import os
from astropy.io import ascii
from donutlib.decamutil import decaminfo

parser = argparse.ArgumentParser(prog='sextoregionDECam')
parser.add_argument("-c", "--catFile", dest="catFile",
                  help="input catalog file")
parser.add_argument("-cut", "--cutString",
                  dest="cutString",default="",
                  help="python code for cuts")
parser.add_argument("-e", "--extName", dest="extName",
                  help="extension Name")

options = parser.parse_args()


def checkCuts(mysexcat,irow,cutString):

    passCuts = False

    # unpack variables into scalars
    flux_auto = mysexcat["FLUX_AUTO"][irow]
    a_image = mysexcat["A_IMAGE"][irow]
    b_image = mysexcat["B_IMAGE"][irow]
    ellipticity = mysexcat["ELLIPTICITY"][irow]
    ave_image = (a_image+b_image)/2.
    flags = mysexcat["FLAGS"][irow]
    iX = mysexcat["X_IMAGE"][irow] 
    iY = mysexcat["Y_IMAGE"][irow]

    # set this to true, so region file still fills all donuts
    goodRegion = True

    # calculate rdecam
    xdecam,ydecam = dinfo.getPosition(options.extName,iX,iY)
    rdecam = numpy.sqrt(xdecam*xdecam + ydecam*ydecam)
    
    # apply cut
    if cutString != "" :
        passCuts = eval(cutString)
    else:
        passCuts = True
    return passCuts


# get decam info
dinfo = decaminfo()
infoDict = dinfo.info()


# cat file name
mysexcat = ascii.read(options.catFile)


# open region file
catBase = options.catFile[:-4]
regionName = catBase + '.reg'
regf = open(regionName,'w')
regf.write("# Region file format: DS9 version 4.1\n")
regf.write("# Filename: \n")
regf.write("global color=green dashlist=8 3 width=1 select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
regf.write("physical\n")

# Loop over Sextractor objects
for irow in range(mysexcat.as_array().shape[0]):

    # apply cuts as desired, based on cutSet
    passCuts = checkCuts(mysexcat,irow,options.cutString)
    if passCuts :

        # the -56 could be found in the Header word DATASECA or B
        # why there isn't an offset in Y too I don't know??
        x_image = mysexcat["X_IMAGE"][irow] - 56
        y_image = mysexcat["Y_IMAGE"][irow]
        a_image = mysexcat["A_IMAGE"][irow]
        b_image = mysexcat["B_IMAGE"][irow]
        theta_image = mysexcat["THETA_IMAGE"][irow]
        regf.write("ellipse(%.5f,%.5f,%.5f,%.5f,%.5f)\n" % (x_image,y_image,a_image,b_image,theta_image))
            
regf.close()




        
