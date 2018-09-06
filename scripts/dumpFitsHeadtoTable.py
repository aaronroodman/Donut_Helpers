#! /usr/bin/env python
import numpy as np
import pandas as pd
from astropy.io import fits as pyfits
import argparse
import glob
import os
import string
from scriptUtil import decodeNumberList


parser = argparse.ArgumentParser(prog='dumpFitsHead')
parser.add_argument("-i", "--inputFileName",
                  dest="inputFileName",
                  default=None,
                  help="input file names, must be in quotes and no tilde allowed")
parser.add_argument("-o", "--outputFile",
                  dest="outputFile",
                  default="temp",
                  help="output root file name")
parser.add_argument("-n", "--imageListString", dest="imageListString",default=None,
                  help="string with image lists, eg. 1:10,20,30:100 ")


# collect the options 
options = parser.parse_args()

# method to fill header tags here
def getTagList(fileName):

    # remove these tags
    headerTagsToDrop = "SIMPLE","BITPIX","NAXIS","NAXIS1","NAXIS2","CTYPE1","CTYPE2","OBJECT","DETSEC","CHECKSUM","DATASUM","TELFOCUS","DATE-OBS","TIME-OBS","DTACQNAM","RECNO","DIMMSEE","TELRA","TELDEC","UPTRTEMP","LWTRTEMP"

    # open fits file
    hdulist = pyfits.open(fileName)
    header = hdulist[0].header
    keylist = list(header.keys())

    # remove some tags
    for dropit in headerTagsToDrop:
        if keylist.__contains__(dropit):
            keylist.remove(dropit)

    hdulist.close()
    return keylist,header

# get all files using glob  --  tilde not allowed!
inputFiles = sorted(glob.glob(os.path.expandvars(options.inputFileName)))

# if imageList is set, use only inputFiles which contain these images...
if options.imageListString!=None:
    useInputFiles = []
    imageList = decodeNumberList(options.imageListString)
    imageStr = []
    for i in imageList:
        imageStr.append(str(i))

    for iFile in inputFiles:
        for iImg in imageStr:
            if iFile.find(iImg)>0:
                useInputFiles.append(iFile)
                break
    inputFiles = useInputFiles


# fill header tag dictionary, use structure of first file in list
aFile = inputFiles[0]
headerTags,header = getTagList(aFile)

# build DataFrame
df = pd.DataFrame(columns=headerTags)

# loop over files
for iFile,aFile in enumerate(inputFiles):

    print(iFile,aFile)
    hdulist = pyfits.open(aFile)
    header = hdulist[0].header

    headList = []
    for aTag in headerTags:
        headList.append(header[aTag])

    df.loc[iFile] = headList

df.to_pickle(options.outputFile+'.pkl')
    




