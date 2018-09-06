#! /usr/bin/env python

# 
# $Rev:: 207                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2015-08-03 11:20:13 -0700 (Mon, 03 Aug 2015)     $: 
#
#
# script to do overscan correction
#
import numpy
import scipy
import pyfits
import argparse
import pdb

parser = argparse.ArgumentParser(prog='overscanDECam')
parser.add_argument("-i", "--inputFile",
                  dest="inputFile",
                  default=None,
                  help="input file name")
parser.add_argument("-o", "--outputFile",
                  dest="outputFile",
                  default=None,
                  help="input file name")
parser.add_argument("-b", "--biasFile",
                  dest="biasFile",
                  default="",
                  help="bias file name")
parser.add_argument("-f", "--flatFile",
                  dest="flatFile",
                  default="",
                  help="flat file name")

parser.add_argument("-n", "--nOverscanCol",
                  dest="nOverscanCol",type=int,
                  default=56,
                  help="number of Overscan Columns (should get from header!!)")


parser.add_argument("-e", "--extName",
                  dest="extName",
                  default="",
                  help="extension Name")



# collect the options 
options = parser.parse_args()


# calculate offset for extension
def calcOffset(locarray):

    # get Mean,Sigma, Mask of pixels below Mean+2.5*sigma
    # use numpy Masked Array!

    # remove zeros (or low values in dead amplifier)
    maskArray = numpy.ma.masked_less(locarray,10.0)

    # cut at 1 sigma in first iteration
    imgMean = maskArray.mean()
    imgStd = maskArray.std()
    maskArray = numpy.ma.masked_greater(maskArray,imgMean+1*imgStd)
    countsNew = (maskArray.mask==False).sum()
    countsOld = -1

    while countsOld!=countsNew:
        countsOld = countsNew
        imgMean = maskArray.mean()
        imgStd = maskArray.std()
        maskArray = numpy.ma.masked_greater(maskArray,imgMean+2.5*imgStd)
        countsNew = (maskArray.mask==False).sum()
        print("calcOffset: ",imgMean,imgStd)
	
    return imgMean




# read in file
# assume only data is in header 1
hduInput = pyfits.open(options.inputFile)

# primary header
headerPrimary = hduInput[0].header

# if no extName is set, assume our data is in Ext 1
if options.extName == "":
    headerExt1 = hduInput[1].header
    dataExt1 = hduInput[1].data
else:
    headerExt1 = hduInput[options.extName].header
    dataExt1 = hduInput[options.extName].data
    

# now overscan correct the data
#
# notice that first and last 6 columns are weird -- so ignore those
# these are PreScan's and should be ignored!!!
#
nEdgeSkip = 6
nOverscanColUsed = float(options.nOverscanCol - nEdgeSkip)
nRows,nCols = dataExt1.shape
nColsHalf = nCols/2

# AmpB

# extract just the desired columns for the overscan and average them
overscanSecB = dataExt1[:,nEdgeSkip:options.nOverscanCol]

# find the median
overscanSecBMedian = numpy.median(overscanSecB,1)

# now explode this column to the full Amp size, with a python trick
overscanSecBFull = overscanSecBMedian.repeat(nColsHalf).reshape(nRows,nColsHalf)


# AmpA

# extract just the desired columns for the overscan and average them
overscanSecA = dataExt1[:,-options.nOverscanCol:-nEdgeSkip]

# find the median
overscanSecAMedian = numpy.median(overscanSecA,1)

# now explode this column to the full Amp size, with a python trick
overscanSecAFull = overscanSecAMedian.repeat(nColsHalf).reshape(nRows,nColsHalf)

# now combine AmpA and AmpB halfs
overscanFull = numpy.column_stack((overscanSecBFull,overscanSecAFull))

# subtract overscan
dataCorr = dataExt1 - overscanFull

# also do bias correction - using a master Bias frame
if options.biasFile !="":

    # read in file
    hduBias = pyfits.open(options.biasFile)

    # extName must be set
    if options.extName == "":
        print("HEY FOOL, must specifiy an extName when subtracting a Bias Frame")
        exit()
    else:
        dataBias = hduBias[options.extName].data
    
    # subtract
    dataCorr = dataCorr - dataBias

# also do flat correction - using a master Flat frame
if options.flatFile !="":

    # read in file
    hduFlat = pyfits.open(options.flatFile)

    # extName must be set
    if options.extName == "":
        print("HEY FOOL, must specifiy an extName when subtracting a Flat Frame")
        exit()
    else:
        dataFlat = hduFlat[options.extName].data
    
    # divide
    dataCorr = dataCorr / dataFlat

# write out main header
hduOutput = pyfits.HDUList()

# make primary header
primaryOutput = pyfits.PrimaryHDU()
hduOutput.append(primaryOutput)

# fill primary header from original file
primaryOutputHeader = primaryOutput.header
for key in headerPrimary:
    val = headerPrimary[key]
    if type(val) != pyfits.header._HeaderCommentaryCards :
        primaryOutputHeader[key] = val

# put data in extension #1
dataOutputHDU = pyfits.ImageHDU(dataCorr)

# fill header
headerOutput = dataOutputHDU.header
for key in headerExt1:
    val = headerExt1[key]
    if type(val) != pyfits.header._HeaderCommentaryCards :
        headerOutput[key] = headerExt1[key]

hduOutput.append(dataOutputHDU)

# write out file
hduOutput.info()
hduOutput.writeto(options.outputFile,clobber=True)



    




