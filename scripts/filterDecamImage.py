#! /usr/bin/env python
#
# Filter image in desired extensions with a postage stamp
#
import numpy
import scipy
from astropy.io import fits as  pyfits
import argparse
from scipy import fftpack

parser = argparse.ArgumentParser(prog='filterDecamImage')
parser.add_argument("-i", "--inputFile",
                  dest="inputFile",
                  default=None,
                  help="input file name")
parser.add_argument("-f", "--filterFile",
                  dest="filterFile",
                  default=None,
                  help="filter file name")
parser.add_argument("-o", "--outputFile",
                  dest="outputFile",
                  default="test",
                  help="output file additional name")


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
	
    return imgMean


# get filter and FFT it
def getFilterFFT(filterFile,shape):

    # read in filter shape
    hdulistFilter = pyfits.open(filterFile)
    filterData = hdulistFilter[0].data    

    # remove the baseline and normalize to 1.0
    # useful in case this is from data...
    baseline = calcOffset(filterData)
    filterData = filterData - baseline                          
    filterDataNorm = filterData/filterData.sum()

    # need to know how big this image is
    nyraw,nxraw = filterData.shape
    
    # place in center of full size array
    filterArr = numpy.zeros(shape,dtype=numpy.float64)
    nyImage,nxImage = shape

    # assume all dimensions are EVEN
    nystart = int(( nyImage/2 - 1 ) - (nyraw/2 ))
    nyend = int(( nyImage/2 - 1 ) + (nyraw/2 ))
    nxstart = int( ( nxImage/2 - 1 ) - (nxraw/2 ) )
    nxend = int(( nxImage/2 - 1 ) + (nxraw/2 ))
    
    filterArr[nystart:nyend,nxstart:nxend] = filterDataNorm
    fftFilter = fftpack.fft2(filterArr)

    return fftFilter


# read in File
hduInput = pyfits.open(options.inputFile)

# assume we have only one extension- #1
headerExt = hduInput[1].header
dataExt = hduInput[1].data
dataCast = numpy.float64(dataExt)

# get the Filter - 
fftFilter = getFilterFFT(options.filterFile,dataExt.shape)

# filter by FFT convolution
fftDataExt = fftpack.fft2(dataCast)
prod = fftDataExt * fftFilter
convImage = fftpack.ifft2(prod)
dataNew = numpy.float32(numpy.abs(fftpack.fftshift(convImage)))

# write out main header
hduOutput = pyfits.HDUList()
primaryOutput = pyfits.PrimaryHDU()
hduOutput.append(primaryOutput)

# write out the Filtered image
dataOutputHDU = pyfits.ImageHDU(dataNew)

# fill output header from input file
headerOutput = dataOutputHDU.header
for key in headerExt:
    val = headerExt[key]
    if type(val) != pyfits.header._HeaderCommentaryCards :
        headerOutput[key] = headerExt[key]

hduOutput.append(dataOutputHDU)

# write out file
hduOutput.info()
hduOutput.writeto(options.outputFile,clobber=True)




    




