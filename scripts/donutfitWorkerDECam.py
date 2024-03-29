#! /usr/bin/env python
#
# Worker code to fit multiple donuts with donutfit
#
import numpy as np
import glob
import os
import shutil
from distutils.dir_util import mkpath
from donutlib.wavefit import wavefit
from donutlib.donutfit import donutfit
from donutlib.decamutil import decaminfo
from donutlib.wavefrontmap import wavefrontmap
from scriptUtil import parsePipelineArguments,dequote
import pdb

def main():

    # get options
    extra_definition = {"inputFileString":["s","i","input file string"],
                        "outputDirectory":["s","o","output directory"]}

    options = parsePipelineArguments(extras=extra_definition)

    # get all the files to fit - have to dequote here since glob won't work on quoted strings
    inputFileString = dequote(options.inputFileString)
    inputFileString = inputFileString.replace("STAR","*")
    inputFileList = glob.glob(os.path.expandvars(inputFileString))

    # expand $env for outputDirectory
    options.outputDirectory = os.path.expandvars(options.outputDirectory)

    # if we are using a scratch disk, copy all the files to the scratch area, and remake the inputFileList
    if options.scratchDir!="" :

        # setup scratch area if necessary
        scratchDirExp = os.path.expandvars(options.scratchDir)
        if not os.path.exists(scratchDirExp):
            mkpath(scratchDirExp)

        # copy stamp files over to scratch
        for i in inputFileList:
            shutil.copy2(i,"%s/." % (scratchDirExp))
        scratchFileString = "%s/*.stamp.fits" % (scratchDirExp)
        inputFileList = glob.glob(os.path.expandvars(scratchFileString))
        outputDir = scratchDirExp
    else:
        outputDir = options.outputDirectory
            
    # decode fixedParamArray strings
    fpaList1 = options.fixedParamArray1
    fpaList2 = options.fixedParamArray2
    fpaList3 = options.fixedParamArray3

    # if we want to use the wavefrontMap, creat the wavefrontMapObject now
    if options.wavefrontMapFile != None:
        wavefrontMapObject = wavefrontmap(file=options.wavefrontMapFile)
    else:
        wavefrontMapObject = None

    # initialize donutfit 
    initDict = {"nZernikeTerms":options.nZernikeTerms,"fixedParamArray1":fpaList1,"fixedParamArray2":fpaList2,"fixedParamArray3":fpaList3,"nFits":options.nFits,"nPixels":options.nPixels,"nbin":options.nbin,"scaleFactor":options.scaleFactor,"pixelOverSample":options.pixelOverSample,"iTelescope":options.iTelescope,"inputrzero":options.rzero,"debugFlag":options.debugFlag,"gain":options.gain,"wavefrontMap":wavefrontMapObject,"deltaWFM":options.deltaWFM,"waveLength":options.waveLength}
    df = donutfit(**initDict)

    # build inputZernikeDict 
    inputZernikeDict = options.inputZernikeDict

    # if options.zern4 isn't equal to 0.0, then modify the inputZernikeDict...
    if options.zern4 != 0.0 :
        for key in inputZernikeDict:
            inputZernikeDict[key][2] = options.zern4

    # form extension list for Science sensors
    extList = []
    dinfo = decaminfo()
    infoDict = dinfo.info()
    for extName in infoDict:
        if not infoDict[extName]["FAflag"]:
            extList.append(extName)
            # add to inputZernikeDict
            inputZernikeDict[extName] = [0.0,0.0,options.zern4,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
            
    # now actually do the fits over all the files
    for imageFile in inputFileList:

        # parse to get the directory and filename
        inputDirectory,inputName = os.path.split(imageFile)

        # strip the .image.fits from inputName to get the baseName
        baseName = inputName[:-11]

        # output file name prefix
        outputName = os.path.join(outputDir,baseName)
        
        fitDict  = {}
        fitDict["inputFile"] = imageFile
        fitDict["outputPrefix"] = outputName
        fitDict["inputrzero"] = options.rzero
        fitDict["inputZernikeDict"] = inputZernikeDict

        # do the fit
        df.setupFit(**fitDict)

        # do additional wavefront mesh fit
        if options.doWave:
            df.gFitFunc.closeFits()

            # now fit an extra component of the wavefront, described by a mesh of points
            winputDict = {"outputPrefix":outputName,"maxIterations":options.maxIterations,"tolerance":options.tolerance,"defineGrid":options.defineGrid,"spacing":options.spacing}

            wfit = wavefit(df,**winputDict)

            # setup initial values to 0 (redundant with defineGrid)
            values = np.zeros(wfit.npar)
            for ipar in range(wfit.npar):
                ix,iy = wfit.coarsegrid[ipar]  # in case we want starting values to vary with x,y
                values[ipar] = 0.0
            wfit.setupCoarseGrid(values)

            # do the fit
            wfit.doFit()

            # get the results
            wfit.outFit()

    # if we are using a scratch disk, copy all the files from the scratch area back to the outputDirectory
    if options.scratchDir!="":        
        # compress the donut.fits files...            
        command = "fpack -D -Y %s/*.donut.fits" % (scratchDirExp)
        os.system(command)
        
        # copy to regular outputDirectory
        scratchFileString = "%s/*.donut.fits*" % (scratchDirExp)
        outputFileList = glob.glob(os.path.expandvars(scratchFileString))
        for i in outputFileList:
            shutil.copy2(i,"%s/." % (options.outputDirectory))

    else:
        # compress the donut.fits files...            
        command = "fpack -D -Y %s/*.donut.fits" % (options.outputDirectory)
        os.system(command)

if __name__ == "__main__":
    main()
