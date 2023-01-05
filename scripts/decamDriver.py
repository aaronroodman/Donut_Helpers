#! /usr/bin/env python
#
# script to do all Donut stuff - for DECam images
# this script is the one that SETS ALL DEFAULTS!!!!
#
import os
import glob
from distutils.dir_util import mkpath
import shutil
import pdb
import time
from donutlib.decamutil import decaminfo
from scriptUtil import parsePipelineArguments

# collect the options

extra_definition = {"iExpoNumber":["i","i","input exposure number"],
#                    "inputFileString":["s","i","input file string"],
                    "outputDirectory":["s","o","output directory"],
                    "date":["i","d","date with format 20140212"],
                    "sequence":["i","seq","image sequence number"],
                    "versionNum":["i","ver","version number"]}

options = parsePipelineArguments(extras=extra_definition)

fititerStr = {1:"first",2:"second",3:"third",4:"fourth"}

# ("File" means the actual file name, "Name" is the name w/o directory)
# inputFileList (globbed - need to extract fileName), outputDirectory
#

# expand $env for outputDirectory
options.outputDirectory = os.path.expandvars(options.outputDirectory)

# make directory if it doesn't exist
if not os.path.exists(options.outputDirectory):
    os.mkdir(options.outputDirectory)

# form extension list
# check first for command line lists, otherwise use defaults
if options.extNamesInput == "None" or options.extNamesInput == None :
    extList = []
    dinfo = decaminfo()
    infoDict = dinfo.info()
    if options.doFandA:
        for extName in infoDict:
            if infoDict[extName]["FAflag"]:
                extList.append(extName)

    if options.doScience:
        for extName in infoDict:
            if not infoDict[extName]["FAflag"]:
                extList.append(extName)

    if options.doPlus:
        for extName in infoDict:
            if infoDict[extName]["Offset"] == 1500.0:
                extList.append(extName)

    if options.doMinus:
        for extName in infoDict:
            if infoDict[extName]["Offset"] == -1500.0:
                extList.append(extName)

                
else:
    extList = options.extNamesInput

if len(extList)==0:
    print("decamDriver: HEY YOU FOOL, you need to specify either -doFandA or -doScience")
    exit

# for building file names
# build raw file name, donut file is DECam_00366777.FS3.0009.stamp.fits  corr is DECam_00366777.FS4.corr.fits
baseName = "DECam_%08d" % (options.iExpoNumber)
imageFile = "%s/desdata/%d/%s.fits.fz" % (options.rootDirectory,options.date,baseName)
imageFile = os.path.expandvars(imageFile)


# do the Snipping?
if not options.noSnip:

    # loop over extensions
    for extName in extList:

        # split into one file per extension
        if not options.noSplitMEF:
            corrFile = os.path.join(options.outputDirectory,baseName + "." + extName + ".corr.fits")
            command = "splitMEF.py -i %s -o %s -e %s " %(imageFile,corrFile,extName)
            print(command)
            if not options.nojobFlag:
                os.system(command)
        else:
            corrFile = imageFile

        # overscan correction (and split into one file per extension)
        if not options.noOverScan:
            corrFile = os.path.join(options.outputDirectory,baseName + "." + extName + ".corr.fits")

            # fill bias and flat files
            biasFile = os.path.expandvars(options.biasFile)
            biasFile = os.path.expanduser(biasFile)

            flatFile = os.path.expandvars(options.flatFile)
            flatFile = os.path.expanduser(flatFile)

            extraOptions = ""
            if biasFile!="":
                extraOptions = "-b %s " % (biasFile)
            if flatFile!="":
                extraOptions =  extraOptions + "-f %s " % (flatFile)

            command = "overscanDECam.py -i %s -o %s -e %s %s " %(imageFile,corrFile,extName,extraOptions)
            print(command)
            if not options.nojobFlag:
                os.system(command)
        else:
            corrFile = imageFile

        # apply donut-shaped filter
        if options.doFilter:
            # filter this file for sextractor 
            filterFile = os.path.join(options.outputDirectory,baseName + "." + extName + ".filter.fits")
            #donutFilterFile = "~/Astrophysics/Donuts/donut-8.7.fits"
            donutFilterFile = options.donutFilterFile
            donutFilterFile = os.path.expanduser(donutFilterFile)
            command = "filterDecamImage.py -i %s -o %s -f %s" %(corrFile,filterFile,donutFilterFile)
            print(command)
            if not options.nojobFlag:
                os.system(command)
        else:
            filterFile = corrFile


        # delete the filterFile
        #if os.path.isfile(tempFile):
        #    os.remove(tempFile)

        # run sextractor
        catBase = "cat.%s.%s" % (baseName,extName)  # no .txt 
        catName = catBase + ".txt"
        catFileBase = os.path.join(options.outputDirectory,catBase)
        catFile = os.path.join(options.outputDirectory,catName)

        #if options.doFilter:
        #    sexConfig = "~/Astrophysics/Sextractor/decamFilter.conf"
        #else:
        #    sexConfig = "~/Astrophysics/Sextractor/decamNoFilter.conf"

        sexConfig = options.sexConfig
        sexConfig = os.path.expanduser(sexConfig)

        command = "sex %s -c %s -CATALOG_NAME %s  -CATALOG_TYPE ASCII_HEAD" % (filterFile,sexConfig,catFile)
        print(command)
        if not options.nojobFlag:
            os.system(command)

        time.sleep(1)

        if (options.doRegion):
            command = "sextoregionDECam.py -c %s -e %s -cut '%s'" % (catFile,extName,options.cutString)
            print(command)
            if not options.nojobFlag:                    
                os.system(command)

        # cut out postage stamps
        command = "sexcatServerDECam.py -c %s -f %s -od %s -o %s -cut '%s' -e %s -nP %d" % (catFileBase,corrFile,options.outputDirectory,baseName,options.cutString,extName,options.nPixels)

        print(command)
        if not options.nojobFlag:
            os.system(command)

        # delete the corrFile
        #if os.path.isfile(corrFile):
        # save for a while                os.remove(corrFile)

        # compress the corrFile
        command = "fpack -D -Y %s" % (corrFile)
        print(command)
        if not options.nojobFlag:
            os.system(command)


# do the Fitting?

# put fit results in a directory v1 (or v2, etc..)
fitoutDir = options.outputDirectory + "/v%d" % (options.versionNum)

if not options.noFit:

    # make directory if it doesn't exist
    if not os.path.exists(fitoutDir):
        os.mkdir(fitoutDir)

    # setup scratch area if desired
    if options.scratchDir!="":
        mkpath(os.path.expandvars(options.scratchDir))

    # fit em all for this file, submit separate jobs if desired!
    if not options.subJobByExt:
        imageFileString = os.path.join(options.outputDirectory,baseName + ".*.stamp.fits") 
        command = "donutfitWorkerDECam.py -i '%s' -o '%s' -d %d -seq %d -ver %d -c %s %s " % (imageFileString,fitoutDir,options.date,options.sequence,options.versionNum,options.configFile,options.remaining)
        print(command)
        if not options.nojobFlag:
            os.system(command)
    else:

        # make directory for log files...
        logdir = "%s/Donuts/logfiles/%d" % (options.rootDirectory,options.date)
        if not os.path.exists(logdir):
            os.mkdir(logdir)

        # loop over extensions, one batch job per extension
        # note the tricky \"\'%s\'\" needed to prevent bsub from expanding wildcards
        for extName in extList:
            imageFileString = os.path.join(options.outputDirectory,baseName) + ".%s.*.stamp.fits" % (extName)
            logFile = "%s/Donuts/logfiles/%d/%s_%s.log" % (options.rootDirectory,options.date,baseName,extName)
            command = "bsub %s -o %s donutfitWorkerDECam.py -i \"\'%s\'\" -o '%s' -d %d -seq %d -ver %d -c %s %s " % (options.queue,logFile,imageFileString,fitoutDir,options.date,options.sequence,options.versionNum,options.configFile,options.remaining)
            print(command)
            if not options.nojobFlag:
                os.system(command)


# do the Analysis?
if not options.noAna:

    # run on all iterations (ie. first,second etc...)
    for i in range(1,options.nFits+1):

        command = "runImageDonutAna.py -d %d -seq %d -exp %d -ver %d -iter %d -c %s %s " % (options.date,options.sequence,options.iExpoNumber,options.versionNum,i,options.configFile,options.remaining)
        print(command)
        if not options.nojobFlag:
            os.system(command)

    # all done - clean up
    # clean up scratch area...
    if options.scratchDir!="":  
        scratchDirExp = os.path.expandvars(options.scratchDir)     
        if not os.path.exists(scratchDirExp):      
            shutil.rmtree(scratchDirExp)
