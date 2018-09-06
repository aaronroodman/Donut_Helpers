#! /usr/bin/env python
#
# script to do all Donut stuff - for DESI images
# this script is the one that SETS ALL DEFAULTS!!!!
#
import os
import glob
from distutils.dir_util import mkpath
import shutil
import pdb
import time
#from donutlib.decamutil import decaminfo
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
    extList = ["GFA0G"]                
else:
    extList = options.extNamesInput

# for building file names
# build raw file name, donut file is PROTODESI_GFA0G_xxxxxxxx.stamp.fits
baseName = "PROTODESI_%s_%08d" % (extList[0],options.iExpoNumber)
imageFile = "%s/%s_0000.fits" % (options.rootDirectory,baseName)
imageFile = os.path.expandvars(imageFile)


# do the Snipping?
if not options.noSnip:

    # loop over extensions
    for extName in extList:

        # run sextractor
        catBase = "cat.%s" % (baseName)  # no .txt 
        catName = catBase + ".txt"
        catFileBase = os.path.join(options.outputDirectory,catBase)
        catFile = os.path.join(options.outputDirectory,catName)

        sexConfig = options.sexConfig
        sexConfig = os.path.expanduser(sexConfig)

        command = "sex %s -c %s -CATALOG_NAME %s  -CATALOG_TYPE ASCII_HEAD" % (imageFile,sexConfig,catFile)
        print(command)
        if not options.nojobFlag:
            os.system(command)

        time.sleep(1)

        if (options.doRegion):
            command = "sextoregionDESI.py -c %s -e %s -cut '%s'" % (catFile,extName,options.cutString)
            print(command)
            if not options.nojobFlag:                    
                os.system(command)

        # cut out postage stamps
        command = "sexcatServerDESI.py -c %s -f %s -od %s -o %s -cut '%s' -e %s -nP %d" % (catFileBase,imageFile,options.outputDirectory,baseName,options.cutString,extName,options.nPixels)

        print(command)
        if not options.nojobFlag:
            os.system(command)

        # compress the corrFile
#        command = "fpack -D -Y %s" % (corrFile)
#        print command
#        if not options.nojobFlag:
#            os.system(command)


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

    # fit em all for this file,
    imageFileString = os.path.join(options.outputDirectory,baseName + "*.fits") 
    command = "donutfitWorkerDESI.py -i '%s' -o '%s' -d %d -seq %d -ver %d -c %s %s " % (imageFileString,fitoutDir,options.date,options.sequence,options.versionNum,options.configFile,options.remaining)
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