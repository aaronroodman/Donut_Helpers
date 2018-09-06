#! /usr/bin/env python
# 
# $Rev:: 203                                                       $:  
# $Author:: roodman        $:  
# $LastChangedDate:: 20#$: 
#

import glob
import os
import shutil
import pickle
import argparse
import pyfits
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import gROOT
from donutlib.donutana import donutana
from collections import OrderedDict
from scriptUtil import parsePipelineArguments
import pdb

# unpack options
extra_definition = {"date":["i","d","date with format 20140212"],
                    "sequence":["i","seq","image sequence number"],
                    "versionNum":["i","ver","version number"],
                    "iteration":["i","iter","iteration number"],
                    "expoNumber":["i","exp","exposure number"]}

options = parsePipelineArguments(extras=extra_definition)

# turn off interactive ROOT graphics
gROOT.SetBatch()

# pick meshName by iteration number
if options.iteration==1:
    meshName = options.meshName1
elif options.iteration==2:
    meshName = options.meshName2
elif options.iteration==3:
    meshName = options.meshName3
else:
    print("runImageDonutAna:  no meshName is setup!")


# Fit to reference mesh
inDict = {"zPointsFile":options.meshDir+"/z4Mesh_"+meshName+".dat", 
          "z5PointsFile":options.meshDir+"/z5Mesh_"+meshName+".dat",
          "z6PointsFile":options.meshDir+"/z6Mesh_"+meshName+".dat",
          "z7PointsFile":options.meshDir+"/z7Mesh_"+meshName+".dat",
          "z8PointsFile":options.meshDir+"/z8Mesh_"+meshName+".dat",
          "z9PointsFile":options.meshDir+"/z9Mesh_"+meshName+".dat",
          "z10PointsFile":options.meshDir+"/z10Mesh_"+meshName+".dat",
          "z11PointsFile":options.meshDir+"/z11Mesh_"+meshName+".dat",
          "z14PointsFile":options.meshDir+"/z14Mesh_"+meshName+".dat",
          "z15PointsFile":options.meshDir+"/z15Mesh_"+meshName+".dat",
          "sensorSet":options.sensorSet,
          "doTrefoil":True,
          "doSpherical":options.doSpherical,
          "doQuadrefoil":options.doQuadrefoil,
          "doRzero":options.doRzero,
          "unVignettedOnly":False,
          "interpMethod":"idw",
          "histFlag":True,
          "donutCutString":options.donutCutString}

da = donutana(**inDict)

# get image directory
imageDir = "%s/Donuts/%ds%d/%d/v%d" % (options.rootDirectory,options.date,options.sequence,options.expoNumber,options.versionNum)
print("runImageDonutAna:  analyzing image directory ",imageDir)

# dictionary for results
anaDict = OrderedDict()

fititerStr = {1:"first",2:"second",3:"third",4:"fourth"}

# if we are using scratchDir get files there, otherwise in the normal directory
if options.scratchDir!="":   
    scratchDirExp = os.path.expandvars(options.scratchDir)     
    scratchFileString = scratchDirExp + "/*.%s.donut.fits*" % (fititerStr[options.iteration])
    donutFileList = glob.glob(os.path.expandvars(scratchFileString))
    # try in non-scratch area too
    if len(donutFileList)==0:
        directory = imageDir + "/*.%s.donut.fits*" % (fititerStr[options.iteration])
        donutFileList = glob.glob(os.path.expandvars(directory))
        
else:
    directory = imageDir + "/*.%s.donut.fits*" % (fititerStr[options.iteration])
    donutFileList = glob.glob(os.path.expandvars(directory))

# print files used
print("runImageDonutAna")
print(donutFileList)

# list of donut headers
donutHeaderList = []

# read in donut fit files and build list of headers
for donutFile in donutFileList:
    hdu = pyfits.open(donutFile)
    header = hdu[0].header
    donutHeaderList.append(header)
        
# now call donutana
anaDict[options.expoNumber] = da.analyzeDonuts(donutHeaderList,doCull=True,cullCut=0.2)
    
# pickle anaDict
# now use specific pickle name...
pickleFileName = imageDir + "/donutana_%s-v%di%d.pickle" % (options.expoNumber,options.versionNum,options.iteration)
pickleFileName = os.path.expandvars(pickleFileName)
print("runImageDonutAna:  writing pickle ",pickleFileName)
pickle.dump( anaDict, open(pickleFileName, "wb" ) )
print("runImageDonutAna:  Done ")
