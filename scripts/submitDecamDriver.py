#! /usr/bin/env python
#
# $Rev:: 222                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2018-07-18 17:26:50 -0700 (Wed, 18 Jul 2018)     $:  
#
#
#  This scripts has NO DEFAULTS of its own!!!
#  it uses the underlying defaults ONLY
#
import os
import pdb
from scriptUtil import decodeNumberList,parsePipelineArguments

# get all command line arguments
extra_definition = {"imageListString":['s',"n","image list string, format eg. 234:299,307,400:450"],
                    "date":["i","d","date with format 20140212"],
                    "sequence":["i","seq","image sequence number"],
                    "versionNum":["i","ver","version number"]}

options = parsePipelineArguments(extras=extra_definition)

if options.queue =='None':
    if options.doFandA or options.doPlus or options.doMinus :
        options.queue = "-W 16:00"
    elif options.doScience:
        options.queue = "-W 36:00"

imageList = decodeNumberList(options.imageListString)

# expand $env for output area directory (ie. 20130829s1)
outputArea = "%s/Donuts/%ds%d" % (options.rootDirectory,options.date,options.sequence)
outputAreaexp = os.path.expandvars(outputArea)

# make directory if it doesn't exist
if not os.path.exists(outputAreaexp):
    os.mkdir(outputAreaexp)

for i in imageList:
    # check to be sure image is present, or we are running with noSnip == True
    #
    imageFile = "%s/desdata/%d/DECam_00%d.fits.fz" % (options.rootDirectory,options.date,i)
    imageFile = os.path.expandvars(imageFile)
    if os.path.exists(imageFile) or options.noSnip:

        # make directory for log files...
        logdir = "%s/Donuts/logfiles/%s" % (options.rootDirectory,options.date)
        if not os.path.exists(logdir):
            os.mkdir(logdir)

        command = "bsub -C 0 %s -o %s/Donuts/logfiles/%d/DECam_00%d.log decamDriver.py -i %d -o %s/%d -c %s -seq %d -ver %d -d %s %s" % (options.queue,options.rootDirectory,options.date,i,i,outputAreaexp,i,options.configFile,options.sequence,options.versionNum,options.date,options.remaining)

        print("Here is the command:")
        print(command)
        
        if options.nojobFlag:
            print(command)
        else:
            os.system(command)
    else:
        print("submitDecamDriver: no file ",imageFile)
        
