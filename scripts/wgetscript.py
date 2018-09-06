#! /usr/bin/env python
#
# $Rev:: 212                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2015-08-17 12:07:18 -0700 (Mon, 17 Aug 2015)     $:  
#
import argparse
import os
from scriptUtil import decodeNumberList

parser = argparse.ArgumentParser(prog='wgetscript')
parser.add_argument("-d", "--date", dest="date",type=int,
                  help="date")
parser.add_argument("-min", "--minImage", dest="minImage",type=int,default=None,
                  help="minimum image no")
parser.add_argument("-max", "--maxImage", dest="maxImage",type=int,default=None,
                  help="maximum image no")
parser.add_argument("-n", "--imageListString", dest="imageListString",default=None,
                  help="string with image lists, eg. 1:10,20,30:100 ")


# unpack options
options = parser.parse_args()

# get on the directory, make it if it isn't there
#
dataDirectory = "/nfs/slac/g/ki/ki06/roodman/desdata/%d" % (options.date)
dataDirectoryExp = os.path.expandvars(dataDirectory)

# make directory if it doesn't exist
if not os.path.exists(dataDirectoryExp):
    os.mkdir(dataDirectoryExp)

# move there!
os.chdir(dataDirectoryExp)

if options.imageListString!=None:
    imageList = decodeNumberList(options.imageListString)
else:
    imageList = list(range(options.minImage,options.maxImage+1))


# find correct location at NCSA - changed again, but all links are here
# obsolete: weblocation = "https://desar2.cosmology.illinois.edu/DESFiles/desardata/DTS/src/%d/src" % (options.date)
weblocation = "https://desar2.cosmology.illinois.edu/DESFiles/desarchive/DTS/raw/%d" % (options.date)

    
for run in imageList:

    command = "wget --no-check-certificate --http-user=roodman --http-password=roo70chips -nc -nd -nH -r -k -p -np -nv --cut-dirs=3 %s/DECam_00%d.fits.fz" % (weblocation,run)

###  old command using anon ftp
###command = "wget ftp://desar.cosmology.illinois.edu/DESFiles/desardata/%s/src/%d/src/DECam_00%d.fits.fz" % (dtsName,options.date,run)
    os.system(command)
