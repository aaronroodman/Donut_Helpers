import argparse
import configparser
import sys
import pdb
import os
import glob

class bunch(object):
    """ convert a dictionary to a Namespace
    """
    def __init__(self, adict):
        self.__dict__.update(adict)


def decodeNumberList(inputString):
    """  This routine decodes a compound list of numbers
    which contains commas and colons.  eg is "141,145,147:150"
    and this should return a list with [141,145,147,148,149,150].
    """

    # this makes a python generator - kinda like an iterator
    ranges = (x.split(":") for x in inputString.split(","))
    # this uses the generator in a double for loop
    theList =  [i for r in ranges for i in range(int(r[0]),int(r[-1]) + 1)]
    return theList


def parsePipelineArguments(extras={}):
    """ Parse the command line, with both a configuration file and explicit options
    """

    # store the command line input
    argv = sys.argv
    argstr = ""
    for sitem in argv[1:]:
        argstr = argstr + sitem + " " 

    # setup the Argument Parser, start with configuration file
    arg_parser = argparse.ArgumentParser(
        description=__doc__,                                          # printed with -h/--help
        formatter_class=argparse.RawDescriptionHelpFormatter,         # Don't mess with format of description
        add_help=False         # Turn off help, so we print all options in response to -h
        )
    arg_parser.add_argument("-c","--configFile",help="Configuration file", metavar="FILE")

    # add extra arguments
    typeDict = {'b':bool,'s':str,'i':int,'f':float,'l':str,'d':str}
    for item in extras:
        keyList = extras[item]
        keyType = keyList[0]
        keyShortCut = keyList[1]
        keyHelp = keyList[2]

        # print "section,item:",section,item,keyHelp,typeDict[keyType]
        arg_parser.add_argument("-%s" % (keyShortCut),"--%s" % (item),help=keyHelp,type=typeDict[keyType])

    # if we have a configuration file, parse it
    cloptions, remaining_options = arg_parser.parse_known_args()

    # add input command line and any extra options to options dictionary
    optionsDict = {"argstr":argstr}
    optionsDict.update(vars(cloptions))

    # if any of the extra options are lists or dictionaries, decode them now...
    for item in extras:
        keyList = extras[item]
        keyType = keyList[0]
        if keyType == 'l' or keyType == 'd' :
            # now convert to a dict or list using Python magic
            commstr = "optionsDict[item] = cloptions.%s" % (item)
            exec(commstr)


    # define sections, variables, types for configuration file parameters
    
    # define sections
    config_definition = {"General":{},"Image":{},"Sextractor":{},"DonutEngine":{},"Analysis":{}}

    # in each section, define variables and setup their information
    # this includes their type:  b is Boolean, i is Integer, f is Float, s is String, l is List, d is Dictionary
    # their default value and the help text
    config_definition["General"] = {"nojobFlag":['b',False,"if True jobs are not submitted"],
                                    "noSnip":['b',False,"if True donuts are not snipped"],
                                    "noFit":['b',False,"if True donuts are not fit"],
                                    "noAna":['b',False,"if True donuts are not analyzed"],
                                    "subJobByExt":['b',False,"if True sub jobs are submitted to fit each extension"],
                                    "queue":['s',"-W 4:00","LSF batch queue"],
                                    "rootDirectory":['s',"/u/ec/roodman/kipacdisk","root level directory"],
                                    "scratchDir":['s',"","scratch directory"],
                                    "printLevel":['i',0,"print level, 0 is none, 1 is some, 2 is lots"],
                                    "debugFlag":['b',False,"debug flag"],
                                    "doFandA":['b',False,"analyze FandA chips"],
                                    "doScience":['b',False,"analyze Science chips"],
                                    "doPlus":['b',False,"analyze Plus FandA chips"],
                                    "doMinus":['b',False,"analyze Minus FandA chips"],
                                    "extNamesInput":['l','None',"list of extentions"]}

    config_definition["Image"] = {"biasFile":['s',"","bias fits file"],
                                  "flatFile":['s',"","flat fits file"],
                                  "noSplitMEF":['b',True,"if True files are not run through splitMEF script"],
                                  "noOverScan":['b',False,"if True files are not overscan corrected"],
                                  "doFilter":['b',False,"apply donut filter"],
                                  "donutFilterFile":['s',"donut-8.7.fits","donut filter file"],
                                  "nOverscanCol":['i',56,"number of overscan columns for DECam fits files"]}

    config_definition["Sextractor"] = {"sexConfig":['s',"~/Astrophysics/Sextractor/decamNoFilter.conf","sex config file"],
                                       "cutString":['s', "flux_auto>1.0e3 and flags==0 and ellipticity<.3 and rdecam<225.0 and ave_image>5","cuts applied to sextractor donut catalog"],
                                       "doRegion":['b',False,"output region files"]}

    config_definition["DonutEngine"] = {"iTelescope":['i',0,"DonutEngine telescope presets, 0 DECam, 1 MosaicII, 2 LSST, 3 Magellan IMACS, 4 Magellan MegaCam"],
                                        "nbin":['i',256,"FFT grid size"],
                                        "nPixels":['i',64,"donut postage stamp pixel size"],
                                        "pixelOverSample":['i',4,"oversampling factor"],
                                        "scaleFactor":['d',2.0,"DonutEngine wavefront scaling factor"],
                                        "rzero":['d',0.125,"default rzero value"],
                                        "zern4":['d',11.0,"default zern4 value"],
                                        "nFits":['i',1,"number of donut fit iterations"],
                                        "nZernikeTerms":['i',11,"Zernike polynomial order"],
                                        "fixedParamArray1":['l',"[0,1,0,0,0,0,0,0,0,0,0,0,1]","fixed/floating parameter mask iteration 1"],
                                        "fixedParamArray2":['l',"[0,1,0,0,0,0,0,0,0,0,0,0,1]","fixed/floating parameter mask iteration 2"],
                                        "fixedParamArray3":['l',"[0,1,0,0,0,0,0,0,0,0,0,0,1]","fixed/floating parameter mask iteration 3"],
                                        "inputZernikeDict":['d',"{'None':[0.0,0.0,0.0],\
                                                                 'FS1':[0.0,0.0,11.0,0.0,0.0,0.0,0.0,0.20,-0.17,-0.08],\
                                                                 'FS2':[0.0,0.0,-11.0,0.0,0.0,0.0,0.0,0.26,-0.01,-0.13],\
                                                                 'FS3':[0.0,0.0,11.0,0.0,0.0,0.0,0.0,0.05,0.25,-0.11],\
                                                                 'FS4':[0.0,0.0,-11.0,0.0,0.0,0.0,0.0,-0.05,0.25,-0.14],\
                                                                 'FN1':[0.0,0.0,-11.0,0.0,0.0,0.0,0.0,0.20,0.17,-0.08],\
                                                                 'FN2':[0.0,0.0,11.0,0.0,0.0,0.0,0.0,0.26,0.01,-0.13],\
                                                                 'FN3':[0.0,0.0,-11.0,0.0,0.0,0.0,0.0,0.05,-0.25,-0.11],\
                                                                 'FN4':[0.0,0.0,11.0,0.0,0.0,0.0,0.0,-0.05,-0.25,-0.14]}",
                                                            "dictionary with mapping from extension to initial parameter values"],
                                        "gain":['d',1.0,"Gain in e/ADU"],
                                        "wavefrontMapFile":['s',None,"Name of wavefront Map File"]
                                        }


    config_definition["Analysis"] = {"donutCutString":['s',"nele>0. and numpy.log10(nele)>3.5 and numpy.sqrt(zern5*zern5+zern6*zern6+zern7*zern7+zern8*zern8)<3.0 and abs(zern4)>5 and abs(zern4)<15","Donut cut string"],
                                     "ndonutscut":['i',10,"Minimum number of donuts per chip"],
                                     "doTrefoil":['b',True,"analyze Trefoil"],
                                     "doSpherical":['b',False,"analyze Spherical"],
                                     "doQuadrefoil":['b',False,"analyze Quadrefoil"],
                                     "doRzero":['b',False,"analyze Rzero"],
                                     "meshDir":['s',"","Zernike mesh directory"],
                                     "meshName1":['s',"","Zernike mesh basename i1"],
                                     "meshName2":['s',"","Zernike mesh basename i2"],
                                     "meshName3":['s',"","Zernike mesh basename i3"],
                                     "sensorSet":['s',"FandAOnly","set of sensors for analysis: FandAOnly, ScienceOnly"]}



    # decode the configuration file, fill configDict, using config_definition
    if cloptions.configFile:
        config_parser = configparser.SafeConfigParser()
        config_parser.optionxform=str  #keeps option name case 
        okfiles = config_parser.read(cloptions.configFile)

        # loop over the sections and decode each variable
        for section in config_definition:
            section_dict = config_definition[section]
            for item in section_dict:
                keyList = section_dict[item]
                keyType = keyList[0]
                keyDefault = keyList[1]
                keyHelp = keyList[2]

                # now get the desired value from the config file
                if keyType == 'i':
                    optionsDict[item] = config_parser.getint(section,item)
                elif keyType == 'f':
                    optionsDict[item] = config_parser.getfloat(section,item)
                elif keyType == 'b':
                    optionsDict[item] = config_parser.getboolean(section,item)
                elif keyType == 's' :
                    optionsDict[item] = config_parser.get(section,item)
                elif keyType == 'l' or keyType == 'd' :
                    optionsDict[item] = config_parser.get(section,item)
    else:
        print("No configuration file specified!")

    # this second parser has all the arguments, meant for overwriting parameters from config file
    morearg_parser = argparse.ArgumentParser(parents=[arg_parser])

    # decode parameters, and add all options to the Argument parser
    typeDict = {'b':bool,'s':str,'i':int,'f':float,'l':str,'d':str}
    for section in config_definition:
        section_dict = config_definition[section]
        for item in section_dict:
            keyList = section_dict[item]
            keyType = keyList[0]
            keyDefault = keyList[1]
            keyHelp = keyList[2]

            # print "section,item:",section,item,keyHelp,typeDict[keyType]
            morearg_parser.add_argument("--%s" % (item),help=keyHelp,type=typeDict[keyType])
            # the logic for parameters is:  config file, command line (defaults in config_definition are NOT used)

    # get all additional command line options
    more_options, moreremaining_options = morearg_parser.parse_known_args(remaining_options)

    # put remaining command line options into a variable - both those used here and not used here
    # so these get passed to next level in cascade of scripts
    remaining = ""
    for sitem in remaining_options:
        remaining = remaining + sitem + " " 
    optionsDict["remaining"] = remaining
    print(remaining)

    # without defaults above, more_options are None, unless set at the command line
    # decode all non-None options
    moreDict = vars(more_options)
    for item in moreDict:
        if moreDict[item]:
            optionsDict[item] = moreDict[item]

    # convert all list or dictionary arguments from strings
    for section in config_definition:
        section_dict = config_definition[section]
        for item in section_dict:
            keyList = section_dict[item]
            keyType = keyList[0]
            keyDefault = keyList[1]
            keyHelp = keyList[2]

            # now get the desired value from the config file
            if keyType == 'l' or keyType == 'd' :
                # convert from a string to a python object
                tempstr = optionsDict[item]
                # now convert to a dict or list using Python magic
                commstr = "optionsDict[item] = %s" % (tempstr)
                exec(commstr)

    # done, return as a namespace
    options = bunch(optionsDict)
                                    
    return options

def cleanDonutsDir(image,date,sequence,rootDirectory):
    ### compress some of the fits files in Donuts image directory
    ###

    outputArea = "%s/Donuts/%ds%d/%d" % (rootDirectory,date,sequence,image)
    outputAreaexp = os.path.expandvars(outputArea)

    # compress corr files
    try:
        os.chdir(outputAreaexp)
    except:
        print("cleanDonutsDir: cannot find %s" % (outputAreaexp))
        return

    # remove .corr.fits.fz if they already exist
    command = "rm *.corr.fits.fz"
    os.system(command)
        
    command = "fpack -D -Y *.corr.fits"
    os.system(command)

    # compress fits files in vN directories
    listvN = glob.glob("v*")
    for vdir in listvN:
        try:
            os.chdir("%s/%s" % (outputAreaexp,vdir))
            command = "fpack -D -Y *.donut.fits"
            os.system(command)
        except:
            print("cleanDonutsDir: cannot find %s/%s" % (outputAreaexp,vdir))
            return

def dequote(s):
    """
    If a string has single or double quotes around it, remove them.
    Make sure the pair of quotes match.
    If a matching pair of quotes is not found, return the string unchanged.
    """
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return s[1:-1]
    return s
