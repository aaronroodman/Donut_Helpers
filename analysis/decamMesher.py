#! /usr/bin/env python
#
# filter input DataFrame with cuts, output new DataFrame and Meshes in Ascii format
#
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from astropy.stats import median_absolute_deviation
import argparse
import copy
import pdb

from collections import OrderedDict
from donutlib.decamutil import decaminfo

# initialization (are these global?)
dinfo = decaminfo()
infoDict = dinfo.info()


def addMAD(df,extList,zernList):
    """ Add nMAD values to the Data Frame
    """
    ndonuts = df.shape[0]

    # list of DataFrames to concatenate
    dfList = []
    
    # create new columns
    for iZ in zernList:
        df['nmad%d' % (iZ)] = 0.0

    if ndonuts>0:
        # loop over terms and sensors
        for iCCD in extList:
            df_sensor = (df.loc[(df.EXTNAME == iCCD)]).copy()
            for iZ in zernList:
                varName = 'ZERN%d' % (iZ)
                theMedian = np.median(df_sensor[varName])
                theMad = median_absolute_deviation(df_sensor[varName])
                df_sensor.loc[:,'nmad%d' % (iZ)] = (df_sensor[varName] - theMedian)/theMad  #gives a warning message but I think this is ok
            dfList.append(df_sensor)
        # make new DataFrame
        newdf = pd.concat(dfList)
    else:
        newdf = df
            
    return newdf

def addNNMAD(df,zernList):
    """ Add nNNMAD values to the Data Frame
    """

    ndonuts = df.shape[0]

    xyArr = df.loc[:,['XDECAM','YDECAM']]
    kdtree = cKDTree(xyArr)
    
    # set k based on the number of entries
    if ndonuts>200:
        nNeighbor = 12
    if ndonuts>30:
        nNeighbor = 6
    elif ndonuts>10:
        nNeighbor = 4
    else:
        nNeighbor = 1  # ie. it will only find itself

    # add columns to the dataFrame
    for iZ in zernList:
        df['nmadNN%d' % (iZ)] = 0.0

    # fill nmadNN
    for idonut in range(ndonuts):
        # find NN neighbors
        xyVal = df.loc[idonut,['XDECAM','YDECAM']]
        d,ind = kdtree.query(xyVal,nNeighbor)
        df_near = df.loc[ind]

        # loop over terms
        for iZ in zernList:
            varName = 'ZERN%d' % (iZ)
            theMedian = np.median(df_near[varName])
            theMad = median_absolute_deviation(df_near[varName])
            df.loc[idonut,'nmadNN%d' % (iZ)] = (df.loc[idonut,varName] - theMedian)/theMad

        # delete df_near object
        #print(idonut,df.loc[idonut,'nmadNN5'],df.shape)
        del df_near

    return 

def addToData(df,extList,zernList):
    """ add some information to the DataFrame
    """

    # add zern4 converted to focal distance
    df['Z4'] = df.ZERN4 * (1500./8.7)
    df['Z4E'] = df.ZERN4E * (1500./8.7)

    # calculate maxMAD and maxMADNN
    newdf = addMAD(df,extList,zernList)
    addNNMAD(newdf,zernList)

    return newdf

def passCuts(df,cutDict):

    # cuts
    neleCut = cutDict["neleCut"]
    sumsqzernCut = cutDict["sumsqzernCut"]
    ncalcallCut = cutDict["ncalcallCut"]
    rdecamCut = cutDict["rdecamCut"]
    absz4LoCut = cutDict["absz4LoCut"]
    absz4HiCut = cutDict["absz4HiCut"]
    madCut = cutDict["nMADCut"]
    madNNCut = cutDict["nMADNNCut"]

    nelePass = np.log10(df.NELE) > neleCut
    sumsqzernPass = np.sqrt(np.power(df.ZERN5,2)+np.power(df.ZERN5,2)+np.power(df.ZERN5,2)+np.power(df.ZERN5,2)) < sumsqzernCut
    ncalcallPass = df.NCALCALL < ncalcallCut
    rdecamPass = np.sqrt(np.power(df.XDECAM,2)+np.power(df.YDECAM,2)) < rdecamCut
    z4LoPass = np.abs(df.ZERN4) > absz4LoCut
    z4HiPass = np.abs(df.ZERN4) < absz4HiCut
    madPass = (df.filter(regex='nmad[0-9]').min(axis=1) > -madCut) & (df.filter(regex='nmad[0-9]').max(axis=1) < madCut)
    madNNPass = (df.filter(regex='nmadNN[0-9]').min(axis=1) > -madNNCut) & (df.filter(regex='nmadNN[0-9]').max(axis=1) < madNNCut)

    allPass = nelePass & sumsqzernPass & ncalcallPass & rdecamPass & z4LoPass & z4HiPass & madPass & madNNPass

    #df['nelePass'] = nelePass
    #df['sumsqzernPass'] = sumsqzernPass
    #df['ncalcallPass'] = ncalcallPass
    #df['rdecamPass'] = rdecamPass
    #df['z4LoPass'] = z4LoPass
    #df['z4HiPass'] = z4HiPass
    #df['madPass'] = madPass
    #df['madNNPass'] = madNNPass
    #df['allPass'] = allPass
    
    #return df
    
    df_pass = df[allPass]
    return df_pass
    
    

def doMesher(optDict):
    """ makes Mesh data files    
    """
    # get cut values
    paramDict = {}
    paramDict.update(optDict)

    # get the DataFrame
    df = pd.read_pickle(paramDict["inputFile"])

    # form extension list
    # check first for command line lists, otherwise use defaults
    if paramDict["extNamesInput"] ==  None:
        extList = []
        if paramDict["doFandA"]:
            for extName in infoDict:
                if infoDict[extName]["FAflag"]:
                    extList.append(extName)

        if paramDict["doScience"]:
            for extName in infoDict:
                if not infoDict[extName]["FAflag"]:
                    extList.append(extName)

        if paramDict["doPlus"]:
            for extName in infoDict:
                if  infoDict[extName]["FAflag"] and infoDict[extName]["Offset"] == 1500.0:
                    extList.append(extName)

    else:
        extList = paramDict["extNamesInput"]

    if len(extList)==0:
        print("decamMesher: HEY YOU FOOL, you need to specify either -doFandA or -doScience")
        exit

    #print("decamMesher extensions: ",extList)

    # loop over files
    resultsDict = {}
    for ifile in paramDict["fileList"]:

        # filter DataFrame by ifile
        df_file = df.loc[(df.IFILE == ifile)]

        # add rows
        df_add = addToData(df_file,extList,paramDict["zernList"])

        # filter with cuts
        df_pass = passCuts(df_add,paramDict)

        # write out the DataFrame
        df_output = df_pass.reset_index()
        df_output.to_pickle(paramDict["outputDirectory"] + "/" + paramDict["outputPrefix"] + "_%d.pkl" % (ifile))

        # write out the Meshes
        df_output['wgt'] = 1.0
        for iZ in paramDict["zernList"]:
            meshFileName = paramDict["outputDirectory"] + "/" + "z%dMesh_%s_%d.dat" % (iZ,paramDict["outputPrefix"],ifile)
            if iZ == 4:
                df_output.to_csv(meshFileName,sep=' ',columns=['EXTNAME','XDECAM','YDECAM','Z%d' % (iZ),'wgt'],float_format='%12.5f',index=False,header=False,quotechar=' ')
            else:
                df_output.to_csv(meshFileName,sep=' ',columns=['EXTNAME','XDECAM','YDECAM','ZERN%d' % (iZ),'wgt'],float_format='%12.5f',index=False,header=False,quotechar=' ')

    # done
    return
        
             
# if calling from the command line, use the argument parser    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='decamMesher')
    parser.add_argument("-i", "--inputFile",
                        dest="inputFile",
                        default=None,
                        help="input root file ")
    parser.add_argument("-f", "--fileList",
                        dest="fileList",nargs='*',type=int,
                        help="list of image file numbers, one mesh per file")
    parser.add_argument("-o", "--outputPrefix",
                        dest="outputPrefix",
                        default=None,
                        help="output File Prefix ")
    parser.add_argument("-odir", "--outputDirectory",
                        dest="outputDirectory",
                        default=None,
                        help="output File Directory ")
    parser.add_argument("-csv", "--csvFile",
                        dest="csvFile",
                        default=None,
                        help="csv File Name ")

    parser.add_argument("-e","--extNamesInput",nargs='*',
                    dest="extNamesInput",
                    default=None,
                    help="List of Extensions - overrides defaults")

    parser.add_argument("-doFandA", "--doFandA",
                        dest="doFandA",action="store_true",default=False,
                        help="do Focus and Alignment chips")

    parser.add_argument("-doScience", "--doScience",
                        dest="doScience",action="store_true",default=False,
                        help="do Science chips")

    parser.add_argument("-doPlus", "--doPlus",
                        dest="doPlus",action="store_true",default=False,
                        help="do F&A +1500 chips")

    parser.add_argument("-nMADCut", "--nMADCut",
                        dest="nMADCut",
                        default=5.0,type=float,
                        help="cut on MAD ")

    parser.add_argument("-nMADNNCut", "--nMADNNCut",
                        dest="nMADNNCut",
                        default=8.0,type=float,
                        help="cut on MADNN")

    parser.add_argument("-ncalcallCut", "--ncalcallCut",
                        dest="ncalcallCut",type=int,
                        default=200,
                        help="cut on ncalcall")
    parser.add_argument("-absz4LoCut", "--absz4LoCut",
                        dest="absz4LoCut",type=float,
                        default=4.0,
                        help="cut on minimum z4")
    parser.add_argument("-absz4HiCut", "--absz4HiCut",
                        dest="absz4HiCut",type=float,
                        default=15.0,
                        help="cut on maximum z4")
    parser.add_argument("-neleCut", "--neleCut",
                        dest="neleCut",
                        default=3.5,type=float,
                        help="cut on log10(nele)")
    parser.add_argument("-sumsqzernCut", "--sumsqzernCut",
                        dest="sumsqzernCut",
                        default=3.0,type=float,
                        help="cut on sum of squares of Zernikes 5 through 8")
    parser.add_argument("-zernList", "--zernList",
                        dest="zernList",nargs='*',
                        default=[4,5,6,7,8,9,10],type=int,
                        help="zernike terms to write out")
    parser.add_argument("-rdecamCut", "--rdecamCut",
                        dest="rdecamCut",
                        default=225.,type=float,
                        help="cut on decam Radius")
    
    # collect the options 
    options = parser.parse_args()
    optDict = vars(options)

    # call it!
    doMesher(optDict)

