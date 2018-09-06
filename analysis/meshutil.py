#
# $Rev:: 205                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2015-05-20 10:29:08 -0700 (Wed, 20 May 2015)     $:
#
import os
import numpy
import numpy.lib.index_tricks as itricks
import scipy
import statsmodels.api as sm
import copy
import pickle
import pdb
import csv
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm,colors
from donutlib.PointMesh import PointMesh
from donutlib.donutana import donutana
from donutlib.decamutil import decaminfo
from scriptUtil import decodeNumberList


maxZernikeTerm = 15


### Utility routines for meshes

def adjustAllMeshes(meshDict,donut_summary,zangleconv=206.264806247):
    """ utility routine to adjust all the meshes in a dict from analyzeDonut
    NOTE that this guy merges INPLACE, and it adjusts w/ 3 DOF for each mesh
    """
    
    for i in range(4,maxZernikeTerm+1):
       meshName = "z%dMesh" % (i)
       if meshName in meshDict:        
           thisMesh = meshDict[meshName]
           thetax = donut_summary["z%dthetax" % (i)]
           thetay = donut_summary["z%dthetay" % (i)]
           delta = donut_summary["z%ddelta" % (i)]
           print("Adjust ",i,thetax,thetay,delta)
           if i==4:
               thisMesh.adjustMesh(thetax,thetay,delta,angleconversion=zangleconv)
           else:
               thisMesh.adjustMesh(thetax,thetay,delta) 


def fitflatMesh(aMesh,method="RLM"):
    """ utility routine to fit one mesh to a flat surface: find its mean and x,y slopes
    2/10/14 rewrite to do a single fit in X and Y 
    """

    # get the points 
    X,Y,Z = aMesh.getXYZpoints()
    delta = numpy.mean(Z)

    # column with just 1
    npoints = Z.shape[0]
    oneColumn = numpy.ones(npoints)

    # get slopes using a linear fit (RLM for a robust fit)
    # formula to fit is ThetaX*y + ThetaY*x + DeltaZ
    aMatrix = numpy.vstack([Y,X,oneColumn]).T

    # make sure there is something to fit
    if npoints>3:
        # use Regression method or use sm.RLM for robust fitting
        if method=="OLS":
            linearModel = sm.OLS(Z,aMatrix)
        elif method=="RLM":
            linearModel = sm.RLM(Z,aMatrix)

        try:
            results = linearModel.fit()
            print("anaMesh: ThetaX = ",results.params[0]," +- ",results.bse[0])
            print("         ThetaY = ",results.params[1]," +- ",results.bse[1])
            print("         DeltaZ = ",results.params[2]," +- ",results.bse[2])
            thetaX = results.params[0]
            thetaY = results.params[1]
            delta = results.params[2]

        except:
            print("anaMesh.py:  fit failed")
            results = None
            thetaX = None
            thetaY = None
            delta = None

    else:
        results = None
        thetaX = None
        thetaY = None
        delta = None

    return {"delta":delta,"thetaX":thetaX,"thetaY":thetaY}


def fitflatAllMeshes(donutDict,method="RLM"):
    """ utility routine to analyze the meshes in the dict: to find the mean and x,y slopes
    """
    anaResults = {}

    for i in range(4,maxZernikeTerm+1):
        meshName = "z%dMesh" % (i)
        if meshName in donutDict:
            aMesh = donutDict[meshName]
            anaResults[meshName] = fitflatMesh(aMesh,method=method)

    return anaResults

       

def mergeAllMeshes(dictOfMeshDicts):
    """ utility routine to merge all the meshes together
    """

    first = True
    for key in list(dictOfMeshDicts.keys()):
        if first:            
            outputDict = {}
            firstMeshDict = dictOfMeshDicts[key]

            for iZ in range(4,maxZernikeTerm+1):
                name = "z%dMesh" % (iZ)
                if name in firstMeshDict:
                    outputDict[name] = copy.deepcopy(firstMeshDict[name])

            first = False            
        else:            
            meshDictOther = dictOfMeshDicts[key]

            for iZ in range(4,maxZernikeTerm+1):
                name = "z%dMesh" % (iZ)
                if name in meshDictOther:
                    thisMesh = outputDict[name]
                    otherMesh = meshDictOther[name]
                    thisMesh.mergeMesh(otherMesh)

    return outputDict


def writeAllMeshes(donutDict,fileName="test",directory=""):        
    """ write all meshes to .dat files
    """

    for i in range(4,maxZernikeTerm+1):
        name = "z%dMesh" % (i)
        if name in donutDict:
            thisMesh = donutDict[name]
            meshName = os.path.join(directory,name+"_"+fileName+".dat")
            thisMesh.writePointsToFile(meshName)


def mkDonutAna(meshName,sensorSet,method,methodVal=None,directory="Meshes",donutCutString="",nInterpGrid=8,doTrefoil=True,doSpherical=False,doQuadrefoil=False):
    """ make a donutana object
    """

    # use whatever mesh files exist!
    fileNameDict = {}
    for iZ in range(4,maxZernikeTerm+1):
        fileName = os.path.join(directory,"z%dMesh_%s.dat" % (iZ,meshName))
        if os.path.lexists(fileName):
            fileNameDict[iZ] = fileName
        else:
            fileNameDict[iZ] = ""

    inDict = {"z4PointsFile":fileNameDict[4], 
              "z5PointsFile":fileNameDict[5], 
              "z6PointsFile":fileNameDict[6], 
              "z7PointsFile":fileNameDict[7], 
              "z8PointsFile":fileNameDict[8], 
              "z9PointsFile":fileNameDict[9], 
              "z10PointsFile":fileNameDict[10], 
              "z11PointsFile":fileNameDict[11], 
              "z12PointsFile":fileNameDict[12], 
              "z13PointsFile":fileNameDict[13], 
              "z14PointsFile":fileNameDict[14], 
              "z15PointsFile":fileNameDict[15], 
              "sensorSet":sensorSet,
              "doTrefoil":doTrefoil,
              "doSpherical":doSpherical,
              "doQuadrefoil":doQuadrefoil,
              "unVignettedOnly":False,
              "interpMethod":method,
              "methodVal":methodVal,
              "nInterpGrid":nInterpGrid,
              "histFlag":True,
              "debugFlag":True,               
              "donutCutString":donutCutString}

    da = donutana(**inDict)
    return da

def getCentralVal(da):
    """ calculate mean,rms at the center of the focal plane
    """

    pts = {}
    pts["N4"] = itricks.mgrid[0.+2.0:33.816-2.0:1j*5,-31.945+2.0:31.945-2.0:1j*10]
    pts["S4"] = itricks.mgrid[-33.816+2.0:0.0-2.0:1j*5,-31.945+2.0:31.945-2.0:1j*10]

    coordList = ["N4","S4"]

    # loop over meshes, get mean,rms from N4 and S4 only
    results = {}

    for iZ in range(4,maxZernikeTerm+1):
        values = []
        for iCoord in coordList:
            if "z%dMesh" % (iZ) in da.meshDict:
                val = da.meshDict["z%dMesh" % (iZ) ].doInterp(iCoord,pts[iCoord][0].flatten(),pts[iCoord][1].flatten())
                values.extend(val)
        varr = numpy.array(values)
        mean = varr.mean()
        rms = numpy.sqrt(varr.var())
        results[iZ] = (mean,rms)

    print(results)
    return results


def getMedian(da,iCoord,iZ):
    """ calculate median of points in a given sensor
    """

    medianVal = -999999.9
    if "z%dMesh" % (iZ) in da.meshDict:
        mesh = da.meshDict["z%dMesh" % (iZ) ]
        x,y,z = mesh.getXYZpoints(iCoord)
        medianVal = numpy.median(z)

    return medianVal

def medianZFandA(da):

    info = da.info
    plusAve = []
    minusAve = []

    for iCoord in list(info.keys()):
        if info[iCoord]['FAflag']:
            if info[iCoord]['Offset'] > 0. :
                medianVal = getMedian(da,[iCoord],4)
                plusAve.append(medianVal)
            elif info[iCoord]['Offset'] < 0. :
                medianVal = getMedian(da,[iCoord],4)
                minusAve.append(medianVal)

    plusArr = numpy.array(plusAve)
    minusArr = numpy.array(minusAve)
    deltaZ = (plusArr.mean() + minusArr.mean())/2.0
    return deltaZ

###
### make MPL plots from a set of Meshes
###

def plotAllMeshes(meshName,sensorSet="ScienceOnly",method="idw",methodVal=(250,1.),nInterpGrid=32,directory="ComboMeshes"):

    from meshutil import mkDonutAna
    from matplotlib import pyplot as plt
    plt.interactive(True)

    minDict = {4:-60.,5:-0.4,6:-0.4,7:-0.25,8:-0.25,9:-0.6,10:-0.6,11:-0.3,14:-0.2,15:-0.2}
    maxDict = {4:60.,5:0.4,6:0.4,7:0.25,8:0.25,9:0.6,10:0.6,11:0.3,14:0.2,15:0.2}
    daS = mkDonutAna(meshName,sensorSet,method,directory=directory,methodVal=methodVal,nInterpGrid=nInterpGrid)

    zernList = [4,5,6,7,8,9,10,11,14,15]
    for iZ in zernList:
        meshN = "z%dMesh" % (iZ)
        if meshN in daS.meshDict:
            f,a,c = daS.meshDict[meshN].plotMeshMPL2D(zmin=minDict[iZ],zmax=maxDict[iZ])
            f.savefig(meshName+"_z%d.png" % (iZ))
            if iZ>4:
                f,a,c = daS.meshDict["z%dMesh" % (iZ)].plotMeshMPL2D(zmin=-0.1,zmax=0.1)
                f.savefig(meshName+"_z%d-zoom.png" % (iZ))
            else:
                f,a,c = daS.meshDict["z%dMesh" % (iZ)].plotMeshMPL2D(zmin=-20,zmax=20)
                f.savefig(meshName+"_z%d-zoom.png" % (iZ))




def plotAllMeshesMPL(meshName,sensorSet="ScienceOnly",method="idw",doPlotInterp=True,doPlotColor=True,writePdf=False,outputName="",methodVal=None,directory="Meshes",donutCutString=""):

    plt.interactive(True)

    da = mkDonutAna(meshName,sensorSet,method,methodVal=methodVal,directory=directory,donutCutString=donutCutString)
    
    outDict = {}
    for iZ in range(4,maxZernikeTerm+1):        
        if "z%dMesh" % (iZ) in da.meshDict:
            outDict["z%d" % (iZ)]  = da.meshDict["z%dMesh" % (iZ)].plotMeshMPL(plotInterp=doPlotInterp,plotColor=doPlotColor,title="")

    if writePdf:
        for iZ in range(4,maxZernikeTerm+1):        
            if "z%dMesh" % (iZ) in da.meshDict:
                outDict["z%d" % (iZ)][0].savefig("z%dMesh_%s.pdf" % (iZ,outputName) )

    return da,outDict


def adjustMeshes(donutDict,anaOut):
    """  method to do custom mesh adjustment, make a deep copy of the input mesh and adjust it!
    7/28/2014 now adjust Spherical and Quatrefoil, but adjust theta x,y for Astigmatism using physical constraint
    note then that we don't adjust the Coma theta x,y
    """


    copyDict = copy.deepcopy(donutDict)

    zMesh = copyDict["z4Mesh"]
    zMesh.adjustMesh(0.0,0.0,anaOut["z4Mesh"]["delta"])

    z5thetaX = 0.5 * (anaOut["z5Mesh"]["thetaX"] + anaOut["z6Mesh"]["thetaY"]) 
    z6thetaY = z5thetaX

    z6thetaX = 0.5 * (anaOut["z6Mesh"]["thetaX"] - anaOut["z5Mesh"]["thetaY"]) 
    z5thetaY = -z6thetaX

    z5Mesh = copyDict["z5Mesh"]
    z5Mesh.adjustMesh(z5thetaX,z5thetaY,anaOut["z5Mesh"]["delta"])

    z6Mesh = copyDict["z6Mesh"]
    z6Mesh.adjustMesh(z6thetaX,z6thetaY,anaOut["z6Mesh"]["delta"])

    print("adjustMesh by : 4",anaOut["z4Mesh"]["delta"])
    print("                5",z5thetaX,z5thetaY,anaOut["z5Mesh"]["delta"])
    print("                6",z6thetaX,z6thetaY,anaOut["z6Mesh"]["delta"])

    # only adjust up to Trefoil - don't adjust Spherical or Quatrefoil
    for i in range(7,maxZernikeTerm+1):
        name = "z%dMesh" % (i)
        if name in copyDict:
            theMesh = copyDict[name]
            theMesh.adjustMesh(0.0,0.0,anaOut[name]["delta"])
            print("                ",i,anaOut[name]["delta"])

    return copyDict


def calcIDWval(nDonuts,radius,fraction=0.25):
    """ radius in arcminutes
    Calculate kNN and radiusBaseline for the inverse distance weighted interpolation
    """

    kNN = int( nDonuts / (3600./numpy.power(radius,2)) ) # number of donuts inside radius
    radiusBaseline = ( 0.015 * radius / (0.27/60.) ) * fraction  # quarter of the radius, in mm
    return (kNN,radiusBaseline)



def washMesh(meshNameA,meshNameB,meshNameOut,pickleOut=None,directoryIn="",directoryOut="",method="idw",methodVal=(20,1.0),sensorSet="ScienceOnly"):
    """ This script does the following washing:

    1) fits meshB against A, and adjusts and culls B (meshNameBr1)
    2) takes meshBr1 and then culls (but does not adjust) the A against it (meshNameAr1)
    3) iterates by adjusting and culling B (meshNameBr2) against the culled A (meshNameAr1)
    4) merges meshNameBr2 and meshNameAr1 (meshNameOut)

    The scripts returns it all, and if desired outputs a pickle with all intermediate information. 
    """

    meshNameAr1 = meshNameA + "-r1"
    meshNameBr1 = meshNameB + "-r1"
    meshNameBr2 = meshNameB + "-r2"

    # collect outputs
    outputDict = {}

    # now clean and adjust mesh B against mesh A
    daA = mkDonutAna(meshNameA,sensorSet,method,methodVal=methodVal,directory=directoryIn)
    daB = mkDonutAna(meshNameB,sensorSet,None,directory=directoryIn)
    meshDictBr1 = daA.analyzeMeshes(daB.meshDict,doCull=True,cullCut=0.20)
    adjustAllMeshes(meshDictBr1,meshDictBr1["donut_summary"])
    writeAllMeshes(meshDictBr1,fileName=meshNameBr1,directory=directoryOut)

    outputDict["anaBr1"] = meshDictBr1

    # now clean mesh A against cleaned/adjusted mesh B    
    daBr1 = mkDonutAna(meshNameBr1,sensorSet,method,methodVal=methodVal,directory=directoryIn)
    daA = mkDonutAna(meshNameA,sensorSet,None,directory=directoryOut)
    meshDictAr1 = daBr1.analyzeMeshes(daA.meshDict,doCull=True,cullCut=0.20)
    writeAllMeshes(meshDictAr1,fileName=meshNameAr1,directory=directoryOut)

    outputDict["anaAr1"] = meshDictAr1

    # now clean and adjust mesh B against mesh Ar1
    daAr1 = mkDonutAna(meshNameAr1,sensorSet,method,methodVal=methodVal,directory=directoryOut)
    daB = mkDonutAna(meshNameB,sensorSet,None,directory=directoryIn)
    meshDictBr2 = daAr1.analyzeMeshes(daB.meshDict,doCull=True,cullCut=0.20)
    adjustAllMeshes(meshDictBr2,meshDictBr2["donut_summary"])
    writeAllMeshes(meshDictBr2,fileName=meshNameBr2,directory=directoryOut)

    outputDict["anaBr2"] = meshDictBr2

    # merge Br2 with Ar1
    dictToMerge = {"Ar1":meshDictAr1,"Br2":meshDictBr2}
    meshOut = mergeAllMeshes(dictToMerge)
    writeAllMeshes(meshOut,fileName=meshNameOut,directory=directoryOut)

    outputDict["meshOut"] = meshOut

    # if desired pickle it all...
    if pickleOut!=None:
        pickle.dump(outputDict,open(os.path.join(directoryOut,pickleOut),"wb"))

    return outputDict



def washnflatMesh(meshNameA,meshNameB,pickleOut=None,directoryIn="",directoryOut="",method="idw",methodVal=(20,1.0),sensorSet="ScienceOnly"):
    """ This script does the following washing:

    1) compares meshA with a Flat shape, and finds adjustments to make the mesh flat
    2) applies these adjustments to make a Flat version of meshA (meshNameAr1)
    3) fits meshB against the flattened version of A, and adjusts and culls B (meshNameBr1)
    4) takes meshBr1 and then culls (but does not adjust) the flattened A against it (meshNameAr2)

    The scripts returns it all, and if desired outputs a pickle with all intermediate information. 
    7/28/2014 - use fitflatAllMeshes now
    """

    meshNameAr1 = meshNameA + "-r1"
    meshNameAr2 = meshNameA + "-r2"
    meshNameBr1 = meshNameB + "-r1"

    # collect outputs
    outputDict = {}

    # flatten mesh A
    daA = mkDonutAna(meshNameA,sensorSet,method,directory=directoryIn)
    anaOut = fitflatAllMeshes(daA.meshDict,method="RLM")    # use Robust fit here!
    meshDictAr1 = adjustMeshes(daA.meshDict,anaOut)
    writeAllMeshes(meshDictAr1,fileName=meshNameAr1,directory=directoryOut)

    # this has the adjustments!
    outputDict["anaAr1"] = anaOut

    # now clean and adjust mesh B against flattened mesh A
    daAr1 = mkDonutAna(meshNameAr1,sensorSet,method,methodVal=methodVal,directory=directoryOut)
    daB = mkDonutAna(meshNameB,sensorSet,None,directory=directoryIn)
    meshDictBr1 = daAr1.analyzeMeshes(daB.meshDict,doCull=True,cullCut=0.20)
    adjustAllMeshes(meshDictBr1,meshDictBr1["donut_summary"])
    writeAllMeshes(meshDictBr1,fileName=meshNameBr1,directory=directoryOut)

    outputDict["anaBr1"] = meshDictBr1

    # now clean mesh A against cleaned/adjusted mesh B    
    daBr1 = mkDonutAna(meshNameBr1,sensorSet,method,methodVal=methodVal,directory=directoryOut)
    daAr1 = mkDonutAna(meshNameAr1,sensorSet,None,directory=directoryOut)
    meshDictAr2 = daBr1.analyzeMeshes(daAr1.meshDict,doCull=True,cullCut=0.20)
    writeAllMeshes(meshDictAr2,fileName=meshNameAr2,directory=directoryOut)

    outputDict["anaAr2"] = meshDictAr2
    
    # check that this mesh is still flat, and save it in the output
    # but don't adjust based on this check
    daAFinal = mkDonutAna(meshNameAr2,sensorSet,method,directory=directoryOut,methodVal=methodVal,nInterpGrid=32)
    anaAFinal = fitflatAllMeshes(daAFinal.meshDict,method="RLM")    # use Robust fit here!

    outputDict["daAFinal"] = daAFinal
    outputDict["anaAFinal"] = anaAFinal

    # if desired pickle it all...
    if pickleOut!=None:
        pickle.dump(outputDict,open(os.path.join(directoryOut,pickleOut),"wb"))

    return outputDict



def combineMeshes(meshNameRef,meshNames,imageList,sensorSet="ScienceOnly",method="idw",methodVal=(20,1.0),directoryIn="",directoryOut="",directoryRef="",pickleOut=None,nameDescr="_All"):
    """ Make a combined mesh by adjusting and culling a series of Meshes against a reference, and then
    merging them together
    """

    daAr2 = mkDonutAna(meshNameRef,sensorSet,method,directory=directoryRef,methodVal=methodVal)
    dictOMeshDict = {}
    for i in imageList:

        print("combineMeshes at image = ",i)
        # now clean and adjust meshes against flattened/culled A
        meshNameB = meshNames + "_%d" % (i)
        daB = mkDonutAna(meshNameB,sensorSet,None,directory=directoryIn)
        meshDictBr1 = daAr2.analyzeMeshes(daB.meshDict,doCull=True,cullCut=0.20)
        adjustAllMeshes(meshDictBr1,meshDictBr1["donut_summary"])
        
        # save this meshDict in a List
        dictOMeshDict[i] = meshDictBr1

    # merge all the Meshes
    mergedMeshes = mergeAllMeshes(dictOMeshDict)
    mergedMeshName = meshNames + nameDescr
    writeAllMeshes(mergedMeshes,fileName=mergedMeshName,directory=directoryOut)

    if pickleOut!=None:
        pickle.dump(dictOMeshDict,open(os.path.join(directoryOut,pickleOut),"wb"))

        
def accumulateMeshes(meshNames,imageList,adjustDict,sensorSet="ScienceOnly",directoryIn="",directoryOut=""):
    """ Make a combined mesh by adjusting a series of Meshes according to a predetermined set of delta,thetax,y
    and then accumulating the result.  No culling is done.
    """

    dictOMeshDict = {}
    for i in imageList:

        # build mesh
        meshName = meshNames + "_%d" % (i)
        da = mkDonutAna(meshName,sensorSet,None,directory=directoryIn)

        # adjust the mesh ala the values in adjustDict
        adjustAllMeshes(da.meshDict,adjustDict[i]["donut_summary"])
        
        # save this meshDict in a List
        dictOMeshDict[i] = da.meshDict

    # merge all the Meshes
    mergedMeshes = mergeAllMeshes(dictOMeshDict)
    mergedMeshName = meshNames + "_Accum"
    writeAllMeshes(mergedMeshes,fileName=mergedMeshName,directory=directoryOut)


def compareEdges(sMesh,fMesh,coords=["FS1","FS2","FS3","FS4","FN1","FN2","FN3","FN4"]):
    """ compare edge values of a Science mesh with an FandA mesh -- at all relevant boundaries and 
    return an array with comparisons of the values at those boundaries
    """
    
    dinfo = decaminfo()
    infoDict = dinfo.infoDict

    frac = 63./64.
    mfrac = 2./64.
    xgap = 3.096
    ygap = 2.450
    ccdHalfSize = 15.0e-3 * 1024. 

    # store all needed Edge information here
    # ( [horizontal xlo,xhi,y, neighbor CCD, deltay] , [vertical ylo,yhi,x , neighbor CCD, deltax] )
    dictOfEdges = {"FN4": 
                   (None, 
                    [infoDict["FN4"]["yCenter"] - ccdHalfSize, 
                     infoDict["FN4"]["yCenter"] + ccdHalfSize, 
                     infoDict["FN4"]["xCenter"] - frac*ccdHalfSize,
                     "N30",-(xgap + mfrac*ccdHalfSize)]),
                   "FN3": 
                   (None, 
                    [infoDict["FN3"]["yCenter"] - ccdHalfSize, 
                     infoDict["FN3"]["yCenter"] + ccdHalfSize, 
                     infoDict["FN3"]["xCenter"] - frac*ccdHalfSize,
                     "N30",-(xgap + mfrac*ccdHalfSize)]),
                   "FS4": 
                   (None, 
                    [infoDict["FS4"]["yCenter"] - ccdHalfSize, 
                     infoDict["FS4"]["yCenter"] + ccdHalfSize, 
                     infoDict["FS4"]["xCenter"] + frac*ccdHalfSize,
                     "S30",+(xgap + mfrac*ccdHalfSize)]),
                   "FS3": 
                   (None, 
                    [infoDict["FS3"]["yCenter"] - ccdHalfSize, 
                     infoDict["FS3"]["yCenter"] + ccdHalfSize, 
                     infoDict["FS3"]["xCenter"] + frac*ccdHalfSize,
                     "S30",+(xgap + mfrac*ccdHalfSize)]),
                   "FN2": 
                   ([infoDict["FN2"]["xCenter"] - ccdHalfSize, 
                     infoDict["FN2"]["xCenter"] + ccdHalfSize, 
                     infoDict["FN2"]["yCenter"] - frac*ccdHalfSize,
                     "N31",-(xgap + mfrac*ccdHalfSize)], 
                    [infoDict["FN2"]["yCenter"] - ccdHalfSize, 
                     infoDict["FN2"]["yCenter"] + ccdHalfSize, 
                     infoDict["FN2"]["xCenter"] - frac*ccdHalfSize,
                     "N28",-(xgap + mfrac*ccdHalfSize)]),
                   "FN1": 
                   ([infoDict["FN1"]["xCenter"] - ccdHalfSize, 
                     infoDict["FN1"]["xCenter"] + ccdHalfSize, 
                     infoDict["FN1"]["yCenter"] - frac*ccdHalfSize,
                     "N28",-(xgap + mfrac*ccdHalfSize)], 
                    [infoDict["FN1"]["yCenter"] - ccdHalfSize, 
                     infoDict["FN1"]["yCenter"] + ccdHalfSize, 
                     infoDict["FN1"]["xCenter"] - frac*ccdHalfSize,
                     "N24",-(xgap + mfrac*ccdHalfSize)]),
                   "FS2": 
                   ([infoDict["FS2"]["xCenter"] - ccdHalfSize, 
                     infoDict["FS2"]["xCenter"] + ccdHalfSize, 
                     infoDict["FS2"]["yCenter"] - frac*ccdHalfSize,
                     "S31",-(xgap + mfrac*ccdHalfSize)], 
                    [infoDict["FS2"]["yCenter"] - ccdHalfSize, 
                     infoDict["FS2"]["yCenter"] + ccdHalfSize, 
                     infoDict["FS2"]["xCenter"] + frac*ccdHalfSize,
                     "S28",+(xgap + mfrac*ccdHalfSize)]),
                   "FS1": 
                   ([infoDict["FS1"]["xCenter"] - ccdHalfSize, 
                     infoDict["FS1"]["xCenter"] + ccdHalfSize, 
                     infoDict["FS1"]["yCenter"] - frac*ccdHalfSize,
                     "S28",-(xgap + mfrac*ccdHalfSize)], 
                    [infoDict["FS1"]["yCenter"] - ccdHalfSize, 
                     infoDict["FS1"]["yCenter"] + ccdHalfSize, 
                     infoDict["FS1"]["xCenter"] + frac*ccdHalfSize,
                     "S24",+(xgap + mfrac*ccdHalfSize)]) }

    values = []
    for faCCD in coords:
        edgeH,edgeV = dictOfEdges[faCCD]

        npoints = 16
        if edgeH != None:
            xlo,xhi,y,scienceCCD,deltay = edgeH
            bump = (xhi-xlo)/(npoints*2.)
            for x in numpy.linspace(xlo+bump,xhi-bump,npoints):
                xval = float(x)
                fValue = fMesh.doInterp(faCCD,[xval],[y])
                sValue = sMesh.doInterp(scienceCCD,[xval],[y+deltay])
                deltaValue = fValue - sValue
                values.append((x,y,fValue,sValue,deltaValue))
        if edgeV != None:
            ylo,yhi,x,scienceCCD,deltax = edgeV
            bump = (yhi-ylo)/(npoints*2.)
            for y in numpy.linspace(ylo+bump,yhi-bump,npoints):
                yval = float(y)
                fValue = fMesh.doInterp(faCCD,[x],[yval])
                sValue = sMesh.doInterp(scienceCCD,[x+deltax],[yval])
                deltaValue = fValue - sValue
                values.append((x,y,fValue,sValue,deltaValue))

    v = numpy.array(values)
    return v

def fitadjustEdges(daS,daF,coords=["FS1","FS2","FS3","FS4","FN1","FN2","FN3","FN4"],mode="all"):
    """  This takes a Science mesh and a FandA mesh, and compares them at the edges, fits the differences
and then adjusts the FandA mesh to match - does this for iZ 5 through 15, whichever are present
For Z4 it adjusts the FandA mesh to have median(extrafocal) = -median(intrafocal)
    
    mode can be 'all' or 'physical' - the latter allows only physical adjustments of the wavefront
             or 'physicalc' - which adds rotation of the Coma surfaces
    """

    summaryDict = {}

    iZ = 4
    summaryDict["z%dthetax" % (iZ)] = 0.0
    summaryDict["z%dthetay" % (iZ)] = 0.0
    summaryDict["z%ddelta" % (iZ)] = medianZFandA(daF)
    print("iZ ", iZ, "  Delta = ",summaryDict["z4delta"])

    for iZ in range(5,maxZernikeTerm+1):
        meshName = "z%dMesh" % (iZ)
        if meshName in daS.meshDict and meshName in daF.meshDict:
            sMesh = daS.meshDict[meshName]
            fMesh = daF.meshDict[meshName]
            compEdge = compareEdges(sMesh,fMesh,coords=coords)
            #print iZ,compEdge

            # fit it to a Robust plane
            zDiff = compEdge[:,4]
            x = compEdge[:,0]
            y = compEdge[:,1]
            npoints = x.shape[0]
            ones = numpy.ones(npoints)
            aMatrix = numpy.vstack([y,x,ones]).T
            linearModel = sm.RLM(zDiff,aMatrix)
            results = linearModel.fit()

            print("iZ ", iZ, "  ThetaX = ",results.params[0]," +- ",results.bse[0])
            print("iZ ", iZ, "  ThetaY = ",results.params[1]," +- ",results.bse[1])
            print("iZ ", iZ, "  DeltaZ = ",results.params[2]," +- ",results.bse[2])

            summaryDict["z%dthetax" % (iZ)] = results.params[0]
            summaryDict["z%dthetay" % (iZ)] = results.params[1]
            summaryDict["z%ddelta" % (iZ)] = results.params[2]

    # make the adjustments physical!
    if mode=='physical' :
       aveZern5ThetaX = 0.5 * ( summaryDict["z5thetax"] + summaryDict["z6thetay"] )
       aveZern6ThetaX = 0.5 * ( summaryDict["z6thetax"] - summaryDict["z5thetay"] )
       summaryDict["z5thetax"] = aveZern5ThetaX
       summaryDict["z6thetay"] = aveZern5ThetaX
       summaryDict["z5thetay"] = -aveZern6ThetaX
       summaryDict["z6thetax"] = aveZern6ThetaX

       for iZ in range(7,maxZernikeTerm+1):
           if "z%dthetax" % (iZ) in summaryDict:
               summaryDict["z%dthetax" % (iZ)] = 0.0
               summaryDict["z%dthetay" % (iZ)] = 0.0

       print(summaryDict)

    elif mode=='physicalc' :

       aveZern5ThetaX = 0.5 * ( summaryDict["z5thetax"] + summaryDict["z6thetay"] )
       aveZern6ThetaX = 0.5 * ( summaryDict["z6thetax"] - summaryDict["z5thetay"] )
       summaryDict["z5thetax"] = aveZern5ThetaX
       summaryDict["z6thetay"] = aveZern5ThetaX
       summaryDict["z5thetay"] = -aveZern6ThetaX
       summaryDict["z6thetax"] = aveZern6ThetaX

       for iZ in range(9,maxZernikeTerm+1):
           if "z%dthetax" % (iZ) in summaryDict:
               summaryDict["z%dthetax" % (iZ)] = 0.0
               summaryDict["z%dthetay" % (iZ)] = 0.0

       print(summaryDict)

    adjustAllMeshes(daF.meshDict,summaryDict)
    return summaryDict

def dumpCsv(meshDict,filename="dadata.csv",expid=None,header=True,meshes=[4,5,6,7,8,9,10,11,14,15]):
    """ dumps a csv file with x,y,z4,z5 etc... from the interpolated points in all the meshes
    """
    # collect meshes
    meshList = []
    for i in meshes:
        meshName = "z%sMesh" % (i)
        meshList.append(meshDict[meshName])

    # now loop through data and stack it up
    dataList = []    
    coordList = meshList[0].coordList    
    for iCoord in coordList:
        try:
            x,y = meshList[0].interpGrids[iCoord]
            dataStack = numpy.vstack((x.flatten(),y.flatten()))
            for aMesh in meshList:
                val = (aMesh.interpValues[iCoord]).flatten()
                dataStack = numpy.vstack((dataStack,val))
            # copy from dataStack to the list
            for i in range(dataStack.shape[1]):
                dataList.append(dataStack[:,i])
        except:
            print("dumpCSv: problem at ",iCoord)

    # write csv file
    csvFile = open(filename,"w")

    fields = OrderedDict([('x',None),('y',None)])
    fields['expid'] = None
    for i in meshes:
        fields['z%s' % (i)] = None
    writerOut = csv.DictWriter(csvFile,fieldnames=fields)
    if header:
        writerOut.writeheader()

    for i in range(len(dataList)):
        outDict = {}
        data = dataList[i]
        outDict['x'] = data[0]
        outDict['y'] = data[1]
        if expid!=None:
            outDict['expid'] = expid
        else:
            outDict['expid'] = 0
        j = 0
        for i in meshes:
            outDict['z%s' % (i)] = data[2+j]
            j = j + 1
        writerOut.writerow(outDict)
                        
    # done
    return dataList            


def dumpCsvPlus(meshDict,filename="dadata.csv",expid=None,header=True,meshes=[4,5,6,7,8,9,10,11,14,15]):
    """ dumps a csv file with x,y,z4,z5 etc... from the interpolated points in all the meshes

    This is a mess - switch to using pandas Data Frames??!!!
    """
    # collect meshes
    meshList = []
    for i in meshes:
        meshName = "z%sMesh" % (i)
        meshList.append(meshDict[meshName])

    # now loop through data and stack it up
    dataList = []    
    coordList = meshList[0].coordList    
    for iCoord in coordList:
        if iCoord != 'S30' and iCoord != 'N30':
            x,y = meshList[0].interpGrids[iCoord]
            dataStack = numpy.vstack((x.flatten(),y.flatten()))
            # zern data first
            for aMesh in meshList:
                val = (aMesh.interpValues[iCoord]).flatten()
                dataStack = numpy.vstack((dataStack,val))
            # MAD values
            for aMesh in meshList:
                val = (aMesh.interpMAD[iCoord]).flatten()
                dataStack = numpy.vstack((dataStack,val))
            # Nentry values
            for aMesh in meshList:
                val = (aMesh.interpNentry[iCoord]).flatten()
                dataStack = numpy.vstack((dataStack,val))
            # copy from dataStack to the list
            for i in range(dataStack.shape[1]):
                dataList.append(dataStack[:,i])

    # write csv file
    csvFile = open(filename,"w")

    fields = OrderedDict([('x',None),('y',None)])
    fields['expid'] = None
    for i in meshes:
        fields['z%s' % (i)] = None
        fields['MAD%s' % (i)] = None
        fields['nEntry%s' % (i)] = None


    writerOut = csv.DictWriter(csvFile,fieldnames=fields)
    if header:
        writerOut.writeheader()

    for i in range(len(dataList)):
        outDict = {}
        data = dataList[i]
        outDict['x'] = data[0]
        outDict['y'] = data[1]
        if expid!=None:
            outDict['expid'] = expid
        else:
            outDict['expid'] = 0
        j = 0
        for i in meshes:
            outDict['z%s' % (i)] = data[2+j]
            j = j + 1
        writerOut.writerow(outDict)
                        
    # done
    return dataList            


    

def compareMeshSequence(meshPrefix,imageStr,nInterpGrid=3):

    """ compare a set of meshes...
    compare individual images to the _All, do just iZ=4,5,6
    use the pickle from comboMeshes which are post-adjustment
    use BMedian now as the method!
    """

    plt.interactive(True)
    
    # probably should use nInterpGrid with a value here to match 2n x n bins/CCD for each image
    meshName1 = meshPrefix + "_All"
    da1 = mkDonutAna(meshName1,"ScienceOnly","bmedian",directory="ComboMeshesv20",nInterpGrid=nInterpGrid)
    dumpCsv(da1.meshDict,filename=meshPrefix+"_All.csv",meshes=[4,5,6,7,8,9,10])

    comboPickle = pickle.load(open("ComboMeshesv20/%s_All.pickle" % (meshPrefix)))

    imageList = decodeNumberList(imageStr)

    for i in imageList:

        print("At image ",i)
        meshDict = comboPickle[i]
        meshName2 = "%s_%d" % (meshPrefix,i)

        iZs = [4,5,6,7,8,9,10]
        minDict = {4:-20.,5:-0.2,6:-0.2,7:-0.075,8:-0.075,9:-0.05,10:-0.05,11:-0.05,14:-0.05,15:-0.05}
        maxDict = {4:20.,5:0.2,6:0.2,7:0.075,8:0.075,9:0.05,10:0.05,11:0.05,14:0.05,15:0.05}
        diffMeshDict = {}
        for iZ in iZs:
            meshName = "z%dMesh" % (iZ)
            if meshName in da1.meshDict and meshName in meshDict:

                # need to rebuild the interpolation for the mesh in the pickle, must start all over from the points!

                aMesh = meshDict[meshName]
                x,y,z = aMesh.getXYZpoints()
                oldGridArray = aMesh.gridArray
                coordList = aMesh.coordList
                gridArray = {}
                for icoord in coordList:
                    # remake as a 2n x n grid in each CCD
                    grid = oldGridArray[icoord]
                    grid[0] = 2*nInterpGrid
                    grid[3] = nInterpGrid
                    gridArray[icoord] = grid
                aNewMesh = PointMesh(coordList,gridArray,aMesh.pointsArray,myMethod='bmedian')

                diffMesh = aNewMesh.subtractMesh(da1.meshDict[meshName])
                diffMeshDict[meshName] = diffMesh

                f,a,c = diffMesh.plotMeshMPL2D(zmin=minDict[iZ],zmax=maxDict[iZ],cmap=cm.jet)
                f.savefig(meshName1+"-minus-"+meshName2+"_bmedian_z%d.png" % (iZ))
                plt.close()

        dumpCsv(diffMeshDict,filename=meshName1+"-minus-"+meshName2+".csv",expid=i,header=False,meshes=[4,5,6,7,8,9,10])

