import os
import pdb
import pickle
import numpy as np
import pandas as pa
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm,colors

from scriptUtil import waitForBatch,decodeNumberList
from meshutil import mkDonutAna,mkDonutAnaDF,writeAllMeshes,combineMeshes,plotAllMeshes,getDonutSummary
from donutlib.PointMesh import PointMesh


def makeindivMeshes(date,sequence,imageStr,iterations=[1,2],version=20):
    """  make Meshes for an out of focus sequence
    """

    imageList = decodeNumberList(imageStr)
    # loop over iterations (usually just 1 and 2)
    for iter in iterations:
        
        meshName = "Science-%ds%d-v%di%d" % (date,sequence,version,iter)

        if iterations==2:
            zstr = "-zernList 4 5 6 7 8 9 10 11 14 15"
        else:
            zstr = " "
            
        # make the individual meshes
        for i in imageList:
            command = "bsub -R rhel60 -W 8:00 -C 0 -o logfiles/mesher-%d.log ~/Astrophysics/Code/analysis/decamMesher.py -i rootfiles/%s.root -o %s -f %d -doScience -nMADCut 20. -nMADNNCut 50.0 %s" % (i,meshName,meshName,i,zstr)
            os.system(command)

            
def buildcomboMeshes(date,sequence,imageStr,iterations=[1,2],version=20):
    """  build a combined mesh, and then make plots for each
    Zernike term of the individual mesh vs. the combined one.  Also make plots for the combined mesh.
    """

    plt.interactive(True)

    imageList = decodeNumberList(imageStr)

    # loop over iterations (usually just 1 and 2)
    for iter in iterations:
        
        meshNames = "Science-%ds%d-v%di%d" % (date,sequence,version,iter)

        ##
        ## build a combo mesh 
        ##

        meshNameAr2 = "Science-20131116s2-v20i1_255171-r2"
        combineMeshes(meshNameAr2,meshNames,imageList,method="idw",methodVal=(20,1.0),directoryIn="Meshesv20",directoryOut="ComboMeshesv20",directoryRef="ComboMeshesv20",pickleOut=meshNames+"_All.pickle")

        # plot combo mesh
        plotAllMeshes(meshNames,methodVal=(250,1.0),nInterpGrid=32,directory="ComboMeshesv20")
        
        # compare individual images to the _All, do all Zernikes
        # use the pickle from comboMeshes which are post-adjustment
        # plot individual files
        
        meshName1 = meshNames + "_All"

        # use nInterpGrid=3 here to match 6x3 bins/CCD for each image
        da1 = mkDonutAna(meshName1,"ScienceOnly","bmedian",directory="ComboMeshesv20",nInterpGrid=3)

        comboPickle = pickle.load(open("ComboMeshesv20/" + meshNames + "_All.pickle"))

        for i in imageList:

            meshDict = comboPickle[i]
            meshName2 = "%s_%d" % (meshNames+i)

            iZs = [4,5,6,7,8,9,10,11]
            minDict = {4:-20.,5:-0.2,6:-0.2,7:-0.075,8:-0.075,9:-0.05,10:-0.05,11:-0.05,14:-0.05,15:-0.05}
            maxDict = {4:20.,5:0.2,6:0.2,7:0.075,8:0.075,9:0.05,10:0.05,11:0.05,14:0.05,15:0.05}
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
                        # remake as a 6x3 grip in each CCD
                        grid = oldGridArray[icoord]
                        grid[0] = 6
                        grid[3] = 3
                        gridArray[icoord] = grid
                    aNewMesh = PointMesh(coordList,gridArray,aMesh.pointsArray,myMethod='bmedian')
            
                    diffMesh = aNewMesh.subtractMesh(da1.meshDict[meshName])
        
                    f,a,c = diffMesh.plotMeshMPL2D(zmin=minDict[iZ],zmax=maxDict[iZ],cmap=cm.jet)
                    f.savefig(meshName1+"-minus-"+meshName2+"_bmedian_z%d.png" % (iZ))
                    plt.close()

                    
def makeIndividualDifferenceMeshes(imgList="284601:284698",meshName1="Science-20140212s2-v22i1_All",meshDirectory="ComboMeshesv22",meshPrefix="Science-20140212s2-v22i1",zernList=None):
    """ compare individual images to the _All,
    use the pickle from comboMeshes which are post-adjustment use BMedian now as the method!
    """

    plt.interactive(False)
    if zernList == None:
        zernList = [4,5,6,7,8,9,10]
        
    # probably should use nInterpGrid=3 here to match 6x3 bins/CCD for each image
    da1 = mkDonutAna(meshName1,"ScienceOnly","bmedian",directory=meshDirectory,nInterpGrid=3)

    # this comboPickle has the points for each Image after the fits to the reference mesh
    comboPickle = pickle.load(open("%s/%s.pickle" % (meshDirectory,meshName1),"rb"))

    imageList = decodeNumberList(imgList)

    for i in imageList:
        print("Image ",i)

        meshDict = comboPickle[i]
        meshName2 = "%s_%d" % (meshPrefix,i)

        iZs = zernList
        minDict = {4:-20.,5:-0.2,6:-0.2,7:-0.075,8:-0.075,9:-0.05,10:-0.05,11:-0.05,14:-0.05,15:-0.05}
        maxDict = {4:20.,5:0.2,6:0.2,7:0.075,8:0.075,9:0.05,10:0.05,11:0.05,14:0.05,15:0.05}
        diffMeshDict = {}

        # dictionary for output dataframe Medianed
        dfDict = OrderedDict()
        # dictionary for individual Donut data
        dfDonutDict = OrderedDict()
        firstOnly = True
    
        for iZ in iZs:
            meshName = "z%dMesh" % (iZ)
            if meshName in da1.meshDict and meshName in meshDict:

                # need to rebuild the interpolation for the mesh in the pickle, must start all over from the points!
            
                aMesh = meshDict[meshName]
                x,y,z = aMesh.getXYZpoints()    # these points are the ones fit to the reference
                oldGridArray = aMesh.gridArray  # this stores limits of the CCDs
                coordList = aMesh.coordList
                gridArray = {}
                for icoord in coordList:
                    # remake as a 6x3 grip in each CCD
                    grid = oldGridArray[icoord]  # use the x,y limits of the CCDs, just change the number of bins
                    grid[0] = 6
                    grid[3] = 3
                    gridArray[icoord] = grid

                # this mesh has the Nentry and MAD values,
                aNewMesh = PointMesh(coordList,gridArray,aMesh.pointsArray,myMethod='bmedian')

                # subtract: this image's mesh - Reference mesh
                refMesh = da1.meshDict[meshName]
                diffMesh = aNewMesh.subtractMesh(refMesh)
                diffMeshDict[meshName] = diffMesh
        
                f,a,c = diffMesh.plotMeshMPL2D(zmin=minDict[iZ],zmax=maxDict[iZ],cmap=cm.jet)
                f.savefig(meshName1+"-minus-"+meshName2+"_bmedian_z%d.png" % (iZ))
                plt.close()

                # fill dataFrame lists
                dfX = []
                dfY = []
                dfMedian = []
                dfDelta = []
                dfNentry = []
                dfMAD = []
                for icoord in coordList:
                    x,y = aNewMesh.interpGrids[icoord]
                    xflat = x.flatten()
                    yflat = y.flatten()
                    difference_values = aNewMesh.doInterp(icoord,xflat,yflat) - refMesh.doInterp(icoord,xflat,yflat)  # same as above!
                    median_values = aNewMesh.interpBMedian[icoord]
                    nentry_values = aNewMesh.interpBNentry[icoord]
                    mad_values = aNewMesh.interpBMAD[icoord]

                try:
                    if median_values.shape[0] > 0:
                        dfX.extend(xflat.tolist())
                        dfY.extend(yflat.tolist())
                        dfMedian.extend(median_values.flatten().tolist())
                        dfDelta.extend(difference_values.flatten().tolist())
                        dfNentry.extend(nentry_values.flatten().tolist())
                        dfMAD.extend(mad_values.flatten().tolist())
                except:
                    print("No median values at ",icoord)

            # build dictionary
            if firstOnly:        
                nEach = len(dfX)        
                dfDict['x'] = dfX
                dfDict['y'] = dfY
                dfDict['ifile'] = (i * np.ones((nEach))).tolist()
                firstOnly = False
            
            dfDict['median_z%d' % (iZ)] = dfMedian
            dfDict['delta_z%d' % (iZ)] = dfDelta
            dfDict['nentry_z%d' % (iZ)] = dfNentry
            dfDict['mad_z%d' % (iZ)] = dfMAD
            
    # write out dfDict as a DataFrame to a pickle file 
    df  = pa.DataFrame(dfDict)
    df.to_pickle(meshName1+"-minus-"+meshName2+".pkl",protocol=4)


def plotMeshes(meshName,meshDirectory):
    da = mkDonutAna(meshName,"ScienceOnly","idw",directory=meshDirectory,methodVal=(250,1.0),nInterpGrid=32)
    iZs = [4,5,6,7,8,9,10,11,14,15]
    minDict = {4:-20.,5:-0.2,6:-0.2,7:-0.125,8:-0.125,9:-0.3,10:-0.3,11:-0.15,14:-0.1,15:-0.1}
    maxDict = {4:20.,5:0.2,6:0.2,7:0.125,8:0.125,9:0.3,10:0.3,11:0.15,14:0.1,15:0.1}
    for iZ in iZs:
        meshN = "z%dMesh" % (iZ)
        if meshN in da.meshDict:
            aMesh = da.meshDict[meshN]
            plt.interactive(False)
            f,a,c = aMesh.plotMeshMPL2D(zmin=minDict[iZ],zmax=maxDict[iZ],cmap=cm.jet)
            plt.interactive(False)
            f.savefig(meshName+"_z%d.png" % (iZ))


def plotMeshesDF(dfDirectory,dfName,zVarPattern="z%dcorr",usedFlag=True,minDict={4:-20.,5:-0.2,6:-0.2,7:-0.125,8:-0.125,9:-0.3,10:-0.3,11:-0.15,14:-0.1,15:-0.1},maxDict={4:20.,5:0.2,6:0.2,7:0.125,8:0.125,9:0.3,10:0.3,11:0.15,14:0.1,15:0.1},interactiveFlag=False,iZs=[4,5,6,7,8,9,10,11,14,15]):
    #iZs = [4,5,6,7,8,9,10,11,14,15]
    if usedFlag:
        da = mkDonutAnaDF(dfDirectory+"/"+dfName+".pkl",iZs,zVarPattern=zVarPattern,sensorSet="ScienceOnly",method="idw",methodVal=(250,1.0),nInterpGrid=32,donutCutString="flagFinal==True")
    else:
        da = mkDonutAnaDF(dfDirectory+"/"+dfName+".pkl",iZs,zVarPattern=zVarPattern,sensorSet="ScienceOnly",method="idw",methodVal=(250,1.0),nInterpGrid=32)
        
    for iZ in iZs:
        meshN = "z%dMesh" % (iZ)
        if meshN in da.meshDict:
            aMesh = da.meshDict[meshN]
            plt.interactive(interactiveFlag)
            f,a,c = aMesh.plotMeshMPL2D(zmin=minDict[iZ],zmax=maxDict[iZ],cmap=cm.jet,interactiveFlag=interactiveFlag)
            plt.interactive(interactiveFlag)
            f.savefig(dfName+"_z%d.png" % (iZ))

def plotMeshesDFPub(dfDirectory,dfName,zVarPattern="z%dcorr",usedFlag=True,minDict={4:-20.,5:-0.2,6:-0.2,7:-0.125,8:-0.125,9:-0.3,10:-0.3,11:-0.15,14:-0.1,15:-0.1},maxDict={4:20.,5:0.2,6:0.2,7:0.125,8:0.125,9:0.3,10:0.3,11:0.15,14:0.1,15:0.1},interactiveFlag=False,iZs=[4,5,6,7,8,9,10,11,14,15]):
    #iZs = [4,5,6,7,8,9,10,11,14,15]
    if usedFlag:
        da = mkDonutAnaDF(dfDirectory+"/"+dfName+".pkl",iZs,zVarPattern=zVarPattern,sensorSet="ScienceOnly",method="idw",methodVal=(250,1.0),nInterpGrid=32,donutCutString="flagFinal==True")
    else:
        da = mkDonutAnaDF(dfDirectory+"/"+dfName+".pkl",iZs,zVarPattern=zVarPattern,sensorSet="ScienceOnly",method="idw",methodVal=(250,1.0),nInterpGrid=32)
        
    for iZ in iZs:
        meshN = "z%dMesh" % (iZ)
        if meshN in da.meshDict:
            aMesh = da.meshDict[meshN]
            plt.interactive(interactiveFlag)
            f,a,c = aMesh.plotMeshMPL2DPub(zmin=minDict[iZ],zmax=maxDict[iZ],cmap=cm.jet,interactiveFlag=interactiveFlag)
            plt.interactive(interactiveFlag)
            f.savefig(dfName+"_z%d_Pub.png" % (iZ))


    
def submitMeshJobs(dateid,seqno,imageListStr,iterList=[2],nosnipStr='True'):
    """ do everything in pipeline - not working yet?"""

    iterName = ['','first','second']
    zernListStr = ['','4 5 6 7 8 9 10','4 5 6 7 8 9 10 11 14 15']
    
    # fit the donuts
    #command = "submitDecamDriver.py -d %d -n %s  -c decam-ScienceExtra-v22.cfg -seq %d -ver 22 --noSnip %s " % (dateid,imageListStr,seqno,nosnipStr)
    #print(command)
    #os.system(command)

    waitForBatch(maxcount=300)

    # collect fits into a pickle
    imageList = decodeNumberList(imageListStr)
    for iFile in imageList:
        for i in iterList:
            command = "bsub -W 2:00 -o dump-%ds%d-%d-%s.log dumpFitsHeadtoTable.py -i '$KIPACDISK/Donuts/%ds%d/%d/v22/DECam*.%s.donut.fits.fz'  -o $KIPACDISK/Donuts/%ds%d/%d/v22/DECam_%d.%s"  % (dateid,seqno,iFile,iterName[i],dateid,seqno,iFile,iterName[i],dateid,seqno,iFile,iFile,iterName[i])
            print(command)
            os.system(command)

    waitForBatch()

    # filter Donuts for meshes
    for iFile in imageList:
        for i in iterList:
            command = "bsub -W 0:60 -o mesher-%ds%d-%d-%s.log decamMesher.py -i $KIPACDISK/Donuts/%ds%d/%d/v22/DECam_%d.%s.pkl -o  Science-%ds%d-v22i%d   -odir Meshesv22 -f %d -doScience -neleCut 5.0 -nMADCut 20. -nMADNNCut 50.0 -zernList %s" % (dateid,seqno,iFile,iterName[i],dateid,seqno,iFile,iFile,iterName[i],dateid,seqno,i,iFile,zernListStr[i])
            print(command)
            os.system(command)    

    waitForBatch()

    #
    # combine, referencing to Zemax
    #
    meshNameAr2 =  "decam_2012-nominalzernike-bysensor"
    meshNames = "Science-%ds%d-v22i2" % (dateid,seqno)
    adjustFlag = [False,False,False,False,False,True,True,False,False,False,False,False,False,False,False,False] # only adjust Astigmatism in Angle
    combineMeshes(meshNameAr2,meshNames,imageList,method="idw",methodVal=(20,1.0),directoryIn="Meshesv22",directoryOut="ComboMeshesZemaxv22",directoryRef="MeshesZemax",methodRef='rbf',pickleOut=meshNames+"_All.pickle",adjustAngleAll=adjustFlag)

    # convert donut_summary info to a DataFrame
    df = getDonutSummary("ComboMeshesZemaxv22/Science-%ds%d-v22i2_All.pickle" % (dateid,seqno))
    df.to_pickle("ComboMeshesZemaxv22/Science-%ds%d-v22i2_All.pkl" % (dateid,seqno))

    #
    # combine, iterating for a 2nd time
    #
    meshNameAr2 = "Science-%ds%d-v22i2_All" % (dateid,seqno)
    meshNames = "Science-%ds%d-v22i2"  % (dateid,seqno)
    adjustFlag = [False,False,False,False,False,True,True,False,False,False,False,False,False,False,False,False] # only adjust Astigmatism in Angle
    combineMeshes(meshNameAr2,meshNames,imageList,method="idw",methodVal=(20,1.0),directoryIn="Meshesv22",directoryOut="ComboMeshesZemaxIteration2v22",directoryRef="ComboMeshesZemaxv22",pickleOut=meshNames+"_All.pickle",adjustAngleAll=adjustFlag)

    # convert donut_summary info to a DataFrame
    df = getDonutSummary("ComboMeshesZemaxIteration2v22/Science-%ds%d-v22i2_All.pickle"  % (dateid,seqno))
    df.to_pickle("ComboMeshesIteration2Zemaxv22/Science-%ds%d-v22i2_All.pkl"  % (dateid,seqno))

    # make plots per image
    meshName1="Science-%ds%d-v22i2_All" % (dateid,seqno)
    meshDirectory="ComboMeshesZemaxIteration2v22"
    meshPrefix="Science-%ds%d-v22i2" % (dateid,seqno)
    zernList= [4,5,6,7,8,9,10,11,14,15]
    makeIndividualDifferenceMeshes(imgList=imageListStr,meshName1=meshName1,meshDirectory=meshDirectory,meshPrefix=meshPrefix,zernList=zernList)

    # make plots for All
    plotMeshesDF("ComboMeshesZemaxIteration2v22","ComboMesh_Science-%ds%d-v22i2_All" % (dateid,seqno) )

    # done
    return
