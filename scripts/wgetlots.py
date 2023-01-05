#! /usr/bin/env python

import argparse
import os
import pdb
import psycopg2
from collections import OrderedDict
from scriptUtil import decodeNumberList

parser = argparse.ArgumentParser(prog='wgetlots')
parser.add_argument("-d", "--mindate", dest="mindate",type=str,default="20120914",
                  help="mindate")
parser.add_argument("-maxd", "--maxdate", dest="maxdate",type=str,default="20201231",
                  help="maxdate")
parser.add_argument("-info", "--info",
                  dest="info",action="store_true",default=False,
                  help="just get info")


# unpack options
options = parser.parse_args()

start_year = options.mindate[0:4]
start_month = options.mindate[4:6]
start_day = options.mindate[6:8]
end_year = options.maxdate[0:4]
end_month = options.maxdate[4:6]
end_day = options.maxdate[6:8] 

def propCode(propid):
    propnum = int(propid[6:11])
    if propnum<=3:
        propcode = 0
    elif propnum>=9990:
        propcode = 1
    else:
        propcode = 2
    return propcode

# build date/image list
dateNimage = []

# query data base
d = psycopg2.connect("")
c = d.cursor()

sqlComm = "SELECT e.id AS expid,round(10000*extract(year from e.date - interval '15:00:00')+100*extract(month from e.date - interval '15:00:00')+extract(day from e.date - interval '15:00:00')) as dateid, e.exptime,e.propid,e.flavor, COALESCE(NULLIF(replace(e.object,',',' '),''),'noobject') AS object, COALESCE(NULLIF(replace(e.program,',',' '),''),'noprogram') AS program, e.utw_temp, e.focus[3] AS hexdz  FROM exposure.exposure AS e WHERE (e.date BETWEEN '%s-%s-%s 20:00:00' AND '%s-%s-%s 10:00:00') AND e.exptime > 0 AND e.flavor = 'object'  ORDER BY e.id;" % (start_year,start_month,start_day,end_year,end_month,end_day)

c.execute(sqlComm)
dbinfo = c.fetchall()
for image in dbinfo:
    try:
        propcode = propCode(image[3])
        if propcode<=1:
            if options.info:
                print(image)
            dateNimage.append((image[0],int(image[1])))
    except:
        print(image)

# write out date and image range...
datedict = OrderedDict()
for expid,dateid in dateNimage:
    if dateid in datedict :
        expMin,expMax = datedict[dateid]
        datedict[dateid] = (min(expid,expMin),max(expid,expMax))
    else:
        datedict[dateid] = (expid,expid)

for key in datedict:
    print(("-d %d -n '%d:%d'" % (key,datedict[key][0],datedict[key][1]))) 

# loop over these guys
for expid,dateid in dateNimage:
    

    # get on the directory, make it if it isn't there
    #
    dataDirectory = "/nfs/slac/g/ki/ki06/roodman/desdata/%d" % (dateid)
    dataDirectoryExp = os.path.expandvars(dataDirectory)

    # make directory if it doesn't exist
    if not os.path.exists(dataDirectoryExp):
        os.mkdir(dataDirectoryExp)

    # move there!
    os.chdir(dataDirectoryExp)

    # DTS area
    dtsName = "DTS"

    # Dec 4, 2013 - not using port :7443 anymore


    # find correct location at NCSA
    if dateid<20140910:
        weblocation = "https://desar2.cosmology.illinois.edu/DESFiles/desardata/DTS/src/%d/src" % (dateid)
    else:
        weblocation = "https://desar2.cosmology.illinois.edu/DESFiles/desarchive/DTS/raw/%d" % (dateid)


    command = "wget --no-check-certificate --http-user=roodman --http-password=FITdonuts123 -nc -nd -nH -r -k -p -np -nv --cut-dirs=3 %s/DECam_00%d.fits.fz" % (weblocation,expid)
    if options.info:
        #print command
        a = 1
    else:
        os.system(command)
