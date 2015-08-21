# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 11:38:56 2015

@author: ISTI_EW

Script to convert the Scream! output directories to mseed directories (more amenable to obspy)
and to create specific event files
"""
import glob, os, obspy, subprocess, sys,shutil, numpy as np
convertFiles=0 #if true convert files
harvestEvents=1 #If true harvest events

###Convert files options
target='data2' #target directory or file
outdir='mseedata' #out directory 
utilPath=os.path.join(os.getcwd(),'gcf2msd.exe') # path to the gcf2msd convert utility

### Convert Files Functions    
def _getSubDirs(target):
    content=glob.glob(os.path.join(target,'*'))
    return [name for name in content if os.path.isdir(name)]

def _gcf2mseed(infile): #call the guralp utility to create temp file, read as obspy stream
    os.system(utilPath+' '+infile+' '+'/o:'+os.path.join(os.getcwd(),'Temp'))
    tr=obspy.read(os.path.join(os.path.join(os.getcwd(),'Temp','*')))
    shutil.rmtree('Temp')
    return tr

###Convert Files block
if convertFiles:
    if not os.path.exists(target): raise Exception('%s does not exists'%target)
    if not os.path.exists(outdir): os.makedirs(outdir)
    subdirs=_getSubDirs(target)
    subTimes=set()
    for subdir in subdirs:
        sdfiles=(glob.glob(os.path.join(subdir,'*')))
        subTimes.update([os.path.basename(x).split('.')[0][:-1] for x in sdfiles])
        
    for stime in subTimes:
        tr=obspy.Stream()
        for subdir in subdirs:
            files=glob.glob(os.path.join(subdir,stime+'*'))
            if len(files)>0:
                tr+=_gcf2mseed(files[0])
        tr.write(os.path.join(outdir,stime+'.msd'),'mseed')
    
### Harvest Events options

events=['2015-07-13T00-54-00','2015-07-13T01-14-00','2015-07-13T02-23-00','2015-07-13T03-04-00','2015-07-13T04-09-00',
        '2015-07-13T04-10-00','2015-07-13T11-11-00','2015-07-13T11-21-00','2015-07-13T16-42-00','2015-07-13T19-32-00',
        '2015-07-13T19-39-00','2015-07-13T20-58-00','2015-07-13T21-18-00','2015-07-14T07-17-00','2015-07-14T10-35-00']
eventDir='EventWaveForms'
condir=outdir #directory where files currently are
timeBefore=60 #time before event in seconds to grab
timeAfter=240 #Time after


### Harvest Events Functions
def _loadTime(time,UTCS,WFfiles): # load continuous data during desired time frame
    dtime=obspy.UTCDateTime(time).timestamp
    difs=dtime-UTCS
    ind=abs(difs).argmin()

    if ind>0 and ind<len(WFfiles)-1:
        toload=WFfiles[ind-1:ind+2]
    elif ind==0:
        toload=WFfiles[ind:ind+2]
    else:
        toload=WFfiles[ind-1:ind]
    return toload

def _gleamUTC(fil):
    bn=os.path.basename(fil)
    ds=bn.split('.')[0]
    year=int(ds[:4])
    mon=int(ds[4:6])
    day=int(ds[6:8])
    hour=int(ds[9:11])
    minu=int(ds[11:13])
    utc=obspy.UTCDateTime(year=year,month=mon,day=day,hour=hour,minute=minu)
    return utc.timestamp
    
### Harvest Block
if harvestEvents:
    if not os.path.exists(eventDir): os.makedirs(eventDir)
    UTCS=np.array([_gleamUTC(x) for x in glob.glob(os.path.join(condir,'*'))]) 
    WFfiles=[x for x in glob.glob(os.path.join(condir,'*'))]
    for event in events:
        ltimes=_loadTime(event,UTCS,WFfiles)
        tr=obspy.Stream()
        for lt in ltimes:
            tr+=obspy.read(lt)
        utc=obspy.UTCDateTime(event)
        utc1=utc-timeBefore
        utc2=utc+timeAfter
        tr.trim(starttime=utc1,endtime=utc2)
        tr.merge()
        for t in tr:
            t.stats.network='NF'
        if not os.path.exists(os.path.join(eventDir,event)): os.makedirs(os.path.join(eventDir,event))
        tr.write(os.path.join(eventDir,event,'NF.'+tr[0].stats.station+'.'+event+'.pkl'),'pickle')
    
    
    
    
    










    
