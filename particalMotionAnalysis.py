# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 15:25:14 2015

@author: ISTI_EW
Perform partical motion analysis
"""
import matplotlib.gridspec as gridspec
import obspy, numpy as np, scipy, matplotlib.pyplot as plt, pandas as pd, glob, os
from obspy.signal.polarization import particle_motion_odr as pma
from mpl_toolkits.mplot3d import Axes3D

## Parameters
picks=pd.read_pickle('EventPicks.pkl')
evedir='EventWaveForms'
filt=[1,12,4,True]
#filt=None
tfp=[0,.75] #range of trace to plot from first arrivial in seconds  
ttp=[-.5,2] #range from first arrivial (pick) to plot

if tfp[0]<ttp[0] or tfp[1]>ttp[1]: raise ValueError('ttp[0] must be less than tfp[0] and ttp[1] must be greater than tfp[1]')


temkey=pd.read_csv('TemplateKey.csv') #read in template key
eventsInDir=[os.path.basename(x) for x in glob.glob(os.path.join(evedir,'*'))]

events=list(set(eventsInDir) & set(temkey.NAME)) #get intersection of events in template key and those found in directory

for event in events:
    print (event)
    print('----------')
    TR=obspy.read(os.path.join(evedir,event,'*'))
    TR.detrend('linear')
    if isinstance(filt,list):
        TR.filter('bandpass',freqmin=filt[0],freqmax=filt[1],corners=filt[2],zerophase=filt[3])
    
    pick=picks[(picks.Station==TR[0].stats.network+'.'+TR[0].stats.station)&(picks.Event==os.path.basename(event))]
    if len(pick)<4:
        print ('%s does not have all required phase picks (Make sure to pick P,Pend,S,Send)'%event )
    else:
        trimp=obspy.UTCDateTime(pick[pick.Phase=='P'].iloc[0].TimeStamp)
        trimpend=obspy.UTCDateTime(pick[pick.Phase=='Pend'].iloc[0].TimeStamp)
        trims=obspy.UTCDateTime(pick[pick.Phase=='S'].iloc[0].TimeStamp)
        trimsend=obspy.UTCDateTime(pick[pick.Phase=='Send'].iloc[0].TimeStamp)
        
        #get trace cooresponding to P phase
        trp=TR.slice(starttime=trimp,endtime=trimpend) 
#        if len(trp)<1:
#            continue
        Ep=trp.select(component='E')
        Np=trp.select(component='N')
        Zp=trp.select(component='Z')
        
        #Get trace cooresponding to S phase
        trs=TR.slice(starttime=trims,endtime=trimsend)
        Es=trs.select(component='E')
        Ns=trs.select(component='N')
        Zs=trs.select(component='Z')        
        
#        elim=[min(E[0]),max(E[0])]
#        nlim=[min(N[0]),max(N[0])]
#        zlim=[min(Z[0]),max(Z[0])]       
#        lims=[min([elim[0],nlim[0],zlim[0]]),max([elim[1],nlim[1],zlim[1]])]
        
        
        ## get slice for total waveform visualization
        trp=TR.slice(starttime=trimp+ttp[0],endtime=trimp+ttp[1]) 
        Etp=trp.select(component='E')
        Ntp=trp.select(component='N')
        Ztp=trp.select(component='Z')
        elimtp=[min(Etp[0]),max(Etp[0])]
        nlimtp=[min(Ntp[0]),max(Ntp[0])]
        zlimtp=[min(Ztp[0]),max(Ztp[0])]       
        limstp=[min([elimtp[0],nlimtp[0],zlimtp[0]]),max([elimtp[1],nlimtp[1],zlimtp[1]])] #absolulte limits
        
        
        Psta=(trimp-trp[0].stats.starttime)
        Pend=(trimpend-trp[0].stats.starttime)    
        Ssta=(trims-trp[0].stats.starttime)
        Send=(trimsend-trp[0].stats.starttime)
        
        vlinestart=(tfp[0]-ttp[0]) #value in samples for drawling line showing where pick was made
        vlinestop=(tfp[1]-ttp[0])
        
        t=obspy.Stream([Z[0],N[0],E[0]])
        out=pma(t)
        
        
        ## Initialize figures and plot
        fig=plt.figure(figsize=[10,10])
        gs=gridspec.GridSpec(6,6)
        
        ax1=plt.subplot(gs[:3,:3])
        ax1.plot(Ep,Np,'b.')
        ax1.plot(Es,Ns,'r.')
        ax1.set_xlabel('North')
        ax1.xaxis.set_label_position('top')
        ax1.set_ylabel('East')
        ax1.yaxis.set_label_position('right')
        xlim=plt.xlim()
        ax1.set_yticks([])
        ax1.set_xticks([])
        ax1.set_xticks
        ax1.set_ylim(limstp)
        ax1.set_xlim(limstp)
        ax1.axvline(c='k')
        ax1.axhline(c='k')
        
        ax2=plt.subplot(gs[3:6,:3])
        ax2.plot(Ep,Zp,'b.')
        ax2.plot(Es,Zs,'r.')
        ax2.set_ylabel('Z')
        ax2.set_yticks([])
        ax2.set_xticks([])
        ax2.set_ylim(limstp)
        ax2.set_xlim(limstp)
        
        ax3=plt.subplot(gs[:3,3:6])
        ax3.plot(Zs,Ns,'r.')
        ax3.plot(Zp,Np,'b.')
        ax3.set_xlabel('Z')
        ax3.xaxis.set_label_position('top')
        ax3.set_yticks([])
        ax3.set_xticks([])
        ax3.set_ylim(limstp)
        ax3.set_xlim(limstp)
        
        #plot traces
                
        
        ax4=plt.subplot(gs[3:4,3:6])
        plt.plot(np.linspace(0,len(Ztp[0])/Ztp[0].stats.sampling_rate,len(Ztp[0])),Ztp[0],'k')
        ax4.set_title('Traces')
        ax4.set_ylabel('Z',labelpad=-.25)
        ax4.yaxis.tick_right()
        ax4.tick_params(labelbottom='off')   
        ax4.set_ylim(limstp)
        ax4.set_yticks([])
        ax4.axvline(Psta,c='b')
        ax4.axvline(Pend,c='b')
        ax4.axvline(Ssta,c='r')
        ax4.axvline(Send,c='r')
 
        ax4.set_xlim([0,len(Ztp[0])/Ztp[0].stats.sampling_rate])
        
        ax5=plt.subplot(gs[4:5,3:6])
        ax5.plot(np.linspace(0,len(Ntp[0])/Ntp[0].stats.sampling_rate,len(Ntp[0])),Ntp[0],'k')
        ax5.set_ylabel('North',labelpad=-.25)
        ax5.yaxis.tick_right()
        ax5.tick_params(labelbottom='off')  
        ax5.set_ylim(limstp)        
        ax5.set_yticks([])
        ax5.axvline(Psta,c='b')
        ax5.axvline(Pend,c='b')
        ax5.axvline(Ssta,c='r')
        ax5.axvline(Send,c='r')
        ax5.set_xlim([0,len(Ntp[0])/Ntp[0].stats.sampling_rate])
        
        ax6=plt.subplot(gs[5:6,3:6])
        ax6.plot(np.linspace(0,len(Etp[0])/Etp[0].stats.sampling_rate,len(Etp[0])),Etp[0],'k')
        ax6.set_ylabel('East',labelpad=-.25)
        ax6.yaxis.tick_right()
        ax6.set_xlabel('Time (seconds)')
        ax6.set_ylim(limstp)
        ax6.set_yticks([])
        ax6.axvline(Psta,c='b')
        ax6.axvline(Pend,c='b')
        ax6.axvline(Ssta,c='r')
        ax6.axvline(Send,c='r')
        ax6.set_xlim([0,len(Etp[0])/Etp[0].stats.sampling_rate])
        
        #tr.plot(fig=fig)
        
#        ax.scatter(E,N,Z)
#        plt.plot(E[0].data,N[0].data,'.')
        #plt.title(str(out))
        fig.savefig(event+'.pdf')
        plt.show()
