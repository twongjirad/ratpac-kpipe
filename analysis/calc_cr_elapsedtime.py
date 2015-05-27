import os,sys
from ROOT import *
# must match with cygen file to have truth info for analysis

pair_list = {}
elaspsed_time_tot = 0
for p in xrange(0,2000):
    cryfile = "/net/nudsk0001/d00/scratch/taritree/cry_gen/0/cry_events_%d.root"%(p)
    #anafile = "/net/nudsk0001/d00/scratch/taritree/crkpipe_ana_wdarknoise_v2/output_cr_analysis_%04d.root"%(p)
    anafile = "/net/nudsk0001/d00/scratch/taritree/crkpipe_ana_wdarknoise/output_cr_analysis_%04d.root"%(p)
    #anafile = "/net/nudsk0001/d00/scratch/taritree/crkpipe_ana_nodarknoise/output_cr_analysis_%04d.root"%(p)
    #anafile = "/net/nudsk0001/d00/scratch/taritree/crkpipe_ana_wdarknoise_v3/output_cr_analysis_%04d.root"%(p)
    rcryfile = TFile( cryfile )
    ranafile = TFile( anafile )
    
    mcdata = ranafile.Get("mcdata")
    nevents = 0
    try:
        nevents = mcdata.GetEntries()
    except:
        nevents = 0

    if nevents>=100:
        cry = rcryfile.Get( "crytree" )
        cry.GetEntry(nevents-1)
        dt = cry.telapsed_sec[0]
        elaspsed_time_tot += dt
        pair_list[p] = { "cryfile":cryfile, "anafile":anafile, "dt":dt }
    rcryfile.Close()
    ranafile.Close()

#print pair_list
print "COMPLETE PAIRS: ",len(pair_list)
print "ELAPSED TIME: ",elaspsed_time_tot,"secs"

