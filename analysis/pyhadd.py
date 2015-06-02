import os,sys
from ROOT import *
# must match with cygen file to have truth info for analysis

pair_list = {}
os.system("rm temp_analysis_f*.root")

NEVENTS_EXPECTED = 1000
ana_merge = "kdarana_merged_wdarknoise_0_200_v4_1.6mhz_4.0sig_pass2.root"

for p in xrange(0,200):

    #anafile = "/net/nudsk0001/d00/scratch/taritree/kpipe_ana_wdarknoise_v4_1.6mhz/kdar_analysis_wdarknoise_v4_1.6mhz_%04d.root"%(p)
    anafile = "/net/nudsk0001/d00/scratch/taritree/kpipe_ana_wdarknoise_v4_1.6mhz/sig4.0/kdar_analysis_wdarknoise_v4_1.6mhz_%04d.root"%(p)
    #anafile = "/net/nudsk0001/d00/scratch/taritree/kpipe_ana_wdarknoise_v4_1.6mhz/sig3.0/kdar_analysis_wdarknoise_v4_1.6mhz_%04d.root"%(p)

    ranafile = TFile( anafile )
    
    mcdata = ranafile.Get("mcdata")
    nevents = 0
    try:
        nevents = mcdata.GetEntries()
    except:
        nevents = 0

    if nevents>=NEVENTS_EXPECTED:
        pair_list[p] = { "anafile":anafile }
    ranafile.Close()

print "COMPLETE PAIRS: ",len(pair_list)

ana_add = ""
ntemp_files = len(pair_list)/100
if len(pair_list)%100!=0:
    ntemp_files+=1
print "Number of temp files: ",ntemp_files
for p in xrange(0,ntemp_files):
    ana_temp = "temp_analysis_f%d.root"%(p)
    ana_add += " "+ana_temp
    ana_addlist = ""
    for n in xrange(0,100):
        fnum = 100*p + n
        if fnum in pair_list:
            f = pair_list[fnum]["anafile"]
            if os.path.exists(f)==True:
                ana_addlist += " "+f
    os.system("hadd %s %s"%(ana_temp,ana_addlist))

os.system("rm %s"%(ana_merge))
os.system( "hadd %s %s"%(ana_merge,ana_add))
os.system("rm temp_analysis_f*.root")

print "FINISHED: Completed files=",len(pair_list)
