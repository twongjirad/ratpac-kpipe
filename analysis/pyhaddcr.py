import os,sys
from ROOT import *
# must match with cygen file to have truth info for analysis

pair_list = {}
os.system("rm temp_cranalysis_f*.root")

cry_merge = "crygen_merged.root"
ana_merge = "crana_merged_wdarknoise_0_2000_v3.root"
for p in xrange(0,2000):
    cryfile = "/net/nudsk0001/d00/scratch/taritree/cry_gen/0/cry_events_%d.root"%(p)
    #anafile = "/net/nudsk0001/d00/scratch/taritree/cr_trg_out/output_cr_analysis_%d.root"%(p)
    #anafile = "/net/nudsk0001/d00/scratch/taritree/crkpipe_ana_nodarknoise/output_cr_analysis_%04d.root"%(p)
    #anafile = "/net/nudsk0001/d00/scratch/taritree/crkpipe_ana_nodarknoise_v2/output_cr_analysis_%04d.root"%(p)
    #anafile = "/net/nudsk0001/d00/scratch/taritree/crkpipe_ana_wdarknoise/output_cr_analysis_%04d.root"%(p)
    #anafile = "/net/nudsk0001/d00/scratch/taritree/crkpipe_ana_wdarknoise_v2/output_cr_analysis_%04d.root"%(p)
    anafile = "/net/nudsk0001/d00/scratch/taritree/crkpipe_ana_wdarknoise_v3/output_cr_analysis_%04d.root"%(p)
    rcryfile = TFile( cryfile )
    ranafile = TFile( anafile )
    
    mcdata = ranafile.Get("mcdata")
    nevents = 0
    try:
        nevents = mcdata.GetEntries()
    except:
        nevents = 0

    if nevents>=100:
        pair_list[p] = { "cryfile":cryfile, "anafile":anafile }
    rcryfile.Close()
    ranafile.Close()

print "COMPLETE PAIRS: ",len(pair_list)

cry_add = ""
ana_add = ""
ntemp_files = len(pair_list)/100
if len(pair_list)%100!=0:
    ntemp_files+=1
print "Number of temp files: ",ntemp_files
for p in xrange(0,ntemp_files):
    ana_temp = "temp_cranalysis_f%d.root"%(p)
    cry_temp = "temp_crygen_f%d.root"%(p)
    cry_add += " "+cry_temp
    ana_add += " "+ana_temp
    cry_addlist = ""
    ana_addlist = ""
    for n in xrange(0,100):
        fnum = 100*p + n
        if fnum in pair_list:
            g = pair_list[fnum]["cryfile"]
            f = pair_list[fnum]["anafile"]
            if os.path.exists(f)==True and os.path.exists(g)==True:
                ana_addlist += " "+f
                cry_addlist += " "+g
    os.system("hadd %s %s"%(ana_temp,ana_addlist))
    #os.system("hadd %s %s"%(cry_temp,cry_addlist))

os.system("rm %s"%(ana_merge))
#os.system( "hadd %s %s"%(cry_merge,cry_add))
os.system( "hadd %s %s"%(ana_merge,ana_add))
os.system("rm temp_cranalysis_f*.root")
