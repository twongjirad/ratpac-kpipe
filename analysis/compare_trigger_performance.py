import os,sys
from ROOT import *

pairs = {0:{"nodark":"/net/nudsk0001/d00/scratch/taritree/kpipe_ana_nodarknoise_v4/kdar_analysis_nodarknoise_v4_0100.root",
            "widark":"/net/nudsk0001/d00/scratch/taritree/kpipe_ana_wdarknoise_v4_1.6mhz/kdar_analysis_wdarknoise_v4_1.6mhz_0100.root"},
         1:{"nodark":"/net/nudsk0001/d00/scratch/taritree/kpipe_ana_nodarknoise_v4/kdar_analysis_nodarknoise_v4_0020.root",
            "widark":"/net/nudsk0001/d00/scratch/taritree/kpipe_ana_wdarknoise_v4_1.6mhz/kdar_analysis_wdarknoise_v4_1.6mhz_0020.root"},
         2:{"nodark":"kdar_analysis_nodarknoise_v4_0020.root",
            "widark":"kdar_analysis_wdarknoise_v4_1.6mhz_0020.root"},
         3:{"nodark":"test_cr_nodark.root",
            "widark":"test_cr_wdark.root"},
         4:{"nodark":"test_kdar_nodark.root",
            #"widark":"test_kdar_100_wdark_wreject.root"},
            "widark":"test_kdar_wdark_wreject.root"},
         }

pair = pairs[4]

wdark = TChain("mcdata")
wdark.Add( pair["widark"] )
wdark.SetAlias("wd","tend-ttrig")

nodark = TChain("mcdata")
nodark.Add( pair["nodark"] )
nodark.SetAlias("nd.wd","nd.tend-nd.ttrig")


wdark.AddFriend( nodark, 'nd' )

denom     = wdark.GetEntries( "mumomv>0 && rv<150" )
nmistakes_more = wdark.GetEntries( "npulses<nd.npulses && mumomv>0 && rv<150" )
nmistakes_less = wdark.GetEntries( "npulses>nd.npulses && mumomv>0 && rv<150" )
nmistakes_tot = wdark.GetEntries( "npulses!=nd.npulses && mumomv>0 && rv<150" )

# ID
wdark.Scan("predark_idpe:nd.idpe:predark_odpe:nd.odpe:rv:zv:mumomv:npulses:nd.npulses:pulsepe:nd.pulsepe:ttrig:nd.ttrig:wd","npulses!=nd.npulses && predark_idpe>0")

# OD
#wdark.Scan("predark_idpe:nd.idpe:predark_odpe:nd.odpe:rv:zv:mumomv:npulses:nd.npulses:pulsepe:nd.pulsepe:ttrig:nd.ttrig:wd:npulses_veto:nd.npulses_veto:pulse_totodpe:nd.pulse_totodpe",
#"predark_odpe>0 && ( (nd.npulses_veto>0 && npulses_veto==0) || (nd.npulses_veto==0 && npulses_veto>0) )")

if denom>0:
    print "[ Denom ]: ",denom
    print "[ Fracton with more ]: ",nmistakes_more," ",float(nmistakes_more)/denom
    print "[ Fracton with less ]: ",nmistakes_less," ",float(nmistakes_less)/denom
    print "[ Fracton different ]: ",nmistakes_tot," ",float(nmistakes_tot)/denom
