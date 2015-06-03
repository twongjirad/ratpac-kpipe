import os,sys
from ROOT import *

pairs = {0:{"nodark":"/net/nudsk0001/d00/scratch/taritree/kpipe_ana_nodarknoise_v4/kdar_analysis_nodarknoise_v4_0100.root",
            "widark":"/net/nudsk0001/d00/scratch/taritree/kpipe_ana_wdarknoise_v4_1.6mhz/kdar_analysis_wdarknoise_v4_1.6mhz_0100.root"},
         1:{"nodark":"/net/nudsk0001/d00/scratch/taritree/kpipe_ana_nodarknoise_v4/kdar_analysis_nodarknoise_v4_0020.root",
            "widark":"/net/nudsk0001/d00/scratch/taritree/kpipe_ana_wdarknoise_v4_1.6mhz/kdar_analysis_wdarknoise_v4_1.6mhz_0020.root"},
         }

pair = pairs[1]

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

wdark.Scan("predark_idpe:nd.idpe:predark_odpe:nd.odpe:rv:zv:mumomv:npulses:nd.npulses:pulsepe:nd.pulsepe:ttrig:nd.ttrig:wd:nd.wd","npulses!=nd.npulses && mumomv>0 && rv<150 && npulses==1 && nd.npulses>1")

print "[ Denom ]: ",denom
print "[ Fracton with more ]: ",nmistakes_more," ",float(nmistakes_more)/denom
print "[ Fracton with less ]: ",nmistakes_less," ",float(nmistakes_less)/denom
print "[ Fracton different ]: ",nmistakes_tot," ",float(nmistakes_tot)/denom
