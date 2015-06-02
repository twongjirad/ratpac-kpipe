import os,sys
from ROOT import *

wdark = TChain("mcdata")
wdark.Add( "test_kdar_wdark_wreject.root" )
wdark.SetAlias("wd","tend-ttrig")

nodark = TChain("mcdata")
nodark.Add( "test_kdar_nodark.root" )


wdark.AddFriend( nodark, 'nd' )

denom     = wdark.GetEntries( "mumomv>0 && rv<150" )
nmistakes_more = wdark.GetEntries( "npulses<nd.npulses && mumomv>0 && rv<150" )
nmistakes_less = wdark.GetEntries( "npulses>nd.npulses && mumomv>0 && rv<150" )
nmistakes_tot = wdark.GetEntries( "npulses!=nd.npulses && mumomv>0 && rv<150" )

wdark.Scan("predark_idpe:nd.idpe:predark_odpe:nd.odpe:rv:zv:mumomv:npulses:nd.npulses:pulsepe:nd.pulsepe:ttrig:nd.ttrig:wd","npulses!=nd.npulses && mumomv>0 && rv<150")

print "[ Denom ]: ",denom
print "[ Fracton with more ]: ",nmistakes_more," ",float(nmistakes_more)/denom
print "[ Fracton with less ]: ",nmistakes_less," ",float(nmistakes_less)/denom
print "[ Fracton different ]: ",nmistakes_tot," ",float(nmistakes_tot)/denom
