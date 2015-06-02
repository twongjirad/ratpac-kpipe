import os,sys
from ROOT import *

wdark = TChain("mcdata")
wdark.Add( "test_kdar_wdark_wreject.root" )

nodark = TChain("mcdata")
nodark.Add( "test_kdar_nodark.root" )


wdark.AddFriend( nodark, 'nd' )


wdark.Scan("predark_idpe:nd.idpe:predark_odpe:nd.odpe:rv:zv:mumomv:npulses:nd.npulses:pulsepe:nd.pulsepe:ttrig:nd.ttrig","npulses!=nd.npulses && mumomv>0 && rv<150")


#nentries = wdark.GetEntries()
#for i in xrange(0,entries):
#    wdark.GetEntry(i)
#    nodark.GetEntry(i)
#    print 
