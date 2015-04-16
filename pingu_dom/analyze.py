import os,sys
import ROOT
ROOT.gSystem.Load('libRATEvent')
from ROOT.RAT import DSReader

ds = DSReader("test.root")
nevents = ds.GetTotal()

hnpe = ROOT.TH1D("hnpe","",50, 0, 500)

for i in xrange(0,nevents):
    r = ds.NextEvent()
    hnpe.Fill( r.GetMC().GetNumPE() )
    mc = r.GetMC()
    
    posv = mc.GetMCParticle(0).GetPosition()
    mumomv = mc.GetMCParticle(0).GetMomentum()
    dir = mumomv.Unit()
    print "event ",i,": pos=",posv.X(),",",posv.Y(),",",posv.Z()," :  dir=",dir.X(),dir.Y(),dir.Z()

c = ROOT.TCanvas("c","c",800,400)
c.Draw()
hnpe.Draw()
c.Update()

raw_input()



