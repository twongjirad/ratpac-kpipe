import os,sys
import ROOT
ROOT.gSystem.Load('libRATEvent')
from ROOT.RAT import DSReader

def fillhist( filename, hnpe ):

    ds = DSReader(filename)
    nevents = ds.GetTotal()

    for i in xrange(0,nevents):
        r = ds.NextEvent()
        mc = r.GetMC()

        posv = mc.GetMCParticle(0).GetPosition()
        mumomv = mc.GetMCParticle(0).GetMomentum()
        dir = mumomv.Unit()

        norm = -1.0*posv.Unit()
        cosv = norm.X()*dir.X() + norm.Y()*dir.Y() + norm.Z()*dir.Z()

        if cosv>0.0:
            hnpe.Fill( r.GetMC().GetNumPE() )
        print "event ",i,": pos=",posv.X(),",",posv.Y(),",",posv.Z()," :  dir=",dir.X(),dir.Y(),dir.Z()," : ",cosv

filename = "test.root"
hnpe = ROOT.TH1D("hnpe","",50, 0, 500)
fillhist( filename, hnpe )

c = ROOT.TCanvas("c","c",800,400)
c.Draw()
hnpe.Draw()
c.Update()

raw_input()



