import os,sys
import ROOT
ROOT.gSystem.Load('libRATEvent')
from ROOT.RAT import DSReader

def fillhist( filename, hnpe, hnpe_cos ):

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
        npe = r.GetMC().GetNumPE()

        if cosv>0.0:
            hnpe.Fill( npe )
            hnpe_cos.Fill( cosv, npe )
        #print "event ",i,": pos=",posv.X(),",",posv.Y(),",",posv.Z()," :  dir=",dir.X(),dir.Y(),dir.Z()," : ",cosv


out = ROOT.TFile("output.root","recreate")
hnpe_cos = ROOT.TH2D( "hnpe_v_cos", "", 50, 0, 1.0, 50, 0, 500 )
hnpe = ROOT.TH1D("hnpe","",500, 0, 500)

files = os.listdir( "/net/nudsk0001/d00/scratch/taritree/pingudom/" )
for f in files:
    if ".root" not in f:
        continue
    filename = "/net/nudsk0001/d00/scratch/taritree/pingudom/"+f
    print filename
    fillhist( filename, hnpe, hnpe_cos )

c = ROOT.TCanvas("c","c",800,400)
c.Draw()
hnpe.Draw()
c.Update()

out.Write()
raw_input()



