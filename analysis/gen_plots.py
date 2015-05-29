import os,sys
from ROOT import *
import array

gStyle.SetOptStat(0)

#folder="kdar_wdark"
folder="kdar_nodark"
os.system("mkdir -p figs/%s/eps"%(folder))

if folder=="kdar_nodark":
    tf = TFile("rootfiles/run100_200_analysis.root")
else:
    tf = TFile("rootfiles/kpipeana_darknoise_300_399.root") # darknoise
mcdata = tf.Get("mcdata")

pe_scale = 4500.0/8500.0
mcdata.SetAlias("pescale","4500.0/8500.0")
mcdata.SetAlias("standard_cuts","pulsepe[0]*pescale>10.0 && pulsepe[1]*pescale>5.0 && abs(pulsez[0]-pulsez[1])<200.0 && pulsepe[0]*pescale>100 && pulsepe[0]*pescale<1300 && pulsepe[1]*pescale<800.0")
mcdata.SetAlias("nobounds_cuts","pulsepe[0]*pescale>10.0 && pulsepe[1]*pescale>5.0 && abs(pulsez[0]-pulsez[1])<200.0")
mcdata.SetAlias("nood_cuts","pulsepe[0]*pescale>10.0 && pulsepe[1]*pescale>5.0 && abs(pulsez[0]-pulsez[1])<200.0 && pulsepe[0]*pescale>600 && pulsepe[0]*pescale<1300 && pulsepe[1]*pescale<800.0")

denom = mcdata.GetEntries("mumomv>0 && rv<150 && abs(zv)<4500")
n_standard = mcdata.GetEntries("standard_cuts && (mumomv>0 && rv<150 && abs(zv)<4500)")
print "Denom: ",denom
print "Standard: ",n_standard," (",float(n_standard)/float(denom),")"

tout = TFile("out.root","RECREATE")

# pe versus total ke
c1 = TCanvas("c1","",800,600)
c1.Draw()
c1.cd().SetLeftMargin(0.15)
c1.cd().SetRightMargin(0.18)
c1.cd().SetTopMargin(0.05)
c1.cd().SetBottomMargin(0.15)
h2_energy_scale = TH2D("h2_energy_scale",";Visible Energy, E_{vis} (MeV);pe in first pulse;fraction of events",50,0,250,50,0,1500)
h2_energy_scale.GetXaxis().SetLabelSize(0.05)
h2_energy_scale.GetXaxis().SetTitleSize(0.06)
h2_energy_scale.GetXaxis().SetTitleOffset(1.1)
h2_energy_scale.GetYaxis().SetLabelSize(0.05)
h2_energy_scale.GetYaxis().SetTitleSize(0.06)
h2_energy_scale.GetYaxis().SetTitleOffset(1.1)
h2_energy_scale.GetZaxis().SetTitleSize(0.05)
mcdata.Draw("pescale*pulsepe[0]:(mumomv+totkeprotonv)>>h2_energy_scale","(mumomv>0 && nobounds_cuts)","COLZ")
maxbx  = array.array('i',[0])
maxby = array.array('i',[0])
maxbz = array.array('i',[0])
h2_energy_scale.GetMaximumBin( maxbx, maxby, maxbz )
maxEvis = h2_energy_scale.GetXaxis().GetBinCenter( maxbx[0] )
maxpe   = h2_energy_scale.GetYaxis().GetBinCenter( maxby[0] )
escale = maxEvis/maxpe
print maxEvis,maxpe,escale,"(MeV/pe)"
print "escale bounds interaction: ",600.0*escale," to ",1300*escale," MeV"
print "escale bounds Michel: ",5.0*escale," to ",800*escale," MeV"
h2_energy_scale.Scale( 1.0/h2_energy_scale.Integral() )
c1.Update()
c1.SaveAs("figs/%s/h2_energy_scale.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h2_energy_scale.eps"%(folder))

# FIRST PULSE PE
c1.cd().SetRightMargin(0.05)
c1.cd().SetLeftMargin(0.15)
h_idpe = TH1D("h_idpe",";pe in prompt pulse; fraction of events",50,0,2000)
h_idpe_zoom = TH1D("h_idpe_zoom",";pe in prompt pulse; fraction of events",50,0,500)
for h in [h_idpe, h_idpe_zoom]:
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(1.2)
#h_idpe.GetZaxis().SetTitleSize(0.05)
mcdata.Draw("pescale*pulsepe[0]>>h_idpe","(mumomv>0 && rv<150 && abs(zv)<4500)","COLZ")
mcdata.Draw("pescale*pulsepe[0]>>h_idpe_zoom","(mumomv>0 && rv<150 && abs(zv)<4500)","COLZ")
idpe_integral = h_idpe.Integral()
h_idpe.Scale( 1.0/idpe_integral )
h_idpe_zoom.Scale( 1.0/idpe_integral )
c1.Update()
h_idpe.Draw()
c1.SaveAs("figs/%s/h_idpe.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_idpe.eps"%(folder))
h_idpe_zoom.Draw()
c1.SaveAs("figs/%s/h_idpe_zoom.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_idpe_zoom.eps"%(folder))


# Michel Pulse PE
c1.cd().SetRightMargin(0.05)
c1.cd().SetLeftMargin(0.15)
h_michel_scale = TH1D("h_michel_scale",";pe in Michel pulse; fraction of events",80,0,800)
h_michel_scale_zoom = TH1D("h_michel_scale_zoom",";pe in Michel pulse; fraction of events",50,0,200)
for h in [ h_michel_scale, h_michel_scale_zoom ]:
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(1.2)
#h_michel_scale.GetZaxis().SetTitleSize(0.05)
mcdata.Draw("pescale*pulsepe[1]>>h_michel_scale","(mumomv>0 && rv<150 && abs(zv)<4500) && npulses>1","COLZ")
mcdata.Draw("pescale*pulsepe[1]>>h_michel_scale_zoom","(mumomv>0 && rv<150 && abs(zv)<4500) && npulses>1","COLZ")
michel_integral = h_michel_scale.Integral()
h_michel_scale.Scale( 1.0/michel_integral )
h_michel_scale_zoom.Scale( 1.0/michel_integral )
c1.Update()
h_michel_scale.Draw()
c1.SaveAs("figs/%s/h_michel_scale.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_michel_scale.eps"%(folder))
h_michel_scale_zoom.Draw()
c1.SaveAs("figs/%s/h_michel_scale_zoom.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_michel_scale_zoom.eps"%(folder))

# OD PE
c1.cd().SetRightMargin(0.05)
c1.cd().SetLeftMargin(0.15)
h_odpe = TH1D("h_odpe",";pe in Michel pulse; fraction of events",50,0,100)
h_odpe.GetXaxis().SetLabelSize(0.05)
h_odpe.GetXaxis().SetTitleSize(0.06)
h_odpe.GetXaxis().SetTitleOffset(1.2)
h_odpe.GetYaxis().SetLabelSize(0.05)
h_odpe.GetYaxis().SetTitleSize(0.06)
h_odpe.GetYaxis().SetTitleOffset(1.2)
#h_odpe.GetZaxis().SetTitleSize(0.05)
if folder=="nodark":
    mcdata.Draw("pescale*odpe>>h_odpe","(mumomv>0 && rv<150 && abs(zv)<4500) && odpe>0","COLZ")
else:
    mcdata.Draw("pescale*pulse_totodpe>>h_odpe","(mumomv>0 && rv<150 && abs(zv)<4500) && pulse_totodpe>0","COLZ")
h_odpe.Scale( 1.0/h_odpe.Integral() )
c1.Update()
c1.SaveAs("figs/%s/h_odpe.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_odpe.eps"%(folder))

raw_input()
