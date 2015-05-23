import os,sys
from ROOT import *
import array

gStyle.SetOptStat(0)

#elasped_time = 1.17161308983*0.5 # seconds
elasped_time = 0.586272565374 # seconds
scale = 1.0/elasped_time
search_window = 250.0e-9 # ns

folder="veto_wdark"
#folder="nodark"
os.system("mkdir -p figs/%s/eps"%(folder))

if folder=="nodark":
    tf = TFile("rootfiles/crana_merged_0_1999.root")
else:
    tf = TFile("crana_merged_wdarknoise_0_1000.root") # darknoise
mcdata = tf.Get("mcdata")

pe_scale = 4500.0/8500.0
mcdata.SetAlias("pescale","4500.0/8500.0")
mcdata.SetAlias("standard_cuts","npulses==2 && pulsepe[1]*pescale>5.0 && abs(pulsez[0]-pulsez[1])<200.0 && pulsepe[0]*pescale>600 && pulsepe[0]*pescale<1300 && pulsepe[1]*pescale<800.0 && npulses_veto==0 && pulse_totodpe<1.0")
mcdata.SetAlias("nood_cuts",    "npulses==2 && pulsepe[1]*pescale>5.0 && abs(pulsez[0]-pulsez[1])<200.0 && pulsepe[0]*pescale>600 && pulsepe[0]*pescale<1300 && pulsepe[1]*pescale<800.0")
mcdata.SetAlias("nobounds_cuts","pulsepe[0]*pescale>10.0 && pulsepe[1]*pescale>5.0 && abs(pulsez[0]-pulsez[1])<200.0")


tout = TFile("out.root","RECREATE")

# CANVAS
c1 = TCanvas("c1","",800,600)
c1.Draw()
c1.cd().SetLeftMargin(0.15)
c1.cd().SetRightMargin(0.18)
c1.cd().SetTopMargin(0.05)
c1.cd().SetBottomMargin(0.15)

# FIRST PULSE PE
c1.cd().SetRightMargin(0.05)
c1.cd().SetLeftMargin(0.15)
h_idpe = TH1D("h_idpe",";pe in prompt pulse; fraction of events",100,0,40000)
h_idpe_zoom = TH1D("h_idpe_zoom",";pe in prompt pulse; fraction of events",100,0,1000)
h_idpe_wcuts = TH1D("h_idpe_wcuts",";pe in prompt pulse; fraction of events",100,0,2000)
for h in [h_idpe,h_idpe_zoom]:
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(1.2)

mcdata.Draw("pescale*pulsepe[0]>>h_idpe","")
c1.Update()
c1.SaveAs("figs/%s/h_idpe.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_idpe.eps"%(folder))

mcdata.Draw("pescale*pulsepe[0]>>h_idpe_wcuts","standard_cuts")
c1.Update()
c1.SaveAs("figs/%s/h_idpe_wcuts.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_idpe_wcuts.eps"%(folder))

mcdata.Draw("pescale*pulsepe[0]>>h_idpe_zoom","")
c1.Update()
c1.SaveAs("figs/%s/h_idpe_zoom.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_idpe_zoom.eps"%(folder))
for h in [h_idpe,h_idpe_zoom,h_idpe_wcuts]:
    h.Scale(scale)

# CR CONTENT
hncr = []
for ncr in ["ncr_muons","ncr_photons","ncr_electrons"]:
    h = TH1D("h_%s"%(ncr),";pe in prompt pulse; fraction of events",10,0,10)
    mcdata.Draw("%s>>h_%s"%(ncr,ncr),"standard_cuts")
    c1.Update()
    c1.SaveAs("figs/%s/h_%s.pdf"%(folder,ncr))


# Michel Pulse PE
c1.cd().SetRightMargin(0.05)
c1.cd().SetLeftMargin(0.15)
h_michel_scale = TH1D("h_michel_scale",";pe in Michel pulse; fraction of events",80,0,800)
h_michel_scale.GetXaxis().SetLabelSize(0.05)
h_michel_scale.GetXaxis().SetTitleSize(0.06)
h_michel_scale.GetXaxis().SetTitleOffset(1.2)
h_michel_scale.GetYaxis().SetLabelSize(0.05)
h_michel_scale.GetYaxis().SetTitleSize(0.06)
h_michel_scale.GetYaxis().SetTitleOffset(1.2)
#h_michel_scale.GetZaxis().SetTitleSize(0.05)
mcdata.Draw("pescale*pulsepe[1]>>h_michel_scale","npulses>1")
h_michel_scale.Scale( 1.0/h_michel_scale.Integral() )
c1.Update()
c1.SaveAs("figs/%s/h_michel_scale.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_michel_scale.eps"%(folder))

# OD PE
c1.cd().SetRightMargin(0.05)
c1.cd().SetLeftMargin(0.15)
h_odpe = TH1D("h_odpe",";pe in veto pulses; pe/bin/sec",50,0,100)
h_odpe.GetXaxis().SetLabelSize(0.05)
h_odpe.GetXaxis().SetTitleSize(0.06)
h_odpe.GetXaxis().SetTitleOffset(1.2)
h_odpe.GetYaxis().SetLabelSize(0.05)
h_odpe.GetYaxis().SetTitleSize(0.06)
h_odpe.GetYaxis().SetTitleOffset(1.2)
#h_odpe.GetZaxis().SetTitleSize(0.05)
if folder=="nodark":
    mcdata.Draw("pescale*odpe>>h_odpe","stadard_cuts && odpe>0")
else:
    mcdata.Draw("pescale*pulse_totodpe>>h_odpe","nood_cuts")
h_odpe.Scale( scale )
c1.Update()
c1.SaveAs("figs/%s/h_odpe.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_odpe.eps"%(folder))

# OD Pulse Z
c1.cd().SetRightMargin(0.05)
c1.cd().SetLeftMargin(0.15)
h_odz = TH1D("h_odz",";pe in veto pulses; pe/bin/sec",100,-4800,4800)
h_odz.GetXaxis().SetLabelSize(0.05)
h_odz.GetXaxis().SetTitleSize(0.06)
h_odz.GetXaxis().SetTitleOffset(1.2)
h_odz.GetYaxis().SetLabelSize(0.05)
h_odz.GetYaxis().SetTitleSize(0.06)
h_odz.GetYaxis().SetTitleOffset(1.2)
#h_odpe.GetZaxis().SetTitleSize(0.05)
if folder=="nodark":
    mcdata.Draw("pescale*odpe>>h_odz","npulses_veto>0 && ")
else:
    mcdata.Draw("pulsez_veto[0]>>h_odz","standard_cuts") # implicit cut, found od pulse
h_odz.Scale( scale )
c1.Update()
c1.SaveAs("figs/%s/h_odz.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_odz.eps"%(folder))
#print "OD Z: ",h_odz.Integral(), " Hz"

# -- RATE --
print "Left over cosmic rate (from h_idpe_wcuts): ",h_idpe_wcuts.Integral()," Hz"
nmuons_excl = mcdata.GetEntries("standard_cuts && ncr_muons>0 && ncr_photons==0 && ncr_electrons==0 && ncr_neutrons==0")
nphotons_excl = mcdata.GetEntries("standard_cuts && ncr_muons==0 && ncr_photons>0 && ncr_electrons==0 && ncr_neutrons==0")
nelectrons_excl = mcdata.GetEntries("standard_cuts && ncr_muons==0 && ncr_photons==0 && ncr_electrons>0 && ncr_neutrons==0")
nneutrons_excl = mcdata.GetEntries("standard_cuts && ncr_muons==0 && ncr_photons==0 && ncr_electrons==0 && ncr_neutrons>0")
print "exclusive: muons=",nmuons_excl*scale,
print "photons=",nphotons_excl*scale,
print "electrons=",nelectrons_excl*scale,
print "neutrons=",nneutrons_excl*scale
nmuons_incl = mcdata.GetEntries("standard_cuts && ncr_muons>0 && (ncr_photons>0 || ncr_electrons>0 || ncr_neutrons>0)")
print "inclusive: muons=",nmuons_incl*scale

raw_input()
