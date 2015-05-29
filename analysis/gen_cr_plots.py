import os,sys
from ROOT import *
import array

gStyle.SetOptStat(0)

elasped_time = 1.17161308983 # seconds
search_window = 250.0e-9 # ns

#folder="veto_wdark"
#folder="veto_nodark"
#folder="veto_wdark_v2"
#folder="veto_nodark_v2"
#folder="veto_wdark_v3"
#folder="veto_wdark_v3_10mhz"
#folder="veto_wdark_v3_10mhz_oldalg"
folder="veto_wdark_v4_10mhz"

os.system("mkdir -p figs/%s/eps"%(folder))

if folder=="veto_nodark":
    #tf = TFile("crana_merged_nodarknoise_0_1999.root")
    #elasped_time = 1.17161308983
    tf = TFile("crana_merged_nodarknoise_0_500_v1.root")
    elasped_time = 0.293815160179
elif folder=="veto_nodark_v2":
    tf = TFile("crana_merged_nodarknoise_0_500_v2.root")
    elasped_time = 0.2920053186
elif folder=="veto_wdark":
    #tf = TFile("crana_merged_wdarknoise_0_1999.root") # darknoise
    tf = TFile("crana_merged_wdarknoise_0_500_v1b.root") # darknoise
    elasped_time = 0.292005318686
elif folder=="veto_wdark_v2":
    tf = TFile("crana_merged_wdarknoise_0_500_v2.root") # darknoise
    elasped_time = 0.2926 ##
elif folder=="veto_wdark_v3":
    tf = TFile("crana_merged_wdarknoise_0_2000_v3.root") # darknoise
    elasped_time = 1.12006 ##
elif folder=="veto_wdark_v3_10mhz":
    tf = TFile("crana_merged_wdarknoise_0_2000_v3_10mhz.root") # darknoise
    elasped_tim = 1.15272276033
elif folder=="veto_wdark_v3_10mhz_oldalg":
    tf = TFile("crana_merged_wdarknoise_0_500_v3_10mhz_oldalg.root") # darknoise
    elasped_tim = 0.293815160179
elif folder=="veto_wdark_v4_10mhz":
    tf = TFile("crana_merged_wdarknoise_0_2000_v4_10mhz.root") # darknoise
    elasped_tim = 1.07652645625
else:
    print "wrong file type",folder
    sys.exit(-1)
mcdata = tf.Get("mcdata")

scale = 1.0/elasped_time
#scale = 1.0
pe_scale = 4500.0/8500.0
mcdata.SetAlias("pescale","4500.0/8500.0")
mcdata.SetAlias("standard_cuts","npulses==2 && abs(pulsez[0]-pulsez[1])<200.0 && pulsepe[0]*pescale>200 && pulsepe[0]*pescale<1300 && pulsepe[1]*pescale>100.0 && pulsepe[1]*pescale<800.0 && pulse_totodpe*pescale<=0")
mcdata.SetAlias("nood_cuts",    "npulses==2 && abs(pulsez[0]-pulsez[1])<200.0 && pulsepe[0]*pescale>200 && pulsepe[0]*pescale<1300 && pulsepe[1]*pescale>100.0 && pulsepe[1]*pescale<800.0")
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

mcdata.Draw("pescale*pulsepe[0]>>h_idpe","npulses>0")
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
ncr_cuts = {"ncr_muons":"ncr_muons>0 && ncr_photons==0 && ncr_electrons==0 && ncr_neutrons==0",
            "ncr_photons":"ncr_muons==0 && ncr_photons>0 && ncr_electrons==0 && ncr_neutrons==0",
            "ncr_electrons":"ncr_muons==0 && ncr_photons==0 && ncr_electrons>0 && ncr_neutrons==0",
            "ncr_neutrons":"ncr_muons==0 && ncr_photons==0 && ncr_electrons==0 && ncr_neutrons>0" }

for ncr in ["ncr_muons","ncr_photons","ncr_electrons","ncr_neutrons"]:
    h = TH1D("h_%s"%(ncr),";pe in prompt pulse; fraction of events",10,0,10)
    mcdata.Draw("%s>>h_%s"%(ncr,ncr),"standard_cuts && %s"%(ncr_cuts[ncr]))
    c1.Update()
    c1.SaveAs("figs/%s/h_%s.pdf"%(folder,ncr))


# Michel Pulse PE
c1.cd().SetRightMargin(0.05)
c1.cd().SetLeftMargin(0.15)
h_michel_scale = TH1D("h_michel_scale",";pe in Michel pulse; fraction of events",50,0,800)
h_michel_scale.GetXaxis().SetLabelSize(0.05)
h_michel_scale.GetXaxis().SetTitleSize(0.06)
h_michel_scale.GetXaxis().SetTitleOffset(1.2)
h_michel_scale.GetYaxis().SetLabelSize(0.05)
h_michel_scale.GetYaxis().SetTitleSize(0.06)
h_michel_scale.GetYaxis().SetTitleOffset(1.2)
#h_michel_scale.GetZaxis().SetTitleSize(0.05)
mcdata.Draw("pescale*pulsepe[1]>>h_michel_scale","npulses==2")
h_michel_scale.Scale( scale )
c1.Update()
c1.SaveAs("figs/%s/h_michel_scale.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_michel_scale.eps"%(folder))

# OD PE
c1.cd().SetLogy(1)
c1.cd().SetRightMargin(0.05)
c1.cd().SetLeftMargin(0.15)
h_odpe = TH1D("h_odpe",";pe in veto pulses; pe/bin/sec",100,0,1000)
h_odpe_zoom = TH1D("h_odpe_zoom",";pe in veto pulses; pe/bin/sec",200,0,400)
h_odpe_cuts = TH1D("h_odpe_wcuts",";pe in veto pulses; pe/bin/sec",200,0,400)
h_pre_odpe = TH1D("h_odpe_pre",";pe in veto pulses; pe/bin/sec",200,0,1000)
h_pulse_odpe = TH2D("h_pulse_odpe",";true pe in veto; total pe in veto pulses",200,0,1000,200,0,1000)
h_missed_odpe = TH1D("h_missed_odpe",";true pe in veto of missed veto events",200,0,1000)
for h in [h_odpe,h_odpe_cuts,h_odpe_zoom,h_pulse_odpe,h_missed_odpe]:
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(1.2)
h_odpe.SetMinimum(2.0)
h_odpe.SetMaximum(3.0e5)
h_pre_odpe.SetMinimum(10.0)
h_pre_odpe.SetMaximum(3.0e5)
h_missed_odpe.SetMinimum(0.5)
h_missed_odpe.SetMaximum(3.0e5)
mcdata.Draw("pulse_totodpe>>h_odpe_wcuts","nood_cuts")
mcdata.Draw("pulse_totodpe>>h_odpe","")
mcdata.Draw("pulse_totodpe>>h_odpe_zoom","")
mcdata.Draw("predark_odpe>>h_odpe_pre","")
mcdata.Draw("predark_odpe>>h_missed_odpe","pulse_totodpe<=0")
mcdata.Draw("pulse_totodpe:predark_odpe>>h_pulse_odpe","","COLZ")
h_odpe.Scale( scale )
h_odpe_cuts.Scale( scale )
h_odpe.Draw()
c1.Update()
c1.SaveAs("figs/%s/h_odpe.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_odpe.eps"%(folder))
h_odpe_cuts.Draw()
c1.SaveAs("figs/%s/h_odpe_wcuts.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_odpe_wcuts.eps"%(folder))
h_odpe_zoom.Draw()
c1.SaveAs("figs/%s/h_odpe_zoom.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_odpe_zoom.eps"%(folder))
h_missed_odpe.Draw()
c1.SaveAs("figs/%s/h_missed_odpe.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_missed_odpe.eps"%(folder))
c1.cd().SetLogy(0)
h_pulse_odpe.Draw("COLZ")
c1.SaveAs("figs/%s/h2_odpe_scale.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_odpe_scale.eps"%(folder))


# OD Pulse Z
c1.cd().SetRightMargin(0.05)
c1.cd().SetLeftMargin(0.15)
h_odz = TH1D("h_odz",";pe in veto pulses; pe/bin/sec",56,-5000,5000)
h_odz.GetXaxis().SetLabelSize(0.05)
h_odz.GetXaxis().SetTitleSize(0.06)
h_odz.GetXaxis().SetTitleOffset(1.2)
h_odz.GetYaxis().SetLabelSize(0.05)
h_odz.GetYaxis().SetTitleSize(0.06)
h_odz.GetYaxis().SetTitleOffset(1.2)
#h_odpe.GetZaxis().SetTitleSize(0.05)
mcdata.Draw("pulsez_veto[0]>>h_odz","nood_cuts")
h_odz.Scale( scale )
c1.Update()
c1.SaveAs("figs/%s/h_odz.pdf"%(folder))
c1.SaveAs("figs/%s/eps/h_odz.eps"%(folder))
#print "OD Z: ",h_odz.Integral(), " Hz"

# -- RATE --
print "Left over cosmic rate (from h_idpe_wcuts): ",h_idpe_wcuts.Integral()," Hz","(",h_idpe_wcuts.GetEntries(),")"
nmuons_excl = mcdata.GetEntries("standard_cuts && ncr_muons>0 && ncr_photons==0 && ncr_electrons==0 && ncr_neutrons==0")
nphotons_excl = mcdata.GetEntries("standard_cuts && ncr_muons==0 && ncr_photons>0 && ncr_electrons==0 && ncr_neutrons==0")
nelectrons_excl = mcdata.GetEntries("standard_cuts && ncr_muons==0 && ncr_photons==0 && ncr_electrons>0 && ncr_neutrons==0")
nneutrons_excl = mcdata.GetEntries("standard_cuts && ncr_muons==0 && ncr_photons==0 && ncr_electrons==0 && ncr_neutrons>0")
tot = nmuons_excl+nphotons_excl+nelectrons_excl+nneutrons_excl
print "exclusive: muons=",nmuons_excl*scale,"(",nmuons_excl/float(tot),")"
print "photons=",nphotons_excl*scale,"(",nphotons_excl/float(tot),")"
print "electrons=",nelectrons_excl*scale,"(",nelectrons_excl/float(tot),")"
print "neutrons=",nneutrons_excl*scale,"(",nneutrons_excl/float(tot),")"
nmuons_incl = mcdata.GetEntries("standard_cuts && ncr_muons>0 && (ncr_photons>0 || ncr_electrons>0 || ncr_neutrons>0)")
print "inclusive: muons=",nmuons_incl*scale

raw_input()
