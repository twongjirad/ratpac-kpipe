import os,sys
import ROOT
ROOT.gSystem.Load('libRATEvent')
from ROOT.RAT import DSReader
#from root_numpy import root2array
#import numpy as np
#import pandas as pd
from math import sqrt,exp

#pmtinfo = pd.DataFrame( root2array('../data/kpipe/PMTINFO.root', 'pmtinfo' ) )
#pmtinfo = pmtinfo.set_index('opdetid')

#inputfile = "output_kpipe_0.root"
inputfile = "test.root"
#reader = DSReader('kpipeout_test.root')
#reader = DSReader('../cry/crkpipe.root')
#reader = DSReader("output_kpipe_cryevents_2.root")

if len(sys.argv)==2:
    inputfile = sys.argv[1]
    
fin = ROOT.TFile( inputfile )
mcdata = fin.Get( "mcdata" )
nevents = mcdata.GetEntries()
twfm = ROOT.vector('double')()
#mcdata.SetBranchAddress('twfm',twfm)
mcdata.SetBranchAddress('twfm_veto',twfm)

out = ROOT.TFile('output_test.root','recreate')

c = ROOT.TCanvas('c','',800,1200)
c.Divide(1,3)

window = 40
tave = 40
darkrate = 1.0e6
#nsipms = 100*100
nsipms = 10
nexp_dark = darkrate*1.0e-9*window*nsipms
staterr = sqrt(window*darkrate*1.0e-9*nsipms)
thresh = nexp_dark + 5.0*staterr
decay_const = 30
use_ave = True
print "Expected dark rate in window: ",nexp_dark,"+/-",staterr
print "thresh: ",thresh

#hz = ROOT.TH1D('hz','',200,-5000,5000)
ht = ROOT.TH1D("htraw",'',10000,0,10000) # 10 usec with 1 ns bins
ht2 = ht.Clone("htwinsum") # hit window
htave = ht.Clone("htave") # hit window ave
hthresh = ht.Clone("hthresh")
hthresh_ave = ht.Clone("hthresh_ave")
hexpect = ht.Clone("hexpect")
hexpect_ave = ht.Clone("hexpect_ave")
for h in [ hthresh, hthresh_ave ]:
    h.SetLineColor(ROOT.kGreen+4)



# =====================================
# FILL HISTOGRAMS
for iev in xrange(0,nevents):
    mcdata.GetEntry(iev)
    for ibin in xrange(0,twfm.size()):
        ht.SetBinContent( ibin+1, twfm.at(ibin) )

    print "===================================="
    print "EVENT ",iev
    print "  pe in hist: ",ht.Integral()
    print "  posv: R=", mcdata.rv," Z=",mcdata.zv
    print "  visible energy: ",mcdata.mumomv+mcdata.totkeprotonv," MeV"," mu=",mcdata.mumomv," proton=",mcdata.totkeprotonv
    print "  full window dark noise ",darkrate*1.0e-9*nsipms*10.0e3," vs. integral=",ht.Integral()

    # Analyze
    pulses = {}
    npulses = 0
    active_pulses = []
    for ibin in xrange(window,ht.GetNbinsX()):

        # sum up events in window
        hits_window = 0
        for i in xrange( ibin-window, ibin ):
            hits_window += ht.GetBinContent(i+1)

        # average hits per bin in window
        ave_window = 0.0
        nbinsum = 0
        for i in xrange( max(ibin-tave,0), ibin+tave ):
            ave_window += ht.GetBinContent(i+1)
            nbinsum += 1
        ave_window /= float( nbinsum )
        ht2.SetBinContent( ibin, hits_window )
        htave.SetBinContent( ibin, ave_window )

        pe_expect = nexp_dark
        pe_expect_ave = nexp_dark/float(window)

        if len(active_pulses)==0:
            # looking for new pulse
            #if ave_window*window>thresh:
            if ( use_ave==False and hits_window>(thresh)) or ( use_ave and ave_window>(thresh/float(window))):
                # start new pulse
                active_pulses.append( npulses )
                ll = float(hits_window)
                if use_ave:
                    ll = ave_window                
                pulses[ npulses ] = { "tstart":ibin, "nhits":hits_window, "end":False, "last_level":ll, "nfalling":0 }
                print "PULSE FOUND (",npulses,"): ",ibin,hits_window
                npulses += 1
            hthresh.SetBinContent( ibin, thresh )
            hthresh_ave.SetBinContent( ibin, thresh/float(window) )
        else:

            # look for secondary peak
            mod_thresh = thresh
            mod_thresh_ave = thresh/float(window)
            onerising = False
            for pulse in active_pulses:
                if "peak" not in pulses[pulse]:
                    # we still have a rising peak. make threshold impossble
                    mod_thresh += 2.0*hits_window
                    mod_thresh_ave += 2.0*ave_window
                    onerising = True
                    pe_expect = hits_window
                    pe_expect_ave = ave_window
                else:
                    expecthits = pulses[pulse]["peak"]*exp( -1.0*( ibin-pulses[pulse]["tpeak"] )/decay_const - 0.0*( ibin-pulses[pulse]["tpeak"] )/120.0 )
                    if nexp_dark<5:
                        mod_thresh += (expecthits + 5.0*sqrt(expecthits))
                    else:
                        mod_thresh += (expecthits + 3.0*sqrt(expecthits))
                    expecthits_ave = (pulses[pulse]["peak_ave"])*exp( -1.0*( ibin-pulses[pulse]["tpeak"] )/decay_const)
                    #expecthits_ave = (pulses[pulse]["peak"])
                    #mod_thresh_ave += expecthits_ave + 3.0*sqrt(expecthits_ave)
                    mod_thresh_ave += expecthits_ave
                    pe_expect += expecthits
                    pe_expect_ave += expecthits_ave
            hthresh.SetBinContent( ibin, mod_thresh )
            hthresh_ave.SetBinContent( ibin, mod_thresh_ave )
            # second pulse thresh
            #if hits_window/float(window)>mod_thresh:
            if onerising==False and ( (use_ave==False and hits_window>mod_thresh) or ( use_ave and ave_window>mod_thresh_ave) ):
                # start new pulse
                active_pulses.append( npulses )
                ll = float(hits_window)
                if use_ave:
                    ll = ave_window
                pulses[ npulses ] = { "tstart":ibin, "nhits":hits_window, "end":False, "last_level":ll, "nfalling":0 }
                print "PULSE FOUND (",npulses,") (overlap): ",ibin,hits_window,mod_thresh," rising=",onerising
                npulses += 1
            #if onefall:
            #    print "   mod thresh: ",mod_thresh, "hit_win/win=",float(hits_window)/float(window)," ave_win=",ave_window


            # in pulse. looking for peak and end
            for pulse in active_pulses:
                if "peak" not in pulses[pulse]:
                    # now looking for peak
                    if (not use_ave and hits_window<pulses[pulse]["last_level"]) or ( use_ave and ave_window<pulses[pulse]["last_level"]):
                        pulses[pulse]["nfalling"] += 1
                    else:
                        pulses[pulse]["nfalling"] = 0
                        if not use_ave:
                            pulses[pulse]["last_level"] = float(hits_window)
                        else:
                            pulses[pulse]["last_level"] = ave_window

                    if pulses[pulse]["nfalling"]>3:
                        print "Found peak: ",ibin-window
                        pulses[pulse]["tpeak"] = ibin-3
                        pulses[pulse]["peak"] = float(hits_window)
                        pulses[pulse]["peak_ave"] = ave_window
                else:
                    if ibin >pulses[pulse]["tpeak"]+8*decay_const:
                        print "End of pulse %d found: "%(pulse),ibin
                        pulses[pulse]["tend"] = pulses[pulse]["tpeak"]+8*decay_const # for now

        for ipulse,pulseinfo in pulses.items():
            if "tend" in pulseinfo and ipulse in active_pulses:
                print "Remove pulse=",ipulse
                active_pulses.remove(ipulse)

        # set expectation
        hexpect.SetBinContent( ibin, pe_expect )
        hexpect_ave.SetBinContent( ibin, pe_expect_ave )
                
    c.cd(1)
    ht.Draw()
    ipulses = pulses.keys()
    ipulses.sort()
    lines = []
    fits = []
    for pulse in ipulses:
        pulseinfo = pulses[pulse]
        s = ROOT.TLine( pulseinfo['tstart'], 0, pulseinfo['tstart'], 50 )
        pend = ROOT.TLine( pulseinfo['tend'], 0, pulseinfo['tend'], 50 )
        if "tpeak" in pulseinfo:
            p = ROOT.TLine( pulseinfo['tpeak'], 0, pulseinfo['tpeak'], 50 )
            p.SetLineColor(ROOT.kMagenta)
            lines.append( s )
            lines.append( p )
            lines.append( pend )
            s.Draw()
            pend.Draw()
            p.Draw()
            f = ROOT.TF1("fit%d"%(pulse),"[0]*TMath::Exp( -((x+%f)-[1])/[2] )+%2f"%(window, nexp_dark/window), pulseinfo['tpeak']-window, pulseinfo['tpeak']+decay_const*8-window )
            f.SetParameter(0, pulseinfo['peak'] )
            f.SetParameter(1, pulseinfo['tpeak']-window )
            f.SetParameter(2, decay_const )
            f.Draw("Rsame")
            fits.append(f)

    c.cd(2)
    ht2.Draw()
    hthresh.Draw("same")
    for l in lines:
        l.Draw()
    hexpect.SetLineColor(ROOT.kRed)
    hexpect.Draw("same")

    c.cd(3)
    htave.Draw()
    hexpect_ave.SetLineColor(ROOT.kRed)
    hexpect_ave.Draw("same")
    hthresh_ave.Draw("same")

    c.Update()
    #if len(pulses)>0:
    raw_input()

raw_input()
