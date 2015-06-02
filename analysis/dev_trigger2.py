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
mcdata.SetBranchAddress('twfm',twfm)
#mcdata.SetBranchAddress('twfm_veto',twfm)

out = ROOT.TFile('output_test.root','recreate')

c = ROOT.TCanvas('c','',1200,1200)
c.Divide(2,3)

window = 20
tave = window
darkrate = 1.6e6
nfalling = window*0.5
nabove = 3
nsipms = 100.0*100
nexp_dark = darkrate*1.0e-9*window*nsipms
staterr = sqrt(nexp_dark)
#nexp_dark = 0.0
#staterr = 1.0
thresh = 3.5*staterr
decay_const = 1.0*45+0.0*67
decay_constants = [ (0.5,30.0), (0.2, 90.0), (0.3, 400.0) ]
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

hveto = {}
hveto["htraw"] = ROOT.TH1D("htraw_veto",'',10000,0,10000) # 10 usec with 1 ns bins
hveto["ht2"] = ht.Clone("htwinsum_veto") # hit window
hveto["htave"] = ht.Clone("htave_veto") # hit window ave
hveto["hthresh"] = ht.Clone("hthresh_veto")
hveto["hthresh_ave"] = ht.Clone("hthresh_ave_veto")
hveto["hexpect"] = ht.Clone("hexpect_veto")
hveto["hexpect_ave"] = ht.Clone("hexpect_ave_veto")

for h in [ hthresh, hthresh_ave, hveto["hthresh"], hveto["hthresh_ave"] ]:
    h.SetLineColor(ROOT.kGreen+4)

# =====================================
# FILL HISTOGRAMS
for iev in xrange(0,nevents):
    mcdata.GetEntry(iev)
    for ibin in xrange(0,twfm.size()):
        ht.SetBinContent( ibin+1, twfm.at(ibin) )
    c.cd(1)
    ht.Draw()
    c.Update()

    print "===================================="
    print "EVENT ",iev
    print "  pe in hist: ",ht.Integral()
    print "  posv: R=", mcdata.rv," Z=",mcdata.zv
    print "  visible energy: ",mcdata.mumomv+mcdata.totkeprotonv," MeV"," mu=",mcdata.mumomv," proton=",mcdata.totkeprotonv
    print "  predark: ID=",mcdata.predark_idpe, " OD=",mcdata.predark_odpe
    print "  full window dark noise ",darkrate*1.0e-9*nsipms*10.0e3," vs. integral=",ht.Integral()
    print "  c++ result: id pulses=",mcdata.npulses," od pulses=",mcdata.npulses_veto
    print "  thresh (ave): ",thresh/float(window)
    for i in xrange(0,mcdata.npulses):
        print " ",i,") t=",mcdata.ttrig[i]," pe=",mcdata.pulsepe[i]

    if iev not in [10,19,37,69,103,112,121,138]:
        continue
    
    # Analyze
    pulses = {}
    npulses = 0
    active_pulses = []
    bins_above = 0
    for ibin in xrange(window,ht.GetNbinsX()):

        # sum up events in window
        hits_window = 0
        for i in xrange( ibin-window, ibin ):
            hits_window += ht.GetBinContent(i+1)

        # average hits per bin in window
        ave_window = 0.0
        nbinsum = 0
        for i in xrange( max(int(ibin-0.5*tave),0), min( int(ibin+0.5*tave),ht.GetNbinsX() ) ):
            ave_window += ht.GetBinContent(i+1)
            nbinsum += 1
        ave_window /= float( nbinsum )
        ht2.SetBinContent( ibin, hits_window )
        htave.SetBinContent( ibin, ave_window )

        #pe_expect = nexp_dark
        #pe_expect_ave = nexp_dark/float(window)
        pe_expect = 0.0
        pe_expect_ave = 0.0

        if len(active_pulses)==0:
            # looking for new pulse
            #if ave_window*window>thresh:
            if ( use_ave==False and hits_window>(thresh)) or ( use_ave and ave_window>(thresh/float(window))):
                if bins_above<nabove:
                    bins_above += 1
                else:
                    # start new pulse
                    active_pulses.append( npulses )
                    ll = float(hits_window)
                    if use_ave:
                        ll = ave_window                
                    pulses[ npulses ] = { "tstart":ibin, "nhits":hits_window, "end":False, "last_level":ll, "nfalling":0,"tlast":ibin }
                    print "PULSE FOUND (",npulses,"): ",ibin,hits_window
                    npulses += 1
                    bins_above = 0
            else:
                # drops below
                bins_above = 0
            hthresh.SetBinContent( ibin, thresh )
            hthresh_ave.SetBinContent( ibin, thresh/float(window) )
        else:

            # look for secondary peak
            mod_thresh = 0.0 #thresh
            mod_thresh_ave = 0.0 #thresh/float(window)
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
                    arg = 0.0
                    for w,dc in decay_constants:
                        arg += w*( ( float(ibin)-pulses[pulse]["tpeak"] )/dc )
                    expecthits = pulses[pulse]["peak"]*exp( -1.0*arg )
                    expect_sig = sqrt(nexp_dark)
                    if expecthits<0:
                        expecthits = 1.0
                    if nexp_dark<5:
                        mod_thresh += (expecthits + 5.0*expect_sig)
                    else:
                        mod_thresh += (expecthits + 3.0*expect_sig)
                    expecthits_ave = pulses[pulse]["peak_ave"]*exp( -1.0*arg )
                    mod_thresh_ave += expecthits_ave + 4.0*(staterr/float(window)) + sqrt(expecthits_ave)
                    pe_expect += expecthits
                    pe_expect_ave += expecthits_ave
            mod_thresh_ave = max( mod_thresh_ave,  thresh/float(window) )
            hthresh.SetBinContent( ibin, mod_thresh )
            hthresh_ave.SetBinContent( ibin, mod_thresh_ave )
            # second pulse thresh
            #if hits_window/float(window)>mod_thresh:
            if onerising==False and ( (use_ave==False and hits_window>mod_thresh) or ( use_ave and ave_window>mod_thresh_ave) ):
                if bins_above<nabove:
                    bins_above += 1
                else:
                    # start new pulse
                    active_pulses.append( npulses )
                    ll = float(hits_window)
                    if use_ave:
                        ll = ave_window
                    pulses[ npulses ] = { "tstart":ibin, "nhits":hits_window, "end":False, "last_level":ll, "nfalling":0,"last_max":0,"tlast":ibin }
                    print "PULSE FOUND (",npulses,") (overlap): ",ibin,hits_window,mod_thresh," rising=",onerising
                    npulses += 1
                    bins_above = 0
            else:
                bins_above = 0
                
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
                            if pulses[pulse]["last_level"]<ave_window:
                                pulses[pulse]["last_level"] = ave_window
                                pulses[pulse]["tlast"] = ibin

                    if pulses[pulse]["nfalling"]>nfalling:
                        print "Found peak: ",ibin-window
                        pulses[pulse]["tpeak"] = pulses[pulse]["tlast"]
                        pulses[pulse]["peak"] = pulses[pulse]["last_level"]*window #float(hits_window) 
                        pulses[pulse]["peak_ave"] = pulses[pulse]["last_level"] # ave win
                else:
                    # LOOK FOR END
                    dt = ibin-pulses[pulse]["tstart"]
                    # below threshold
                    if ave_window<thresh/float(window):
                        if ibin-pulses[pulse]["tstart"]<50 or pulses[pulse]["nfalling"]<nfalling:
                            print "pulse too small: dt=",dt
                            pulses[pulse]["reject"] = True
                        
                    # fall below expectation
                    elif dt>50 and pe_expect_ave < (thresh/float(window))*1.05:
                        print "End of pulse %d found via expectation (dt="%(pulse),ibin-pulses[pulse]["tstart"], "): ",ibin
                        pulses[pulse]["tend"] = ibin+decay_const
                        
                    elif ibin >pulses[pulse]["tpeak"]+12*decay_const:
                        print "End of pulse %d found via itme out: "%(pulse),ibin
                        pulses[pulse]["tend"] = pulses[pulse]["tpeak"]+10*decay_const # for now
                    #if len(active_pulses)==1 and ave_window>pulses[pulse]["peak_ave"]:
                    #    print "adjust peak: ",ave_window, pulses[pulse]["peak_ave"]
                    #    pulses[pulse]["peak_ave"] = ave_window
                    #    pulses[pulse]["peak"] = ave_window*window
                    #    pulses[pulse]["tpeak"] = ibin

        for ipulse,pulseinfo in pulses.items():
            if "tend" in pulseinfo and ipulse in active_pulses:
                print "Remove pulse (OK) =",ipulse
                active_pulses.remove(ipulse)
            if "reject" in pulseinfo and ipulse in active_pulses:
                print "Remove pulse (REJECTED) =",ipulse
                active_pulses.remove(ipulse)
                pulses.pop(ipulse,None)

        # set expectation
        hexpect.SetBinContent( ibin, pe_expect )
        hexpect_ave.SetBinContent( ibin, pe_expect_ave )

    print "Final number of pulses: ",len(pulses)
    pulseids = pulses.keys()
    pulseids.sort()
    for id, pulse in pulses.items():
        print id,pulse
                
    c.cd(1)
    ht.Draw()
    ipulses = pulses.keys()
    ipulses.sort()
    lines = []
    fits = []
    for pulse in ipulses:
        pulseinfo = pulses[pulse]
        s = ROOT.TLine( pulseinfo['tstart'], 0, pulseinfo['tstart'], 50 )
        try:
            pend = ROOT.TLine( pulseinfo['tend'], 0, pulseinfo['tend'], 50 )
        except:
            end = ht.GetNbinsX()
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

    c.cd(3)
    ht2.Draw()
    hthresh.Draw("same")
    for l in lines:
        l.Draw()
    hexpect.SetLineColor(ROOT.kRed)
    hexpect.Draw("same")

    c.cd(5)
    htave.Draw()
    hexpect_ave.SetLineColor(ROOT.kRed)
    hexpect_ave.Draw("same")
    hthresh_ave.Draw("same")
    for l in lines:
        l.Draw()
    c.Update()
    #if len(pulses)>0:
    raw_input()

raw_input()
