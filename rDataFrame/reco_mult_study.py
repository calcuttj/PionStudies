import ROOT as RT
from argparse import ArgumentParser as ap
from math import acos
from array import array
RT.gROOT.SetBatch(1)
RT.gStyle.SetOptStat(0)


parser = ap()
parser.add_argument( "-m", type=str, help='Input MC file' )
parser.add_argument( "-d", type=str, help='Input data file' )
parser.add_argument( "-o", type=str, help='Output file', default='reco_mult_study.root' )

args = parser.parse_args()

fMC = RT.TFile(args.m)
fData = RT.TFile(args.d)
tMC = fMC.Get("pduneana/beamana")
tData = fData.Get("pduneana/beamana")

tMC.Draw("@reco_daughter_PFP_ID.size()>>hMC(10, 0, 10)", "selection_ID < 4")
tData.Draw("@reco_daughter_PFP_ID.size()>>hData(10, 0, 10)", "selection_ID < 4")

hMC = RT.gDirectory.Get("hMC")
hData = RT.gDirectory.Get("hData")

tMC.Draw("@reco_daughter_PFP_ID.size()>>hMCCut(10, 0, 10)", "selection_ID == 3")
tData.Draw("@reco_daughter_PFP_ID.size()>>hDataCut(10, 0, 10)", "selection_ID == 3")

hMCCut = RT.gDirectory.Get("hMCCut")
hDataCut = RT.gDirectory.Get("hDataCut")

rMC = hMCCut.Clone()
rMC.Divide(hMC)
rData = hDataCut.Clone()
rData.Sumw2()
hData.Sumw2()
rData.Divide(hData)

fout = RT.TFile(args.o, "RECREATE")
fout.cd()
rMC.SetLineColor(RT.kRed)
rData.SetMarkerStyle(20)
rData.Sumw2()
rData.SetLineColor(RT.kBlack)
rData.SetMarkerColor(RT.kBlack)

c1 = RT.TCanvas()
c1.SetTicks()
rMC.SetMinimum(0.)
rMC.SetMaximum(1.1)
rMC.SetTitle(";Multiplicity (PFPs);Fraction Rejected")
rMC.SetTitleSize(.05, "XY");
rMC.SetTitleSize(.05)
rMC.Draw("hist")
rData.Draw("e1 same")
leg = RT.TLegend()
leg.AddEntry(rMC, "MC", "l")
leg.AddEntry(rData, "Data", "lp")
leg.Draw("same")
c1.Write("cut_ratio")
####


tMC.Draw("n_track_daughters>>hMCTrack(10, 0, 10)", "selection_ID < 4")
tData.Draw("n_track_daughters>>hDataTrack(10, 0, 10)", "selection_ID < 4")

hMCTrack = RT.gDirectory.Get("hMCTrack")
hDataTrack = RT.gDirectory.Get("hDataTrack")

tMC.Draw("n_track_daughters>>hMCTrackCut(10, 0, 10)", "selection_ID == 3")
tData.Draw("n_track_daughters>>hDataTrackCut(10, 0, 10)", "selection_ID == 3")

hMCTrackCut = RT.gDirectory.Get("hMCTrackCut")
hDataTrackCut = RT.gDirectory.Get("hDataTrackCut")

rMCTrack = hMCTrackCut.Clone()
rMCTrack.Divide(hMCTrack)
rDataTrack = hDataTrackCut.Clone()
rDataTrack.Sumw2()
hDataTrack.Sumw2()
rDataTrack.Divide(hDataTrack)

fout.cd()
rMCTrack.SetLineColor(RT.kRed)
rDataTrack.SetMarkerStyle(20)
rDataTrack.Sumw2()
rDataTrack.SetLineColor(RT.kBlack)
rDataTrack.SetMarkerColor(RT.kBlack)

rMCTrack.SetMinimum(0.)
rMCTrack.SetMaximum(1.1)
rMCTrack.SetTitle(";Multiplicity (Tracks);Fraction Rejected")
rMCTrack.SetTitleSize(.05, "XY");
rMCTrack.SetTitleSize(.05)
rMCTrack.Draw("hist")
rDataTrack.Draw("e1 same")
leg.Draw("same")
c1.Write("cut_ratio_track")


###
tMC.Draw("n_shower_daughters>>hMCShower(10, 0, 10)", "selection_ID < 4")
tData.Draw("n_shower_daughters>>hDataShower(10, 0, 10)", "selection_ID < 4")

hMCShower = RT.gDirectory.Get("hMCShower")
hDataShower = RT.gDirectory.Get("hDataShower")

tMC.Draw("n_shower_daughters>>hMCShowerCut(10, 0, 10)", "selection_ID == 3")
tData.Draw("n_shower_daughters>>hDataShowerCut(10, 0, 10)", "selection_ID == 3")

hMCShowerCut = RT.gDirectory.Get("hMCShowerCut")
hDataShowerCut = RT.gDirectory.Get("hDataShowerCut")

rMCShower = hMCShowerCut.Clone()
rMCShower.Divide(hMCShower)
rDataShower = hDataShowerCut.Clone()
rDataShower.Sumw2()
hDataShower.Sumw2()
rDataShower.Divide(hDataShower)

fout.cd()
rMCShower.SetLineColor(RT.kRed)
rDataShower.SetMarkerStyle(20)
rDataShower.Sumw2()
rDataShower.SetLineColor(RT.kBlack)
rDataShower.SetMarkerColor(RT.kBlack)

rMCShower.SetMinimum(0.)
rMCShower.SetMaximum(1.1)
rMCShower.SetTitle(";Multiplicity (Showers);Fraction Rejected")
rMCShower.SetTitleSize(.05, "XY");
rMCShower.SetTitleSize(.05)
rMCShower.Draw("hist")
rDataShower.Draw("e1 same")
leg.Draw("same")
c1.Write("cut_ratio_shower")
####

hMCShower.SetLineColor(RT.kRed)
hDataShower.SetMarkerStyle(20)
hDataShower.SetLineColor(RT.kBlack)
hDataShower.SetMarkerColor(RT.kBlack)
hMCShower.Scale(tData.GetEntries("selection_ID < 4")/tMC.GetEntries("selection_ID < 4"))
hMCShower.Draw("hist")
hDataShower.Draw("e1 same")
leg.Draw("same")
c1.Write("cut_shower_rate")

fout.Close()
