import ROOT as RT
from argparse import ArgumentParser as ap
from math import acos
from array import array
RT.gROOT.SetBatch(1)
RT.gStyle.SetOptStat(0)


parser = ap()
parser.add_argument( "-m", type=str, help='Input MC file' )
parser.add_argument( "-d", type=str, help='Input data file' )
parser.add_argument( "-o", type=str, help='Output file', default='reco_pos_study.root' )

args = parser.parse_args()

fMC = RT.TFile(args.m)
fData = RT.TFile(args.d)
tMC = fMC.Get("pduneana/beamana")
tData = fData.Get("pduneana/beamana")

tMC.Draw("reco_daughter_allTrack_startZ>>hMC(25, 0, 250)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")
tData.Draw("reco_daughter_allTrack_startZ>>hData(25, 0, 250)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")

hMC = RT.gDirectory.Get("hMC")
hData = RT.gDirectory.Get("hData")

tMC.Draw("reco_daughter_allTrack_startZ>>hMCCut(25, 0, 250)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4 && !is_secondary_pion")
tData.Draw("reco_daughter_allTrack_startZ>>hDataCut(25, 0, 250)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4 && !is_secondary_pion")

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
rMC.SetTitle("Secondary Tracks;Reconstructed Start Z (cm);Fraction (Not Pionlike)")
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
tMC.Draw("reco_daughter_allTrack_startZ>>hMCPion(25, 0, 250)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4 && is_secondary_pion")
tData.Draw("reco_daughter_allTrack_startZ>>hDataPion(25, 0, 250)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4 && is_secondary_pion")

hMCPion = RT.gDirectory.Get("hMCPion")
hDataPion = RT.gDirectory.Get("hDataPion")

rMCPion = hMCPion.Clone()
rMCPion.Divide(hMC)
rDataPion = hDataPion.Clone()
rDataPion.Sumw2()
rDataPion.Divide(hData)

fout.cd()
rMCPion.SetLineColor(RT.kRed)
rDataPion.SetMarkerStyle(20)
rDataPion.Sumw2()
rDataPion.SetLineColor(RT.kBlack)
rDataPion.SetMarkerColor(RT.kBlack)

c1.cd()
rMCPion.SetMinimum(0.)
rMCPion.SetMaximum(1.1)
rMCPion.SetTitle("Secondary Tracks;Reconstructed Start Z (cm);Fraction (Pionlike)")
rMCPion.SetTitleSize(.05, "XY");
rMCPion.SetTitleSize(.05)
rMCPion.Draw("hist")
rDataPion.Draw("e1 same")
leg.Draw("same")
c1.Write("cut_pion")

############
tMC.Draw("reco_daughter_allTrack_len>>hMCLen(25, 0, 250)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")
tData.Draw("reco_daughter_allTrack_len>>hDataLen(25, 0, 250)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")

hMCLen = RT.gDirectory.Get("hMCLen")
hDataLen = RT.gDirectory.Get("hDataLen")

tMC.Draw("reco_daughter_allTrack_len>>hMCLenCut(25, 0, 250)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4 && !is_secondary_pion")
tData.Draw("reco_daughter_allTrack_len>>hDataLenCut(25, 0, 250)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4 && !is_secondary_pion")

hMCLenCut = RT.gDirectory.Get("hMCLenCut")
hDataLenCut = RT.gDirectory.Get("hDataLenCut")

rMCLen = hMCLenCut.Clone()
rMCLen.Divide(hMCLen)
rDataLen = hDataLenCut.Clone()
rDataLen.Sumw2()
hDataLen.Sumw2()
rDataLen.Divide(hDataLen)

rMCLen.SetLineColor(RT.kRed)
rDataLen.SetMarkerStyle(20)
rDataLen.Sumw2()
rDataLen.SetLineColor(RT.kBlack)
rDataLen.SetMarkerColor(RT.kBlack)

print("cd")
fout.cd()
c2 = RT.TCanvas()
c2.SetTicks()
print("done")
rMCLen.SetMinimum(0.)
rMCLen.SetMaximum(1.1)
rMCLen.SetTitle("Secondary Tracks;Reconstructed Length (cm);Fraction (Not Pionlike)")
rMCLen.SetTitleSize(.05, "XY");
rMCLen.SetTitleSize(.05)
print("drawing")
rMCLen.Draw("hist")
rDataLen.Draw("e1 same")
leg.Draw("same")

print("writing")
c2.Write("len_ratio")
print("wrote")

hMCLen.SetLineColor(RT.kRed)
hDataLen.SetLineColor(RT.kBlack)
hDataLen.SetMarkerColor(RT.kBlack)
hDataLen.SetMarkerStyle(20)
hMCLen.Scale(tData.GetEntries("selection_ID < 4")/tMC.GetEntries("selection_ID < 4"))
hMCLen.SetTitle("Secondary Tracks;Reconstructed Length (cm);Count")
hMCLen.SetTitleSize(.05, "XY")
hMCLen.SetTitleSize(.05)
hDataLen.SetTitle("Secondary Tracks;Reconstructed Length (cm);Count")
hDataLen.SetTitleSize(.05, "XY")
hDataLen.SetTitleSize(.05)

print(hMCLen.GetMaximum(), hDataLen.GetMaximum())
if hMCLen.GetMaximum() > hDataLen.GetMaximum():
  hMCLen.Draw("hist")
  hDataLen.Draw("e1 same")
else:
  hDataLen.Draw("e1")
  hMCLen.Draw("hist same")
leg.Draw("same")
c2.Write("len")

#c1.cd()
#hMC.SetLineColor(RT.kRed)
#hData.SetLineColor(RT.kBlack)
#hData.SetMarkerColor(RT.kBlack)
#hData.SetMarkerStyle(20)
#hMC.Scale(tData.GetEntries("selection_ID < 4")/tMC.GetEntries("selection_ID < 4"))
#hMC.SetTitle("Secondary Tracks;Reconstructed #theta;Count")
#hMC.SetTitleSize(.05, "XY")
#hMC.SetTitleSize(.05)
#hData.SetTitle("Secondary Tracks;Reconstructed #theta;Count")
#hData.SetTitleSize(.05, "XY")
#hData.SetTitleSize(.05)
#if hMC.GetMaximum() > hData.GetMaximum():
#  hMC.Draw("hist")
#  hData.Draw("e1 same")
#else:
#  hData.Draw("e1")
#  hMC.Draw("hist same")
#leg.Draw("same")
#c1.Write("rate")

fout.Close()
