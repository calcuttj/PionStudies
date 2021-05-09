import ROOT as RT
from argparse import ArgumentParser as ap
from math import acos
from array import array
RT.gROOT.SetBatch(1)
RT.gStyle.SetOptStat(0)


parser = ap()
parser.add_argument( "-m", type=str, help='Input MC file' )
parser.add_argument( "-d", type=str, help='Input data file' )
parser.add_argument( "-o", type=str, help='Output file', default='reco_angular_study.root' )
parser.add_argument( "-n", type=str, help='', default='45')

args = parser.parse_args()

fMC = RT.TFile(args.m)
fData = RT.TFile(args.d)
tMC = fMC.Get("pduneana/beamana")
tData = fData.Get("pduneana/beamana")

tMC.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi()>>hMC("+args.n+", 0, 180)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")
tData.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi()>>hData("+args.n+", 0, 180)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")

hMC = RT.gDirectory.Get("hMC")
hData = RT.gDirectory.Get("hData")

tMC.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi()>>hMCCut("+args.n+", 0, 180)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4 && !is_secondary_pion")
tData.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi()>>hDataCut("+args.n+", 0, 180)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4 && !is_secondary_pion")

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
rMC.SetTitle("Secondary Tracks;Reconstructed #theta;Fraction (Not Pionlike)")
rMC.SetTitleSize(.05, "XY");
rMC.SetTitleSize(.05)
rMC.Draw("hist")
rData.Draw("e1 same")
leg = RT.TLegend()
leg.AddEntry(rMC, "MC", "l")
leg.AddEntry(rData, "Data", "lp")
leg.Draw("same")
c1.Write("cut_ratio")
rMC.Write()


####
tMC.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi()>>hMCCut2("+args.n+", 0, 180)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4 && is_secondary_pion")
tData.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi()>>hDataCut2("+args.n+", 0, 180)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4 && is_secondary_pion")

hMCCut2 = RT.gDirectory.Get("hMCCut2")
hDataCut2 = RT.gDirectory.Get("hDataCut2")

rMC2 = hMCCut2.Clone()
rMC2.Divide(hMC)
rData2 = hDataCut2.Clone()
rData2.Sumw2()
rData2.Divide(hData)

fout.cd()
rMC2.SetLineColor(RT.kRed)
rData2.SetMarkerStyle(20)
rData2.Sumw2()
rData2.SetLineColor(RT.kBlack)
rData2.SetMarkerColor(RT.kBlack)

c1.cd()
rMC2.SetMinimum(0.)
rMC2.SetMaximum(1.1)
rMC2.SetTitle("Secondary Tracks;Reconstructed #theta;Fraction (Pionlike)")
rMC2.SetTitleSize(.05, "XY");
rMC.SetTitleSize(.05)
rMC2.Draw("hist")
rData2.Draw("e1 same")
leg.Draw("same")
c1.Write("cut_ratio2")

############

c1.cd()
hMC.SetLineColor(RT.kRed)
hData.SetLineColor(RT.kBlack)
hData.SetMarkerColor(RT.kBlack)
hData.SetMarkerStyle(20)
hMC.Scale(tData.GetEntries("selection_ID < 4")/tMC.GetEntries("selection_ID < 4"))
hMC.SetTitle("Secondary Tracks;Reconstructed #theta;Count")
hMC.SetTitleSize(.05, "XY")
hMC.SetTitleSize(.05)
hData.SetTitle("Secondary Tracks;Reconstructed #theta;Count")
hData.SetTitleSize(.05, "XY")
hData.SetTitleSize(.05)
if hMC.GetMaximum() > hData.GetMaximum():
  hMC.Draw("hist")
  hData.Draw("e1 same")
else:
  hData.Draw("e1")
  hMC.Draw("hist same")
leg.Draw("same")
c1.Write("rate")

##################################
tMC.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi():reco_beam_endZ>>hMC2D(20, 0, 240, 30, 0, 180)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")
tData.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi():reco_beam_endZ>>hData2D(20, 0, 240, 30, 0, 180)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")

hMC2D = RT.gDirectory.Get("hMC2D")
hData2D = RT.gDirectory.Get("hData2D")

tMC.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi():reco_beam_endZ>>hMC2DCut(20, 0, 240, 30, 0, 180)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID == 3")
tData.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi():reco_beam_endZ>>hData2DCut(20, 0, 240, 30, 0, 180)", "reco_daughter_allTrack_Theta > -999 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID == 3")

hMC2DCut = RT.gDirectory.Get("hMC2DCut")
hData2DCut = RT.gDirectory.Get("hData2DCut")

rMC2DCut = hMC2DCut.Clone()
rMC2DCut.Divide(hMC2D)
rData2DCut = hData2DCut.Clone()
rData2DCut.Divide(hData2D)

rMC2DCut.Write()
rData2DCut.Write()

###
fout.Close()
