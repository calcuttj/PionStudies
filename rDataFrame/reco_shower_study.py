import ROOT as RT
from argparse import ArgumentParser as ap
from math import acos
from array import array
RT.gROOT.SetBatch(1)
RT.gStyle.SetOptStat(0)


parser = ap()
parser.add_argument( "-m", type=str, help='Input MC file' )
parser.add_argument( "-d", type=str, help='Input data file' )
parser.add_argument( "-o", type=str, help='Output file', default='reco_shower_study.root' )

args = parser.parse_args()

fMC = RT.TFile(args.m)
fData = RT.TFile(args.d)
tMC = fMC.Get("pduneana/beamana")
tData = fData.Get("pduneana/beamana")

tMC.Draw("reco_daughter_allShower_startZ>>hMC(50, 0, 300)", "reco_daughter_allShower_ID != -999 && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4")
tData.Draw("reco_daughter_allShower_startZ>>hData(50, 0, 300)", "reco_daughter_allShower_ID != -999 && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4")

hMC = RT.gDirectory.Get("hMC")
hData = RT.gDirectory.Get("hData")

tMC.Draw("reco_daughter_allShower_startZ>>hMCCut(50, 0, 300)", "reco_daughter_allShower_ID != -999 && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4 && !is_pi0_shower")
tData.Draw("reco_daughter_allShower_startZ>>hDataCut(50, 0, 300)", "reco_daughter_allShower_ID != -999 && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4 && !is_pi0_shower")

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
rMC.SetTitle("Secondary Showers;Reconstructed Start Z (cm);Fraction (Not #pi^{0}-like)")
rMC.SetTitleSize(.05, "XY");
rMC.SetTitleSize(.05)
rMC.Draw("hist")
rData.Draw("e1 same")
leg = RT.TLegend()
leg.AddEntry(rMC, "MC", "l")
leg.AddEntry(rData, "Data", "lp")
leg.Draw("same")
c1.Write("cut_ratio")

##
tMC.Draw("acos(reco_daughter_allShower_dirZ)*180/TMath::Pi()>>hMCTheta(45, 0, 180)",
         "reco_daughter_allShower_ID != -999 && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4")
tData.Draw("acos(reco_daughter_allShower_dirZ)*180/TMath::Pi()>>hDataTheta(45, 0, 180)",
                "reco_daughter_allShower_ID != -999 && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4")

hMCTheta = RT.gDirectory.Get("hMCTheta")
hDataTheta = RT.gDirectory.Get("hDataTheta")

tMC.Draw("acos(reco_daughter_allShower_dirZ)*180/TMath::Pi()>>hMCThetaCut(45, 0, 180)",
         "reco_daughter_allShower_ID != -999 && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4 && !is_pi0_shower")
tData.Draw("acos(reco_daughter_allShower_dirZ)*180/TMath::Pi()>>hDataThetaCut(45, 0, 180)",
                "reco_daughter_allShower_ID != -999 && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4 && !is_pi0_shower")


hMCThetaCut = RT.gDirectory.Get("hMCThetaCut")
hDataThetaCut = RT.gDirectory.Get("hDataThetaCut")

rMCTheta = hMCThetaCut.Clone()
rMCTheta.Divide(hMCTheta)
rDataTheta = hDataThetaCut.Clone()
rDataTheta.Sumw2()
hDataTheta.Sumw2()
rDataTheta.Divide(hDataTheta)

rMCTheta.SetLineColor(RT.kRed)
rDataTheta.SetMarkerStyle(20)
rDataTheta.Sumw2()
rDataTheta.SetLineColor(RT.kBlack)
rDataTheta.SetMarkerColor(RT.kBlack)

rMCTheta.SetMinimum(0.)
rMCTheta.SetMaximum(1.1)
rMCTheta.SetTitle("Secondary Showers;Reconstructed #theta;Fraction (Not #pi^{0}-like)")
rMCTheta.SetTitleSize(.05, "XY");
rMCTheta.SetTitleSize(.05)
rMCTheta.Draw("hist")
rDataTheta.Draw("e1 same")
leg.Draw("same")
c1.Write("cut_ratio_theta")



####
#tMC.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi()>>hMCCut2(45, 0, 180)", "reco_daughter_allShower_ID != -999 && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4 && is_secondary_pion")
#tData.Draw("reco_daughter_allTrack_Theta*180./TMath::Pi()>>hDataCut2(45, 0, 180)", "reco_daughter_allShower_ID != -999 && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4 && is_secondary_pion")
#
#hMCCut2 = RT.gDirectory.Get("hMCCut2")
#hDataCut2 = RT.gDirectory.Get("hDataCut2")
#
#rMC2 = hMCCut2.Clone()
#rMC2.Divide(hMC)
#rData2 = hDataCut2.Clone()
#rData2.Sumw2()
#rData2.Divide(hData)
#
#fout.cd()
#rMC2.SetLineColor(RT.kRed)
#rData2.SetMarkerStyle(20)
#rData2.Sumw2()
#rData2.SetLineColor(RT.kBlack)
#rData2.SetMarkerColor(RT.kBlack)
#
#c1.cd()
#rMC2.SetMinimum(0.)
#rMC2.SetMaximum(1.1)
#rMC2.SetTitle("Secondary Tracks;Reconstructed #theta;Fraction (Pionlike)")
#rMC2.SetTitleSize(.05, "XY");
#rMC.SetTitleSize(.05)
#rMC2.Draw("hist")
#rData2.Draw("e1 same")
#leg.Draw("same")
#c1.Write("cut_ratio2")
#
#############
#
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
