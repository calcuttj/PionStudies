import ROOT as RT
from argparse import ArgumentParser as ap
from array import array

RT.gROOT.SetBatch(1)
RT.gStyle.SetOptStat(0)

def draw_vertical_line(x, color):
  #l = RT.TLine(l, 0., l, hDataEndZ.GetMaximum())
  l = RT.TLine()
  lm = RT.gPad.GetLeftMargin()
  rm = 1. - RT.gPad.GetRightMargin()
  tm = 1. - RT.gPad.GetTopMargin()
  bm = RT.gPad.GetBottomMargin()
  xndc = (rm - lm)*((x - RT.gPad.GetUxmin())/(RT.gPad.GetUxmax() - RT.gPad.GetUxmin())) + lm
  l.SetLineColor(color)
  l.DrawLineNDC(xndc, bm, xndc, tm)
  

parser = ap()
parser.add_argument( "-m", type=str, help='MC file' )
parser.add_argument( "-d", type=str, help='Data file' )
parser.add_argument( "-o", type=str, help='Output file', default='full_study.root')
parser.add_argument( "--zoom_endZ", type=bool, help='', default = False)

args = parser.parse_args()

fMC = RT.TFile(args.m, "OPEN")

t = fMC.Get("pduneana/beamana")

fData = RT.TFile(args.d, "OPEN")
tData = fData.Get("pduneana/beamana")

scale = tData.GetEntries() / t.GetEntries()

hEndZs = []
sEndZ = RT.THStack()
hCuts = []
sCut = RT.THStack()

colors = [
 602, 433, 837, 419,
 403, 795, 630, 629,
 893, 619, 397, 632
]


hMichelScoreBeams = []
sMichelScoreBeams = RT.THStack()
lBeam = RT.TLegend()

#tBeams = ["Abs", "Cex", "Inel", "UpstreamInt", "Muons", "PastFV", "Other"]
tBeams = ["Absorption", "Charge Exchange", "BG Interaction", "Upstream Interaction", "Muons", "Past FV", "Other"]

mc_200_220 = 0.
mc_220_500 = 0.
bins = [i*2. for i in range(0, 112)]
bins += [234. + i*12 for i in range(0, 24) ]
#print(bins)
if args.zoom_endZ:
  print('zooming')
  hDataEndZ = RT.TH1D("hDataEndZ", "", len(bins)-1, array("d", bins))
  tData.Draw("reco_beam_endZ>>hDataEndZ", "selection_ID < 5")
else:
  tData.Draw("reco_beam_endZ>>hDataEndZ(500, 0, 500)", "selection_ID < 5")
  hDataEndZ = RT.gDirectory.Get("hDataEndZ")
hDataEndZ.SetMarkerStyle(20)
hDataEndZ.SetMarkerColor(RT.kBlack)
hDataEndZ.SetLineColor(RT.kBlack)
bin_200 = hDataEndZ.GetXaxis().FindBin(201.)
data_200_220 = hDataEndZ.Integral(bin_200, bin_200 + 1)
data_220_500 = hDataEndZ.Integral(bin_200+2, 500)


hDataFlow = RT.TH1D("DataFlow", "", 3, 0, 3)
hDataFlow.SetBinContent(1, tData.GetEntries())
hDataFlow.SetBinContent(2, tData.GetEntries("selection_ID < 6"))
hDataFlow.SetBinContent(3, tData.GetEntries("selection_ID < 5"))
hDataFlow.SetMarkerStyle(20)
hDataFlow.SetMarkerColor(RT.kBlack)
hDataFlow.SetLineColor(RT.kBlack)
sFlows = RT.THStack()

for i in range(1, 8):
  if args.zoom_endZ:
    h = RT.TH1D("hMCEndZ" + str(i), "", len(bins)-1, array("d", bins))
    t.Draw("reco_beam_endZ>>hMCEndZ" + str(i),
           "new_interaction_topology == " + str(i) + " && selection_ID < 5")
    hEndZs.append(h)
  else:
    t.Draw("reco_beam_endZ>>hMCEndZ" + str(i) + "(500, 0, 500)",
           "new_interaction_topology == " + str(i) + " && selection_ID < 5")
    hEndZs.append(RT.gDirectory.Get("hMCEndZ" + str(i)))
  hEndZs[-1].SetFillColor(colors[i-1])
  hEndZs[-1].SetLineColor(colors[i-1])
  hEndZs[-1].Scale(scale)
  mc_200_220 += hEndZs[-1].Integral(bin_200, bin_200 + 1)
  mc_220_500 += hEndZs[-1].Integral(bin_200 + 2, 500)
  sEndZ.Add(hEndZs[-1])

  t.Draw("reco_daughter_PFP_michelScore_collection>>hMichelScoreBeam" + str(i) + "(25, 0, 1)",
         "new_interaction_topology == " + str(i) + " && selection_ID < 3 && reco_daughter_PFP_trackScore_collection < .3")
  hMichelScoreBeams.append(RT.gDirectory.Get("hMichelScoreBeam" + str(i)))
  hMichelScoreBeams[-1].SetFillColor(colors[i-1])
  hMichelScoreBeams[-1].SetLineColor(colors[i-1])
  hMichelScoreBeams[-1].Scale(scale)
  sMichelScoreBeams.Add(hMichelScoreBeams[-1])

  lBeam.AddEntry(hMichelScoreBeams[-1], tBeams[i-1], "f")

  t.Draw("selection_ID>>hMCCut" + str(i) + "(2, 5, 7)", "selection_ID > 4 && new_interaction_topology == " + str(i))
  hCuts.append(RT.gDirectory.Get("hMCCut" + str(i)))
  hCuts[-1].SetFillColor(colors[i-1])
  hCuts[-1].SetLineColor(colors[i-1])
  hCuts[-1].Scale(scale)
  sCut.Add(hCuts[-1])

  hFlow = RT.TH1D("hFlow" + str(i), "", 3, 0, 3)
  hFlow.SetBinContent(1, t.GetEntries("new_interaction_topology == " + str(i)))
  hFlow.SetBinContent(2, t.GetEntries("selection_ID < 6 && new_interaction_topology == " + str(i)))
  hFlow.SetBinContent(3, t.GetEntries("selection_ID < 5 && new_interaction_topology == " + str(i)))
  hFlow.SetFillColor(colors[i-1])
  hFlow.SetLineColor(colors[i-1])
  hFlow.Scale(scale)
  sFlows.Add(hFlow)
   

fout = RT.TFile(args.o, "RECREATE")

cFlow = RT.TCanvas("cFlow", "cFlow")
sFlows.Draw("hist")
sFlows.GetHistogram().SetTitle(";;")
sFlows.GetHistogram().GetXaxis().SetBinLabel(1, "All")
sFlows.GetHistogram().GetXaxis().SetBinLabel(2, "Pandora")
sFlows.GetHistogram().GetXaxis().SetBinLabel(3, "Beam Quality")
hDataFlow.Draw("same pe1")
cFlow.Write()

tData.Draw("selection_ID>>hDataCut(2, 5, 7)", "selection_ID > 4")
hDataCut = RT.gDirectory.Get("hDataCut")
print(hDataCut.GetBinContent(1), hDataCut.GetBinContent(2))
hDataCut.SetMarkerStyle(20)
hDataCut.SetMarkerColor(RT.kBlack)
hDataCut.SetLineColor(RT.kBlack)

cCut = RT.TCanvas("cCut", "cCut")
cCut.SetTicks()
sCut.Draw("hist")
sCut.GetHistogram().GetXaxis().SetBinLabel(1, "Fails Beam Cuts")
sCut.GetHistogram().GetXaxis().SetBinLabel(2, "No Beam Track")
sCut.GetHistogram().GetXaxis().SetLabelSize(.04)
hDataCut.Draw("same pe1")
lBeam.AddEntry(hDataCut, "Data", "lp")
lBeam.Draw("same")
cCut.Write()

cEndZ = RT.TCanvas("cEndZ", "cEndZ")
cEndZ.SetTicks()
sEndZ.Draw("hist")
sEndZ.GetHistogram().SetTitle(";Reconstructed End Z (cm)")
if (hDataEndZ.GetMaximum() > sEndZ.GetHistogram().GetMaximum()):
  hDataEndZ.SetTitle(";Reconstructed End Z (cm)")
  hDataEndZ.Draw("pe1")
  sEndZ.Draw("hist same")
  hDataEndZ.Draw("pe1 same")
  #lLow = RT.TLine(222., 0., 222., hDataEndZ.GetMaximum())
else:
  sEndZ.Draw("hist")
  hDataEndZ.Draw("same pe1")
  #lLow = RT.TLine(222., 0., 222., sEndZ.GetHistogram().GetMaximum())

draw_vertical_line(222., RT.kBlack)
#lLow.SetLineColor(RT.kBlack)
#lLow.Draw("same")
lBeam.Draw("same")

cEndZ.Write()


tData.Draw("reco_daughter_PFP_michelScore_collection>>hDataMichelScore(25, 0, 1)",
       "selection_ID < 3 && reco_daughter_PFP_trackScore_collection < .3")
hDataMichelScore = RT.gDirectory.Get("hDataMichelScore")
hDataMichelScore.SetMarkerStyle(20)
hDataMichelScore.SetMarkerColor(RT.kBlack)
hDataMichelScore.SetLineColor(RT.kBlack)


cMichelScoreBeam = RT.TCanvas("cMichelScoreBeam", "cMichelScoreBeam")
cMichelScoreBeam.SetTicks()
sMichelScoreBeams.Draw("hist")
hDataMichelScore.Draw("same pe1")
lBeam.Draw("same")
cMichelScoreBeam.Write()

hBeamXSamples = []
hBeamYSamples = []
hBeamZSamples = []
hBeamCosSamples = []

sBeamXSamples = RT.THStack()
sBeamYSamples = RT.THStack()
sBeamZSamples = RT.THStack()
sBeamCosSamples = RT.THStack()

tBeamBTs = ["Primary #pi^{+}", "Primary #mu^{+}", "Cosmic", 
            "Beam Background", "Upstream Interaction", "Other"]
lBeamBT = RT.TLegend()


hBeamTrackScores = []
sBeamTrackScores = RT.THStack()
for i in range(1, 7):
  t.Draw("reco_beam_startX - beam_inst_X>>hBeamXSample" + str(i) + "(60, -30, 30)",
  #t.Draw("reco_beam_calo_startX - beam_inst_X>>hBeamXSample" + str(i) + "(60, -30, 30)",
         "beam_backtrack == " + str(i)) 
  hBeamXSamples.append(RT.gDirectory.Get("hBeamXSample" + str(i)))
  hBeamXSamples[-1].SetFillColor(colors[i-1])
  hBeamXSamples[-1].SetLineColor(colors[i-1])
  hBeamXSamples[-1].Scale(scale)
  sBeamXSamples.Add(hBeamXSamples[-1])

  t.Draw("reco_beam_startY - beam_inst_Y>>hBeamYSample" + str(i) + "(120, -30, 30)",
  #t.Draw("reco_beam_calo_startY - beam_inst_Y>>hBeamYSample" + str(i) + "(120, -30, 30)",
         "beam_backtrack == " + str(i)) 
  hBeamYSamples.append(RT.gDirectory.Get("hBeamYSample" + str(i)))
  hBeamYSamples[-1].SetFillColor(colors[i-1])
  hBeamYSamples[-1].SetLineColor(colors[i-1])
  hBeamYSamples[-1].Scale(scale)
  sBeamYSamples.Add(hBeamYSamples[-1])

  t.Draw("reco_beam_startZ - beam_inst_Z>>hBeamZSample" + str(i) + "(100, 0, 50)",
  #t.Draw("reco_beam_calo_startZ - beam_inst_Z>>hBeamZSample" + str(i) + "(100, 0, 50)",
         "beam_backtrack == " + str(i)) 
  hBeamZSamples.append(RT.gDirectory.Get("hBeamZSample" + str(i)))
  hBeamZSamples[-1].SetFillColor(colors[i-1])
  hBeamZSamples[-1].SetLineColor(colors[i-1])
  hBeamZSamples[-1].Scale(scale)
  sBeamZSamples.Add(hBeamZSamples[-1])

  lBeamBT.AddEntry(hBeamXSamples[-1], tBeamBTs[i-1], "f")

  t.Draw("reco_beam_trackDirX*beam_inst_dirX + reco_beam_trackDirY*beam_inst_dirY + reco_beam_trackDirZ*beam_inst_dirZ>>hBeamCosSample" + str(i) + "(100, 0, 1)",
         "beam_backtrack == " + str(i)) 
  hBeamCosSamples.append(RT.gDirectory.Get("hBeamCosSample" + str(i)))
  hBeamCosSamples[-1].SetFillColor(colors[i-1])
  hBeamCosSamples[-1].SetLineColor(colors[i-1])
  hBeamCosSamples[-1].Scale(scale)
  sBeamCosSamples.Add(hBeamCosSamples[-1])

  hFlow = RT.TH1D("hFlow2" + str(i), "", 3, 0, 3)
  hFlow.SetBinContent(1, t.GetEntries("beam_backtrack == " + str(i)))
  hFlow.SetBinContent(2, t.GetEntries("selection_ID < 6 && beam_backtrack == " + str(i)))
  hFlow.SetBinContent(3, t.GetEntries("selection_ID < 5 && beam_backtrack == " + str(i)))
  hFlow.SetFillColor(colors[i-1])
  hFlow.SetLineColor(colors[i-1])
  hFlow.Scale(scale)
  sFlows.Add(hFlow)

  t.Draw("reco_beam_PFP_trackScore_collection>>hBeamTrackScores" + str(i) + "(102, -.02, 1)", "reco_beam_PFP_ID != -999 && reco_beam_PFP_trackScore_collection > 0. && beam_backtrack == " + str(i)) 
  hBeamTrackScores.append(RT.gDirectory.Get("hBeamTrackScores"+ str(i)))
  hBeamTrackScores[-1].SetBinContent(1, t.GetEntries("reco_beam_PFP_ID == -999 && beam_backtrack == " + str(i)))
  hBeamTrackScores[-1].SetBinContent(2, t.GetEntries("reco_beam_PFP_trackScore_collection < 0. && reco_beam_PFP_ID != -999 && beam_backtrack == " + str(i)))
  hBeamTrackScores[-1].SetFillColor(colors[i-1])
  hBeamTrackScores[-1].SetLineColor(colors[i-1])
  hBeamTrackScores[-1].Scale(scale)
  sBeamTrackScores.Add(hBeamTrackScores[-1])
 
tData.Draw("reco_beam_PFP_trackScore_collection>>hDataBeamTrackScores(102, -.02, 1)", "reco_beam_PFP_ID != -999 && reco_beam_PFP_trackScore_collection > 0.")
hDataBeamTrackScores = RT.gDirectory.Get("hDataBeamTrackScores")
hDataBeamTrackScores.SetBinContent(1, tData.GetEntries("reco_beam_PFP_ID == -999"))
hDataBeamTrackScores.SetBinContent(2, tData.GetEntries("reco_beam_PFP_trackScore_collection < 0. && reco_beam_PFP_ID != -999"))
hDataBeamTrackScores.SetMarkerColor(RT.kBlack)
hDataBeamTrackScores.SetLineColor(RT.kBlack)
hDataBeamTrackScores.SetMarkerStyle(20)
cBeamTrackScores = RT.TCanvas("cBeamTrackScores", "cBeamTrackScores")
cBeamTrackScores.SetTicks()
hDataBeamTrackScores.Draw("pe1")
sBeamTrackScores.Draw("hist same")
hDataBeamTrackScores.Draw("pe1 same")
lBeamBT.Draw()
cBeamTrackScores.Write()


fout.cd()
tData.Draw("reco_beam_calo_startX - beam_inst_X>>hDataBeamX(30, -30, 30)")
hDataBeamX = RT.gDirectory.Get("hDataBeamX")
hDataBeamX.SetMarkerColor(RT.kBlack)
hDataBeamX.SetLineColor(RT.kBlack)
hDataBeamX.SetMarkerStyle(20)
cBeamXSamples = RT.TCanvas("cBeamX", "cBeamX")
cBeamXSamples.SetTicks()
sBeamXSamples.Draw("hist")
if (sBeamXSamples.GetHistogram().GetMaximum() > hDataBeamX.GetMaximum()):
  sBeamXSamples.Draw("hist")
  sBeamXSamples.GetHistogram().SetTitle(";Beam #DeltaX (cm)")
  hDataBeamX.Draw("same pe1")
else:
  hDataBeamX.Draw("pe1")
  sBeamXSamples.Draw("same hist")
  hDataBeamX.SetTitle(";Beam #DeltaX (cm)")
  hDataBeamX.Draw("same pe1")
lBeamBT.AddEntry(hDataBeamX, "Data", "lp")
lBeamBT.Draw("same")

lLow = RT.TLine(0., 0., 0., max([sBeamXSamples.GetHistogram().GetMaximum(), hDataBeamX.GetMaximum()]))
lHigh = RT.TLine(10., 0., 10., max([sBeamXSamples.GetHistogram().GetMaximum(), hDataBeamX.GetMaximum()]))
#lLow.Draw("same")
#lHigh.Draw("same")
draw_vertical_line(0., RT.kBlack)
draw_vertical_line(10., RT.kBlack)

lLowMC = RT.TLine(-2., 0., -2., max([sBeamXSamples.GetHistogram().GetMaximum(), hDataBeamX.GetMaximum()]))
lHighMC = RT.TLine(2., 0., 2., max([sBeamXSamples.GetHistogram().GetMaximum(), hDataBeamX.GetMaximum()]))
#lLowMC.SetLineColor(RT.kRed)
#lLowMC.Draw("same")
#lHighMC.SetLineColor(RT.kRed)
#lHighMC.Draw("same")
draw_vertical_line(-2., RT.kRed)
draw_vertical_line(2., RT.kRed)

cBeamXSamples.Write()

tData.Draw("reco_beam_calo_startY - beam_inst_Y>>hDataBeamY(120, -30, 30)")
hDataBeamY = RT.gDirectory.Get("hDataBeamY")
hDataBeamY.SetMarkerColor(RT.kBlack)
hDataBeamY.SetLineColor(RT.kBlack)
hDataBeamY.SetMarkerStyle(20)
cBeamYSamples = RT.TCanvas("cBeamY", "cBeamY")
cBeamYSamples.SetTicks()
#sBeamYSamples.Draw("hist")
#sBeamYSamples.GetHistogram().SetTitle(";Beam #DeltaY (cm)")
#print(sBeamYSamples.GetHistogram().GetMaximum(), hDataBeamY.GetMaximum())
hDataBeamY.SetTitle(";Beam #DeltaY (cm)")
hDataBeamY.Draw("pe1")
sBeamYSamples.Draw("same hist")
hDataBeamY.Draw("pe1 same")
lBeamBT.Draw("same")

lLow = RT.TLine(-5., 0., -5., max([sBeamYSamples.GetHistogram().GetMaximum(), hDataBeamY.GetMaximum()]))
lHigh = RT.TLine(10., 0., 10., max([sBeamYSamples.GetHistogram().GetMaximum(), hDataBeamY.GetMaximum()]))
#lLow.Draw("same")
#lHigh.Draw("same")
draw_vertical_line(-5., RT.kBlack)
draw_vertical_line(10., RT.kBlack)

lLowMC = RT.TLine(-1.5, 0., -1.5, max([sBeamYSamples.GetHistogram().GetMaximum(), hDataBeamY.GetMaximum()]))
lHighMC = RT.TLine(2., 0., 2., max([sBeamYSamples.GetHistogram().GetMaximum(), hDataBeamY.GetMaximum()]))
#lLowMC.SetLineColor(RT.kRed)
#lLowMC.Draw("same")
#lHighMC.SetLineColor(RT.kRed)
#lHighMC.Draw("same")
draw_vertical_line(-1.5, RT.kRed)
draw_vertical_line(2., RT.kRed)

cBeamYSamples.Write()

tData.Draw("reco_beam_calo_startZ - beam_inst_Z>>hDataBeamZ(30, 0, 50)")
hDataBeamZ = RT.gDirectory.Get("hDataBeamZ")
hDataBeamZ.SetMarkerColor(RT.kBlack)
hDataBeamZ.SetLineColor(RT.kBlack)
hDataBeamZ.SetMarkerStyle(20)
cBeamZSamples = RT.TCanvas("cBeamZ", "cBeamZ")
cBeamZSamples.SetTicks()
sBeamZSamples.Draw("hist")
sBeamZSamples.GetHistogram().SetTitle(";Beam #DeltaZ (cm)")
hDataBeamZ.SetTitle(";Beam #DeltaZ (cm)")
hDataBeamZ.Draw("pe1")
sBeamZSamples.Draw("hist same")
hDataBeamZ.Draw("same pe1")
lBeamBT.Draw("same")

lLow = RT.TLine(30., 0., 30., max([sBeamZSamples.GetHistogram().GetMaximum(), hDataBeamZ.GetMaximum()]))
lHigh = RT.TLine(35., 0., 35., max([sBeamZSamples.GetHistogram().GetMaximum(), hDataBeamZ.GetMaximum()]))
#lLow.Draw("same")
#lHigh.Draw("same")
draw_vertical_line(30., RT.kBlack)
draw_vertical_line(35., RT.kBlack)

lLowMC = RT.TLine(28.5, 0., 28.5, max([sBeamZSamples.GetHistogram().GetMaximum(), hDataBeamZ.GetMaximum()]))
lHighMC = RT.TLine(31., 0., 31., max([sBeamZSamples.GetHistogram().GetMaximum(), hDataBeamZ.GetMaximum()]))
#lLowMC.SetLineColor(RT.kRed)
#lLowMC.Draw("same")
#lHighMC.SetLineColor(RT.kRed)
#lHighMC.Draw("same")
draw_vertical_line(28.5, RT.kRed)
draw_vertical_line(31., RT.kRed)
cBeamZSamples.Write()

tData.Draw("reco_beam_trackDirX*beam_inst_dirX + reco_beam_trackDirY*beam_inst_dirY + reco_beam_trackDirZ*beam_inst_dirZ>>hDataBeamCos(100, 0, 1)")
hDataBeamCos = RT.gDirectory.Get("hDataBeamCos")
hDataBeamCos.SetMarkerColor(RT.kBlack)
hDataBeamCos.SetLineColor(RT.kBlack)
hDataBeamCos.SetMarkerStyle(20)
cBeamCosSamples = RT.TCanvas("cBeamCos", "cBeamCos")
cBeamCosSamples.SetTicks()
sBeamCosSamples.Draw("hist")
sBeamCosSamples.GetHistogram().SetTitle(";Beam cos#theta")
sBeamCosSamples.GetHistogram().SetMaximum(max([sBeamCosSamples.GetHistogram().GetMaximum(), hDataBeamCos.GetMaximum()]))
sBeamCosSamples.Draw("hist")
hDataBeamCos.Draw("same pe1")
lBeamBT.Draw("same")
#lLow = RT.TLine(.93, 0., .93, max([sBeamCosSamples.GetHistogram().GetMaximum(), hDataBeamCos.GetMaximum()]))
#lLow.SetLineColor(RT.kBlack)
#lLow.Draw("same")
#lLowMC = RT.TLine(.97, 0., .97, max([sBeamCosSamples.GetHistogram().GetMaximum(), hDataBeamCos.GetMaximum()]))
#lLowMC.SetLineColor(RT.kRed)
#lLowMC.Draw("same")
draw_vertical_line(.93, RT.kBlack)
draw_vertical_line(.97, RT.kRed)
cBeamCosSamples.Write()

sBeamCosSamples.GetHistogram().GetXaxis().SetRangeUser(.9, 1.)
sBeamCosSamples.Draw("hist")
hDataBeamCos.Draw("same pe1")
lBeamBT.Draw("same")
draw_vertical_line(.3, RT.kBlack)
draw_vertical_line(.7, RT.kRed)
cBeamCosSamples.Write("cBeamCosSamples_zoom")


hTrackScores = []
sTrackScores = RT.THStack()

hMichelScores = []
sMichelScores = RT.THStack()

lTrackScores = RT.TLegend()
tTrackScores = [
  "Self", "Cosmic", "#pi", "#mu", "p", "#gamma", "Nucleus",
  "Daughter+", "Daughter++", "Michel", "#pi^{0} #gamma", "Other"]

hTruncatedMeans = []
sTruncatedMeans = RT.THStack()

hChi2PreTruncs = []
sChi2PreTruncs = RT.THStack()

hChi2PostTruncs = []
sChi2PostTruncs = RT.THStack()

for i in range(1, 13):
  t.Draw("reco_daughter_PFP_trackScore_collection>>hTrackScore" + str(i) + "(50, 0, 1)",
         "daughter_categories == " + str(i) + " && selection_ID < 4")
  hTrackScores.append(RT.gDirectory.Get("hTrackScore" + str(i)))
  hTrackScores[-1].SetFillColor(colors[i-1])
  hTrackScores[-1].SetLineColor(colors[i-1])
  hTrackScores[-1].Scale(scale)
  sTrackScores.Add(hTrackScores[-1])
  lTrackScores.AddEntry(hTrackScores[-1], tTrackScores[i-1], "f")

  t.Draw("reco_daughter_PFP_michelScore_collection>>hMichelScore" + str(i) + "(25, 0, 1)",
         "daughter_categories == " + str(i) + " && selection_ID < 4 && reco_daughter_PFP_trackScore_collection < .3")
  hMichelScores.append(RT.gDirectory.Get("hMichelScore" + str(i)))
  hMichelScores[-1].SetFillColor(colors[i-1])
  hMichelScores[-1].SetLineColor(colors[i-1])
  hMichelScores[-1].Scale(scale)
  sMichelScores.Add(hMichelScores[-1])

  t.Draw("reco_daughter_allTrack_truncLibo_dEdX_pos>>hTruncatedMean1" + str(i) + "(100, 0, 10)",
         "daughter_categories == " + str(i) + " && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")
  hTruncatedMeans.append(RT.gDirectory.Get("hTruncatedMean1" + str(i)))
  hTruncatedMeans[-1].SetFillColor(colors[i-1])
  hTruncatedMeans[-1].SetLineColor(colors[i-1])
  hTruncatedMeans[-1].Scale(scale)
  hTruncatedMeans[-1].AddBinContent(1, hTruncatedMeans[-1].GetBinContent(0))
  sTruncatedMeans.Add(hTruncatedMeans[-1])

  t.Draw("reco_daughter_allTrack_Chi2_proton/reco_daughter_allTrack_Chi2_ndof>>hChi2PreTrunc1" + str(i) + "(100, 0, 500)",
         "daughter_categories == " + str(i) + " && reco_daughter_allTrack_Chi2_ndof > 0 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")
  hChi2PreTruncs.append(RT.gDirectory.Get("hChi2PreTrunc1" + str(i)))
  hChi2PreTruncs[-1].SetFillColor(colors[i-1])
  hChi2PreTruncs[-1].SetLineColor(colors[i-1])
  hChi2PreTruncs[-1].Scale(scale)
  sChi2PreTruncs.Add(hChi2PreTruncs[-1])

  t.Draw("reco_daughter_allTrack_Chi2_proton/reco_daughter_allTrack_Chi2_ndof>>hChi2PostTrunc1" + str(i) + "(25, 0, 500)",
         "daughter_categories == " + str(i) + " && reco_daughter_allTrack_Chi2_ndof > 0 && reco_daughter_PFP_trackScore_collection > .3 && ((reco_daughter_allTrack_truncLibo_dEdX_pos > 2.8 && reco_daughter_allTrack_truncLibo_dEdX_pos < 3.4) || reco_daughter_allTrack_truncLibo_dEdX_pos < 0.5) && selection_ID < 4")
  hChi2PostTruncs.append(RT.gDirectory.Get("hChi2PostTrunc1" + str(i)))
  hChi2PostTruncs[-1].SetFillColor(colors[i-1])
  hChi2PostTruncs[-1].SetLineColor(colors[i-1])
  hChi2PostTruncs[-1].Scale(scale)
  sChi2PostTruncs.Add(hChi2PostTruncs[-1])


tData.Draw("reco_daughter_PFP_trackScore_collection>>hDataTrackScore(50, 0, 1)",
           "selection_ID < 4")
hDataTrackScore = RT.gDirectory.Get("hDataTrackScore")
hDataTrackScore.SetMarkerStyle(20)
hDataTrackScore.SetMarkerColor(RT.kBlack)
hDataTrackScore.SetLineColor(RT.kBlack)

cTrackScores = RT.TCanvas("cTrackScores", "cTrackScores")
cTrackScores.SetTicks()
sTrackScores.Draw("hist")
sTrackScores.GetHistogram().SetTitle(";CNN Track Score");
hDataTrackScore.Draw("pe1 same")
lTrackScores.Draw("same") ## add data to leg
draw_vertical_line(.3, RT.kBlack)
#lLow = RT.TLine(.3, 0., .3, sTrackScores.GetHistogram().GetMaximum())
#lLow.SetLineColor(RT.kBlack)
#lLow.Draw("same")
cTrackScores.Write()


cMichelScores = RT.TCanvas("cMichelScores", "cMichelScores")
cMichelScores.SetTicks()
sMichelScores.Draw("hist")
hDataMichelScore.Draw("pe1 same")
lTrackScores.Draw("same") ## add data to leg
cMichelScores.Write()

tData.Draw("reco_daughter_allTrack_truncLibo_dEdX_pos>>hDataTruncatedMean(100, 0, 10)",
           "reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")
hDataTruncatedMean = RT.gDirectory.Get("hDataTruncatedMean")
hDataTruncatedMean.SetMarkerStyle(20)
hDataTruncatedMean.SetMarkerColor(RT.kBlack)
hDataTruncatedMean.SetLineColor(RT.kBlack)
hDataTruncatedMean.AddBinContent(1, hDataTruncatedMean.GetBinContent(0))

cTruncatedMeans = RT.TCanvas("cTruncatedMeans1", "cTruncatedMeans1")
cTruncatedMeans.SetTicks()
sTruncatedMeans.Draw("hist")
sTruncatedMeans.GetHistogram().SetTitle(";Truncated Mean dE/dX (MeV/cm)")
hDataTruncatedMean.Draw("same pe1")
lTrackScores.Draw("same")
cTruncatedMeans.Write()

tData.Draw("reco_daughter_allTrack_Chi2_proton/reco_daughter_allTrack_Chi2_ndof>>hDataChi2PreTrunc(100, 0, 500)",
           "reco_daughter_allTrack_Chi2_ndof > 0 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")
hDataChi2PreTrunc = RT.gDirectory.Get("hDataChi2PreTrunc")
hDataChi2PreTrunc.SetMarkerStyle(20)
hDataChi2PreTrunc.SetMarkerColor(RT.kBlack)
hDataChi2PreTrunc.SetLineColor(RT.kBlack)

tData.Draw("reco_daughter_allTrack_Chi2_proton/reco_daughter_allTrack_Chi2_ndof>>hDataChi2PostTrunc(25, 0, 500)",
           "reco_daughter_allTrack_Chi2_ndof > 0 && reco_daughter_PFP_trackScore_collection > .3 && ((reco_daughter_allTrack_truncLibo_dEdX_pos > 2.8 && reco_daughter_allTrack_truncLibo_dEdX_pos < 3.4) || reco_daughter_allTrack_truncLibo_dEdX_pos < 0.5) && selection_ID < 4")
hDataChi2PostTrunc = RT.gDirectory.Get("hDataChi2PostTrunc")
hDataChi2PostTrunc.SetMarkerStyle(20)
hDataChi2PostTrunc.SetMarkerColor(RT.kBlack)
hDataChi2PostTrunc.SetLineColor(RT.kBlack)

cChi2PreTruncs = RT.TCanvas("cChi2PreTruncs1", "cChi2PreTruncs1")
cChi2PreTruncs.SetTicks()
sChi2PreTruncs.Draw("hist")
hDataChi2PreTrunc.Draw("pe1 same")
lTrackScores.AddEntry(hDataChi2PreTrunc, "Data", "lp")
lTrackScores.Draw("same")
cChi2PreTruncs.Write()

cChi2PostTruncs = RT.TCanvas("cChi2PostTruncs1", "cChi2PostTruncs1")
cChi2PostTruncs.SetTicks()
sChi2PostTruncs.Draw("hist")
sChi2PostTruncs.GetHistogram().SetTitle(";PID #chi^{2}/ndof")
hDataChi2PostTrunc.Draw("same pe1")
lTrackScores.Draw("same")
cChi2PostTruncs.Write()



hTruncatedMeans = []
sTruncatedMeans = RT.THStack()

hTrackScores = []
sTrackScores = RT.THStack()

lTrackScores = RT.TLegend()
tTrackScores = [
  "#pi", "#mu", "p", "#gamma",
  "Nucleus", "e", "Other"
]


hChi2PreTruncs = []
sChi2PreTruncs = RT.THStack()

hChi2PostTruncs = []
sChi2PostTruncs = RT.THStack()

hShowerEnergies = []
sShowerEnergies = RT.THStack()
hShowerDists = []
sShowerDists = RT.THStack()

hShowerDists = []
sShowerDists = RT.THStack()

hMichelScorePDGs = []
sMichelScorePDGs = RT.THStack()

for i in range(1, 8):
  t.Draw("reco_daughter_allTrack_truncLibo_dEdX_pos>>hTruncatedMean" + str(i) + "(100, 0, 10)",
         "daughter_PDGs_types == " + str(i) + " && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")
  hTruncatedMeans.append(RT.gDirectory.Get("hTruncatedMean" + str(i)))
  hTruncatedMeans[-1].SetFillColor(colors[i-1])
  hTruncatedMeans[-1].SetLineColor(colors[i-1])
  hTruncatedMeans[-1].Scale(scale)
  hTruncatedMeans[-1].AddBinContent(1, hTruncatedMeans[-1].GetBinContent(0))
  sTruncatedMeans.Add(hTruncatedMeans[-1])

  t.Draw("reco_daughter_PFP_trackScore_collection>>hTrackScorePDG" + str(i) + "(50, 0, 1)",
         "daughter_PDGs_types == " + str(i))
  hTrackScores.append(RT.gDirectory.Get("hTrackScorePDG" + str(i)))
  hTrackScores[-1].SetFillColor(colors[i-1])
  hTrackScores[-1].SetLineColor(colors[i-1])
  hTrackScores[-1].Scale(scale)
  sTrackScores.Add(hTrackScores[-1])
  lTrackScores.AddEntry(hTrackScores[-1], tTrackScores[i-1], "f")

  t.Draw("reco_daughter_allTrack_Chi2_proton/reco_daughter_allTrack_Chi2_ndof>>hChi2PreTrunc" + str(i) + "(100, 0, 500)",
         "daughter_PDGs_types == " + str(i) + " && reco_daughter_allTrack_Chi2_ndof > 0 && reco_daughter_PFP_trackScore_collection > .3 && selection_ID < 4")
  hChi2PreTruncs.append(RT.gDirectory.Get("hChi2PreTrunc" + str(i)))
  hChi2PreTruncs[-1].SetFillColor(colors[i-1])
  hChi2PreTruncs[-1].SetLineColor(colors[i-1])
  hChi2PreTruncs[-1].Scale(scale)
  sChi2PreTruncs.Add(hChi2PreTruncs[-1])

  t.Draw("reco_daughter_allTrack_Chi2_proton/reco_daughter_allTrack_Chi2_ndof>>hChi2PostTrunc" + str(i) + "(25, 0, 500)",
         "daughter_PDGs_types == " + str(i) + " && reco_daughter_allTrack_Chi2_ndof > 0 && reco_daughter_PFP_trackScore_collection > .3 && reco_daughter_allTrack_truncLibo_dEdX_pos > 2.8 && reco_daughter_allTrack_truncLibo_dEdX_pos < 3.4 && selection_ID < 4")
  hChi2PostTruncs.append(RT.gDirectory.Get("hChi2PostTrunc" + str(i)))
  hChi2PostTruncs[-1].SetFillColor(colors[i-1])
  hChi2PostTruncs[-1].SetLineColor(colors[i-1])
  hChi2PostTruncs[-1].Scale(scale)
  sChi2PostTruncs.Add(hChi2PostTruncs[-1])

  t.Draw("reco_daughter_allShower_energy>>hShowerEnergy" + str(i) + "(100, 0, 500)",
         "daughter_PDGs_types == " + str(i) + " && reco_daughter_PFP_trackScore_collection > 0. && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4")
  hShowerEnergies.append(RT.gDirectory.Get("hShowerEnergy" + str(i)))
  hShowerEnergies[-1].SetFillColor(colors[i-1])
  hShowerEnergies[-1].SetLineColor(colors[i-1])
  hShowerEnergies[-1].Scale(scale)
  sShowerEnergies.Add(hShowerEnergies[-1])

  t.Draw("shower_dists>>hShowerDist" + str(i) + "(100, 0, 500)",
         "daughter_PDGs_types == " + str(i) + " && reco_daughter_PFP_trackScore_collection > 0. && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4 && reco_daughter_allShower_ID > 0")
  hShowerDists.append(RT.gDirectory.Get("hShowerDist" + str(i)))
  hShowerDists[-1].SetFillColor(colors[i-1])
  hShowerDists[-1].SetLineColor(colors[i-1])
  hShowerDists[-1].Scale(scale)
  sShowerDists.Add(hShowerDists[-1])

  t.Draw("reco_daughter_PFP_michelScore_collection>>hMichelScorePDG" + str(i) + "(25, 0, 1)",
         "daughter_PDGs_types == " + str(i) + " && selection_ID < 4 && reco_daughter_PFP_trackScore_collection < .3")
  hMichelScorePDGs.append(RT.gDirectory.Get("hMichelScorePDG" + str(i)))
  hMichelScorePDGs[-1].SetFillColor(colors[i-1])
  hMichelScorePDGs[-1].SetLineColor(colors[i-1])
  hMichelScorePDGs[-1].Scale(scale)
  sMichelScorePDGs.Add(hMichelScorePDGs[-1])

cTruncatedMeans = RT.TCanvas("cTruncatedMeans", "cTruncatedMeans")
cTruncatedMeans.SetTicks()
sTruncatedMeans.Draw("hist")
sTruncatedMeans.GetHistogram().SetTitle(";Truncated Mean dE/dX (MeV/cm)")
hDataTruncatedMean.Draw("same pe1")
lTrackScores.Draw("same")
cTruncatedMeans.Write()

cTrackScores = RT.TCanvas("cTrackScoresPDG", "cTrackScoresPDG")
cTrackScores.SetTicks()
sTrackScores.Draw("hist")
sTrackScores.GetHistogram().SetTitle(";CNN Track Score")
hDataTrackScore.Draw("pe1 same")
lTrackScores.Draw("same")
draw_vertical_line(.3, RT.kBlack)
cTrackScores.Write()


cMichelScorePDGs = RT.TCanvas("cMichelScoresPDG", "cMichelScoresPDG")
cMichelScorePDGs.SetTicks()
sMichelScorePDGs.Draw("hist")
lTrackScores.Draw("same")
cMichelScorePDGs.Write()


cChi2PreTruncs = RT.TCanvas("cChi2PreTruncs", "cChi2PreTruncs")
cChi2PreTruncs.SetTicks()
sChi2PreTruncs.Draw("hist")
hDataChi2PreTrunc.Draw("pe1 same")
lTrackScores.AddEntry(hDataChi2PreTrunc, "Data", "lp")
lTrackScores.Draw("same")
cChi2PreTruncs.Write()

cChi2PostTruncs = RT.TCanvas("cChi2PostTruncs", "cChi2PostTruncs")
cChi2PostTruncs.SetTicks()
sChi2PostTruncs.Draw("hist")
sChi2PostTruncs.GetHistogram().SetTitle(";PID #chi^{2}/ndof")
hDataChi2PostTrunc.Draw("same pe1")
lTrackScores.Draw("same")
cChi2PostTruncs.Write()


tData.Draw("reco_daughter_allShower_energy>>hDataShowerEnergy(100, 0, 500)",
           "reco_daughter_PFP_trackScore_collection > 0. && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4")
hDataShowerEnergy = RT.gDirectory.Get("hDataShowerEnergy")
hDataShowerEnergy.SetMarkerStyle(20)
hDataShowerEnergy.SetMarkerColor(RT.kBlack)
hDataShowerEnergy.SetLineColor(RT.kBlack)

cShowerEnergies = RT.TCanvas("cShowerEnergies", "cShowerEnergies")
cShowerEnergies.SetTicks()
sShowerEnergies.Draw("hist")
sShowerEnergies.GetHistogram().SetTitle(";Shower Energy (MeV)")
hDataShowerEnergy.Draw("same pe1")
lTrackScores.Draw("same")
cShowerEnergies.Write()

tData.Draw("shower_dists>>hDataShowerDist(100, 0, 500)",
         "reco_daughter_PFP_trackScore_collection > 0. && reco_daughter_PFP_trackScore_collection < .3 && selection_ID < 4 && reco_daughter_allShower_ID > 0")
hDataShowerDist = RT.gDirectory.Get("hDataShowerDist")
hDataShowerDist.SetMarkerStyle(20)
hDataShowerDist.SetMarkerColor(RT.kBlack)
hDataShowerDist.SetLineColor(RT.kBlack)

cShowerDists = RT.TCanvas("cShowerDists", "cShowerDists")
cShowerDists.SetTicks()
sShowerDists.Draw("hist")
sShowerDists.GetHistogram().SetTitle(";Shower Distance (cm)")
hDataShowerDist.Draw("same pe1")
lTrackScores.Draw("same")
cShowerDists.Write()




hBeamTrackScores = []
sBeamTrackScores = RT.THStack()
hDataBeamTrackScoresList = []
sDataBeamTrackScores = RT.THStack()

beam_types = [13, -999, 11]
lBeamTypes = RT.TLegend()
titles = ["Track", "None", "Shower"]
for i in range(0, 3):
  t.Draw("reco_beam_PFP_trackScore_collection>>hBeamTrackScores2" + str(i) + "(101, -.01, 1)", "reco_beam_PFP_ID != -999 && reco_beam_type == " + str(beam_types[i])) 
  hBeamTrackScores.append(RT.gDirectory.Get("hBeamTrackScores2"+ str(i)))
  hBeamTrackScores[-1].SetBinContent(1, t.GetEntries("reco_beam_PFP_ID == -999 && reco_beam_type == " + str(beam_types[i])))
  hBeamTrackScores[-1].SetFillColor(colors[i-1])
  hBeamTrackScores[-1].SetLineColor(colors[i-1])
  hBeamTrackScores[-1].Scale(scale)
  sBeamTrackScores.Add(hBeamTrackScores[-1])

  tData.Draw("reco_beam_PFP_trackScore_collection>>hDataBeamTrackScores2" + str(i) + "(101, -.01, 1)", "reco_beam_PFP_ID != -999 && reco_beam_type == " + str(beam_types[i])) 
  hDataBeamTrackScoresList.append(RT.gDirectory.Get("hDataBeamTrackScores2"+ str(i)))
  hDataBeamTrackScoresList[-1].SetBinContent(1, tData.GetEntries("reco_beam_PFP_ID == -999 && reco_beam_type == " + str(beam_types[i])))
  hDataBeamTrackScoresList[-1].SetFillColor(colors[i-1])
  hDataBeamTrackScoresList[-1].SetLineColor(colors[i-1])
  sDataBeamTrackScores.Add(hDataBeamTrackScoresList[-1])

  lBeamTypes.AddEntry(hBeamTrackScores[-1], titles[i], "f")

cBeamTrackScoresBeamTypes = RT.TCanvas("cBeamTrackScoresBeamTypes", "cBeamTrackScoresBeamTypes")
cBeamTrackScoresBeamTypes.SetTicks()
hDataBeamTrackScores.Draw("pe1")
hDataBeamTrackScores.SetTitle(";Track Score")
sBeamTrackScores.Draw("hist same")
hDataBeamTrackScores.Draw("pe1 same")
lBeamTypes.Draw()
cBeamTrackScoresBeamTypes.Write()

sDataBeamTrackScores.Write("sDataBeamTrackScores")

fout.Close()
