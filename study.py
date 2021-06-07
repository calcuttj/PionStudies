import ROOT as RT#from ROOT import * 
import sys
from defcuts import defcuts, testcuts, testcuts_FS, ang_pos_test_cut, ang_cut_str, pos_cut_str
from array import array



RT.gROOT.SetBatch(1)

f = RT.TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

base_cut = " && reco_beam_type == 13 && (true_beam_PDG == -13 || true_beam_PDG == 211)"

fout = RT.TFile( sys.argv[2], "RECREATE" )

cuts = defcuts()

lenhists = dict()
endZhists = dict()
endYhists = dict()
startZhists = dict()
startXhists = dict()
startYhists = dict()
trackDirZhists = dict()
trackDirXhists = dict()
trackDirYhists = dict()
caloDirZhists = dict()
caloDirXhists = dict()
caloDirYhists = dict()
caloZhists = dict()
caloXhists = dict()
caloYhists = dict()
deltaXhists = dict()
deltaYhists = dict()
coshists = dict()
cosSCEhists = dict()
beamXhists = dict()
beamYhists = dict()
beam_dirXhists = dict()
beam_dirYhists = dict()
beam_dirZhists = dict()
chi2hists   = dict()
cnnhists   = dict()
cnn_collectionhists = dict()

names = [
  "PrimaryPion",
  "PrimaryMuon",
  #"PrimaryProton",
  #"PrimaryElectron",
  "Cosmic",
  "PrimaryBeamNotTrig",
  #"UpstreamIntToPiPlus",
  #"UpstreamIntToPiMinus",
  "UpstreamIntToPion",
  "UpstreamIntToProton",
  "UpstreamIntToKaon",
  "UpstreamIntToNuc",
  #"PionInelToPiPlus",
  #"PionInelToPiMinus",
  #"PionInelToProton",
  #"PionInelToKaon",
  #"PionInelToNuc",
  #"UpstreamInt",
  "Decay",
  #"Other"
]

title_map = {
  "PrimaryPion": "Primary #pi^{+}",
  "PrimaryMuon": "Primary #mu",
  "Cosmic":"Cosmic",
  "PrimaryBeamNotTrig": "Extra Primary",
  #"UpstreamIntToPiPlus": "#pi^{+} From Interaction",
  #"UpstreamIntToPiMinus": "#pi^{-} From Interaction",
  "UpstreamIntToPion": "#pi^{#pm} From Interaction",
  "UpstreamIntToProton": "p From Interaction",
  "UpstreamIntToKaon": "K From Interaction",
  "UpstreamIntToNuc": "Nuc From Interaction",
  "Decay":"Decay Product"
}

tree.Draw( "reco_beam_len>>lenhist(40,0,500.)", "1 " + base_cut )
lenhist = RT.gDirectory.Get("lenhist")
#print("Total:", lenhist.Integral())
#for name,cut in cuts.iteritems():

cut_total = 0
for name in names:
  cut = cuts[name]
  ##print name, cut

  tree.Draw( "reco_beam_len>>len"+name+"(40,0.,500.)", cut  + base_cut )
  lenhists[name] = RT.gDirectory.Get("len"+name)

  ##print(name+":",lenhists[name].Integral())
  cut_total = cut_total + lenhists[name].Integral()

  tree.Draw( "reco_beam_endZ>>endZ"+name+"(40,0.,500.)", cut  + base_cut )
  endZhists[name] = RT.gDirectory.Get("endZ"+name)

  tree.Draw( "reco_beam_endY>>endY"+name+"(40,0.,500.)", cut  + base_cut )
  endYhists[name] = RT.gDirectory.Get("endY"+name)

  tree.Draw( "reco_beam_startX>>startX"+name+"(40,-100.,100.)", cut  + base_cut )
  startXhists[name] = RT.gDirectory.Get("startX"+name)

  tree.Draw( "reco_beam_startY>>startY"+name+"(40,380.,500.)", cut  + base_cut )
  startYhists[name] = RT.gDirectory.Get("startY"+name)

  tree.Draw( "reco_beam_startZ>>startZ"+name+"(50,0.,50.)", cut  + base_cut )
  startZhists[name] = RT.gDirectory.Get("startZ"+name)

  tree.Draw( "reco_beam_trackDirX>>trackDirX"+name+"(40, -.75, .5)", cut  + base_cut )
  trackDirXhists[name] = RT.gDirectory.Get("trackDirX"+name)

  tree.Draw( "reco_beam_trackDirY>>trackDirY"+name+"(40, -.75, .25)", cut  + base_cut )
  trackDirYhists[name] = RT.gDirectory.Get("trackDirY"+name)

  tree.Draw( "reco_beam_trackDirZ>>trackDirZ"+name+"(20, .75, 1.)", cut  + base_cut )
  trackDirZhists[name] = RT.gDirectory.Get("trackDirZ"+name)

  tree.Draw( "reco_beam_calo_startDirX[0]>>caloDirX"+name+"(40, -.75, .5)", cut  + base_cut )
  caloDirXhists[name] = RT.gDirectory.Get("caloDirX"+name)

  tree.Draw( "reco_beam_calo_startDirY[0]>>caloDirY"+name+"(40, -.75, .25)", cut  + base_cut )
  caloDirYhists[name] = RT.gDirectory.Get("caloDirY"+name)

  tree.Draw( "reco_beam_calo_startDirZ[0]>>caloDirZ"+name+"(20, .75, 1.)", cut  + base_cut )
  caloDirZhists[name] = RT.gDirectory.Get("caloDirZ"+name)

  tree.Draw( "reco_beam_calo_startX[0]>>caloX"+name+"(40, -100., 100.)", cut  + base_cut )
  caloXhists[name] = RT.gDirectory.Get("caloX"+name)

  tree.Draw( "reco_beam_calo_startY[0]>>caloY"+name+"(40, 380., 500.)", cut  + base_cut )
  caloYhists[name] = RT.gDirectory.Get("caloY"+name)

  tree.Draw( "reco_beam_calo_startZ[0]>>caloZ"+name+"(20, -20., 20.)", cut  + base_cut )
  caloZhists[name] = RT.gDirectory.Get("caloZ"+name)

  tree.Draw( "data_BI_X>>beamX"+name+"(40,-100.,100.)", cut  + base_cut )
  beamXhists[name] = RT.gDirectory.Get("beamX"+name)

  tree.Draw( "data_BI_Y>>beamY"+name+"(40,380.,500.)", cut  + base_cut )
  beamYhists[name] = RT.gDirectory.Get("beamY"+name)

  tree.Draw( "data_BI_dirX>>beam_dirX"+name+"(100, -.3, -.1)", cut  + base_cut )
  beam_dirXhists[name] = RT.gDirectory.Get("beam_dirX"+name)

  tree.Draw( "data_BI_dirY>>beam_dirY"+name+"(100, -.3, -.1)", cut  + base_cut )
  beam_dirYhists[name] = RT.gDirectory.Get("beam_dirY"+name)

  tree.Draw( "data_BI_dirZ>>beam_dirZ"+name+"(100, .9, 1.)", cut  + base_cut )
  beam_dirZhists[name] = RT.gDirectory.Get("beam_dirZ"+name)

  #tree.Draw( "reco_beam_startX - (true_beam_startX + -1.*true_beam_startZ*(true_beam_startDirX/true_beam_startDirZ))>>deltaX"+name+"(50,-100.,100.)", cut  + base_cut )
  #deltaXhists[name] = RT.gDirectory.Get("deltaX"+name)

  #tree.Draw( "reco_beam_startY - (true_beam_startY + -1.*true_beam_startZ*(true_beam_startDirY/true_beam_startDirZ))>>deltaY"+name+"(50,-100.,100.)", cut  + base_cut )
  #deltaYhists[name] = RT.gDirectory.Get("deltaY"+name)

  tree.Draw( "reco_beam_startY - data_BI_Y>>deltaY"+name+"(50,-50.,50.)", cut  + base_cut )
  deltaYhists[name] = RT.gDirectory.Get("deltaY"+name)

  tree.Draw( "reco_beam_startX - data_BI_X>>deltaX"+name+"(50,-50.,50.)", cut  + base_cut )
  deltaXhists[name] = RT.gDirectory.Get("deltaX"+name)

  tree.Draw( "(reco_beam_trackDirX*data_BI_dirX + reco_beam_trackDirY*data_BI_dirY + reco_beam_trackDirZ*data_BI_dirZ)>>cos"+name+"(50, .75, 1.)", cut  + base_cut )
  coshists[name] = RT.gDirectory.Get("cos"+name)

  tree.Draw( "(reco_beam_calo_startDirX[0]*data_BI_dirX + reco_beam_calo_startDirY[0]*data_BI_dirY + reco_beam_calo_startDirZ[0]*data_BI_dirZ)>>cosSCE"+name+"(50, .75, 1.)", cut  + base_cut )
  cosSCEhists[name] = RT.gDirectory.Get("cosSCE"+name)

 # tree.Draw( "(reco_beam_calo_startDirX[1]*data_BI_dirX + reco_beam_calo_startDirY[1]*data_BI_dirY + reco_beam_calo_startDirZ[1]*data_BI_dirZ)>>cosSCE1"+name+"(50, .75, 1.)", cut  + base_cut )
 # cosSCE1hists[name] = RT.gDirectory.Get("cosSCE1"+name)

 # tree.Draw( "(reco_beam_calo_startDirX[2]*data_BI_dirX + reco_beam_calo_startDirY[2]*data_BI_dirY + reco_beam_calo_startDirZ[2]*data_BI_dirZ)>>cosSCE2"+name+"(50, .75, 1.)", cut  + base_cut )
 # cosSCE2hists[name] = RT.gDirectory.Get("cosSCE2"+name)

  tree.Draw( "reco_beam_Chi2_proton/reco_beam_Chi2_ndof>>chi2"+name+"(100,0,400)", cut + base_cut)
  chi2hists[name] = RT.gDirectory.Get("chi2"+name)

  tree.Draw( "reco_beam_PFP_trackScore>>cnn"+name+"(100,0,1.)", cut + base_cut)
  cnnhists[name] = RT.gDirectory.Get("cnn"+name)

  tree.Draw( "reco_beam_PFP_trackScore_collection>>cnn_collection"+name+"(100,0,1.)", cut + " && ( true_beam_PDG == 211 || true_beam_PDG == -13)")
  cnn_collectionhists[name] = RT.gDirectory.Get("cnn_collection"+name)


  fout.cd()
  #lenhists[name].Write()
  ##print


#print(cut_total)

colors = {
  "PrimaryMuon": (RT.kOrange+1),
  "PrimaryPion": (RT.kBlue-4),
  #"PrimaryProton": (RT.kTeal+2),
  #"PrimaryElectron": (RT.kViolet-7),
  "PrimaryBeamNotTrig": (RT.kRed+2),
  "Cosmic": (RT.kSpring-8),
  #"UpstreamIntToPiPlus": (RT.kViolet-3),
  #"UpstreamIntToPiMinus": (RT.kCyan-2),
  "UpstreamIntToPion": (RT.kViolet-3),
  "UpstreamIntToProton": (RT.kTeal),
  "UpstreamIntToKaon": (RT.kYellow-6),
  "UpstreamIntToNuc": (RT.kOrange+10),
  #"PionInelToPiPlus",
  #"PionInelToPiMinus",
  #"PionInelToProton",
  #"PionInelToKaon",
  #"PionInelToNuc",
  #"UpstreamInt",
  "Other":(RT.kMagenta-10),
  "Decay":(RT.kOrange-7)
}

leg = RT.TLegend(.15,.5,.45,.8)


lenstack = RT.THStack("lenstack","")
endZstack = RT.THStack("endZstack","")
endYstack = RT.THStack("endYstack","")
startXstack = RT.THStack("startXstack","")
startYstack = RT.THStack("startYstack","")
startZstack = RT.THStack("startZstack","")
trackDirXstack = RT.THStack("trackDirXstack","")
trackDirYstack = RT.THStack("trackDirYstack","")
trackDirZstack = RT.THStack("trackDirZstack","")
caloDirXstack = RT.THStack("caloDirXstack","")
caloDirYstack = RT.THStack("caloDirYstack","")
caloDirZstack = RT.THStack("caloDirZstack","")
caloXstack = RT.THStack("caloXstack","")
caloYstack = RT.THStack("caloYstack","")
caloZstack = RT.THStack("caloZstack","")
deltaXstack = RT.THStack("deltaXstack","")
deltaYstack = RT.THStack("deltaYstack","")
cosstack = RT.THStack("cosstack","")
cosSCEstack = RT.THStack("cosSCEstack","")
beamXstack = RT.THStack("beamXstack","")
beamYstack = RT.THStack("beamYstack","")
beam_dirXstack = RT.THStack("beam_dirXstack","")
beam_dirYstack = RT.THStack("beam_dirYstack","")
beam_dirZstack = RT.THStack("beam_dirZstack","")
chi2stack   = RT.THStack("chi2stack","")
cnnstack   = RT.THStack("cnnstack","")
cnn_collectionstack   = RT.THStack("cnn_collectionstack","")
for name in names:


  lenhist = lenhists[name]
  lenhist.SetFillColor(colors[name]) 
  lenhist.SetLineColor(colors[name]) 
  lenstack.Add(lenhist)

  endZhist = endZhists[name]
  endZhist.SetFillColor(colors[name]) 
  endZhist.SetLineColor(colors[name]) 
  endZstack.Add(endZhist)

  endYhist = endYhists[name]
  endYhist.SetFillColor(colors[name]) 
  endYhist.SetLineColor(colors[name]) 
  endYstack.Add(endYhist)

  startXhist = startXhists[name]
  startXhist.SetFillColor(colors[name]) 
  startXhist.SetLineColor(colors[name]) 
  startXstack.Add(startXhist)

  startYhist = startYhists[name]
  startYhist.SetFillColor(colors[name]) 
  startYhist.SetLineColor(colors[name]) 
  startYstack.Add(startYhist)

  startZhist = startZhists[name]
  startZhist.SetFillColor(colors[name]) 
  startZhist.SetLineColor(colors[name]) 
  startZstack.Add(startZhist)
  leg.AddEntry( startZhist, title_map[name], "f")

  trackDirXhist = trackDirXhists[name]
  trackDirXhist.SetFillColor(colors[name]) 
  trackDirXhist.SetLineColor(colors[name]) 
  trackDirXstack.Add(trackDirXhist)

  trackDirYhist = trackDirYhists[name]
  trackDirYhist.SetFillColor(colors[name]) 
  trackDirYhist.SetLineColor(colors[name]) 
  trackDirYstack.Add(trackDirYhist)

  trackDirZhist = trackDirZhists[name]
  trackDirZhist.SetFillColor(colors[name]) 
  trackDirZhist.SetLineColor(colors[name]) 
  trackDirZstack.Add(trackDirZhist)

  caloDirXhist = caloDirXhists[name]
  caloDirXhist.SetFillColor(colors[name]) 
  caloDirXhist.SetLineColor(colors[name]) 
  caloDirXstack.Add(caloDirXhist)

  caloDirYhist = caloDirYhists[name]
  caloDirYhist.SetFillColor(colors[name]) 
  caloDirYhist.SetLineColor(colors[name]) 
  caloDirYstack.Add(caloDirYhist)

  caloDirZhist = caloDirZhists[name]
  caloDirZhist.SetFillColor(colors[name]) 
  caloDirZhist.SetLineColor(colors[name]) 
  caloDirZstack.Add(caloDirZhist)

  caloXhist = caloXhists[name]
  caloXhist.SetFillColor(colors[name]) 
  caloXhist.SetLineColor(colors[name]) 
  caloXstack.Add(caloXhist)

  caloYhist = caloYhists[name]
  caloYhist.SetFillColor(colors[name]) 
  caloYhist.SetLineColor(colors[name]) 
  caloYstack.Add(caloYhist)

  caloZhist = caloZhists[name]
  caloZhist.SetFillColor(colors[name]) 
  caloZhist.SetLineColor(colors[name]) 
  caloZstack.Add(caloZhist)

  deltaXhist = deltaXhists[name]
  deltaXhist.SetFillColor(colors[name]) 
  deltaXhist.SetLineColor(colors[name]) 
  deltaXstack.Add(deltaXhist)

  deltaYhist = deltaYhists[name]
  deltaYhist.SetFillColor(colors[name]) 
  deltaYhist.SetLineColor(colors[name]) 
  deltaYstack.Add(deltaYhist)

  coshist = coshists[name]
  coshist.SetFillColor(colors[name]) 
  coshist.SetLineColor(colors[name]) 
  cosstack.Add(coshist)

  cosSCEhist = cosSCEhists[name]
  cosSCEhist.SetFillColor(colors[name]) 
  cosSCEhist.SetLineColor(colors[name]) 
  cosSCEstack.Add(cosSCEhist)

  chi2hist = chi2hists[name]
  chi2hist.SetFillColor(colors[name]) 
  chi2hist.SetLineColor(colors[name]) 
  chi2stack.Add(chi2hist)

  cnnhist = cnnhists[name]
  cnnhist.SetFillColor(colors[name]) 
  cnnhist.SetLineColor(colors[name]) 
  cnnstack.Add(cnnhist)

  cnn_collectionhist = cnn_collectionhists[name]
  cnn_collectionhist.SetFillColor(colors[name]) 
  cnn_collectionhist.SetLineColor(colors[name]) 
  cnn_collectionstack.Add(cnn_collectionhist)

  beamXhist = beamXhists[name]
  beamXhist.SetFillColor(colors[name]) 
  beamXhist.SetLineColor(colors[name]) 
  beamXstack.Add(beamXhist)

  beamYhist = beamYhists[name]
  beamYhist.SetFillColor(colors[name]) 
  beamYhist.SetLineColor(colors[name]) 
  beamYstack.Add(beamYhist)

  beam_dirXhist = beam_dirXhists[name]
  beam_dirXhist.SetFillColor(colors[name]) 
  beam_dirXhist.SetLineColor(colors[name]) 
  beam_dirXstack.Add(beam_dirXhist)

  beam_dirYhist = beam_dirYhists[name]
  beam_dirYhist.SetFillColor(colors[name]) 
  beam_dirYhist.SetLineColor(colors[name]) 
  beam_dirYstack.Add(beam_dirYhist)

  beam_dirZhist = beam_dirZhists[name]
  beam_dirZhist.SetFillColor(colors[name]) 
  beam_dirZhist.SetLineColor(colors[name]) 
  beam_dirZstack.Add(beam_dirZhist)

leg.Write("leg")

lenstack.Write()
endZstack.Write()
endYstack.Write()
startZstack.Write()
startYstack.Write()
startXstack.Write()
trackDirZstack.Write()
trackDirYstack.Write()
trackDirXstack.Write()
caloDirZstack.Write()
caloDirYstack.Write()
caloDirXstack.Write()
caloZstack.Write()
caloYstack.Write()
caloXstack.Write()
deltaYstack.Write()
deltaXstack.Write()
cosstack.Write()
cosSCEstack.Write()
beamYstack.Write()
beamXstack.Write()
beam_dirYstack.Write()
beam_dirZstack.Write()
beam_dirXstack.Write()
chi2stack.Write()
cnn_collectionstack.Write()
cnnstack.Write()



#### looking at the final state of the primary pions ####


colors = {
  "PrimaryPionDecay": (RT.kTeal+3),
  #"PrimaryPionFastScint": (RT.kRed-4),
  "PrimaryPionInteract": (RT.kBlue-4),
  #"PrimaryPion": (RT.kBlue-4),
  "PrimaryMuon": (RT.kOrange+1),
  #"PrimaryProton": (RT.kTeal+2),
  #"PrimaryElectron": (RT.kViolet-7),
  "PrimaryBeamNotTrig": (RT.kRed+2),
  "Cosmic": (RT.kSpring-8),
  #"UpstreamIntToPiPlus": (RT.kViolet-3),
  #"UpstreamIntToPiMinus": (RT.kCyan-2),
  "UpstreamIntToPion": (RT.kViolet-3),
  "UpstreamIntToProton": (RT.kTeal),
  "UpstreamIntToKaon": (RT.kYellow-6),
  "UpstreamIntToNuc": (RT.kOrange+10),
  #"PionInelToPiPlus",
  #"PionInelToPiMinus",
  #"PionInelToProton",
  #"PionInelToKaon",
  #"PionInelToNuc",
  #"UpstreamInt",
  "Other":(RT.kMagenta-10),
  "Decay":(RT.kOrange-7)
}

leg = RT.TLegend(.6,.6,.85,.85)
names = [
  "PrimaryPionDecay",
  #"PrimaryPionFastScint",
  "PrimaryPionInteract",
  #"PrimaryPion",
  "PrimaryMuon",
  #"PrimaryProton",
  #"PrimaryElectron",
  "Cosmic",
  "PrimaryBeamNotTrig",
  #"UpstreamIntToPiPlus",
  #"UpstreamIntToPiMinus",
  "UpstreamIntToPion",
  "UpstreamIntToProton",
  "UpstreamIntToKaon",
  "UpstreamIntToNuc",
  #"PionInelToPiPlus",
  #"PionInelToPiMinus",
  #"PionInelToProton",
  #"PionInelToKaon",
  #"PionInelToNuc",
  #"UpstreamInt",
  "Decay",
  #"Other"
]

titles = [
  #"Primary #pi^{+} Decay",
  #"Primary #pi^{+} Decay at Rest",
  #"Primary #pi^{+} Interacting",
  "Primary #pi^{+}",
  "Primary #mu",
  #"PrimaryProton",
  #"PrimaryElectron",
  "Cosmic",
  "Primary Beam Not Trig",
  "\"Upstream\" Int #rightarrow #pi^{+}",
  "\"Upstream\" Int #rightarrow #pi^{-}",
  "\"Upstream\" Int #rightarrow p",
  "\"Upstream\" Int #rightarrow K",
  "\"Upstream\" Int #rightarrow Nucleus",
#  "Pion #rightarrow #pi^{+}",
#  "Pion #rightarrow #pi^{-}",
#  "Pion #rightarrow p",
#  "Pion #rightarrow K",
#  "Pion #rightarrow Nucleus",
#  "Pion #rightarrow #gamma",
#  "\"Upstream\" Interaction",
  "Decay Product",
  #"Other"
]

lenhists = dict()
endZhists = dict()
endYhists = dict()

endP_hists = dict()
endP_stack = RT.THStack("endPstack","")

for name in names:
  cut = cuts[name]
  #print name, cut

  tree.Draw( "reco_beam_len>>len_FS_"+name+"(40,0.,500.)", cut  + base_cut )
  lenhists[name] = RT.gDirectory.Get("len_FS_"+name)
  fout.cd()
  #lenhists[name].Write()
  #print

  tree.Draw( "reco_beam_true_byHits_endP>>endP_"+name+"(40,0.,1.2)", cut  + base_cut)
  endP_hists[name] = RT.gDirectory.Get("endP_"+name)
  #endP_hists[name].Write()

  endP_hists[name].SetFillColor(colors[name]) 
  endP_hists[name].SetLineColor(colors[name]) 
  endP_stack.Add(endP_hists[name])

  tree.Draw( "reco_beam_endZ>>endZ_FS_"+name+"(40,0.,500.)", cut  + base_cut )
  endZhists[name] = RT.gDirectory.Get("endZ_FS_"+name)

  tree.Draw( "reco_beam_endY>>endY_FS_"+name+"(40,0.,500.)", cut  + base_cut )
  endYhists[name] = RT.gDirectory.Get("endY_FS_"+name)


endP_stack.Write()

leg = RT.TLegend(.6,.6,.85,.85)
lenstack = RT.THStack("lenstack_FS","")
endZstack = RT.THStack("endZstack_FS","")
endYstack = RT.THStack("endYstack_FS","")

for name,title in zip(names,titles):

  lenhist = lenhists[name]
  lenhist.SetFillColor(colors[name]) 
  lenhist.SetLineColor(colors[name]) 
  lenstack.Add(lenhist)
  leg.AddEntry( lenhist, title, "f")

  endZhist = endZhists[name]
  endZhist.SetFillColor(colors[name]) 
  endZhist.SetLineColor(colors[name]) 
  endZstack.Add(endZhist)

  endYhist = endYhists[name]
  endYhist.SetFillColor(colors[name]) 
  endYhist.SetLineColor(colors[name]) 
  endYstack.Add(endYhist)

leg.Write("leg_FS")

lenstack.Write()
endZstack.Write()
endYstack.Write()


endhists = dict()

for name in names:
  cut = cuts[name]
  #print name, cut

  tree.Draw( "true_beam_endZ:reco_beam_endZ>>reco_beam_endZ"+name+"(40,0.,500.,40,0.,500.)", cut  + base_cut, "colz" )
  endhists[name] = RT.gDirectory.Get("endZ"+name)
  fout.cd()
  #endhists[name].Write()
  #print



##### Z pos #####

z_pos_dir = fout.mkdir( "z_pos_dir", "track z position" )
z_pos_dir.cd()

zhists = dict()
zstack = RT.THStack("zstack","")

for name in names:
  cut = cuts[name]
  #print name, cut
  tree.Draw( "reco_beam_startZ>>startZ_"+name+"(50,0.,50.)", cut  + base_cut )

  zhists[name] = RT.gDirectory.Get("startZ_"+name)
  #zhists[name].Write()

  zhists[name].SetFillColor(colors[name]) 
  zhists[name].SetLineColor(colors[name]) 
  zstack.Add(zhists[name])

  #print

zstack.Write()
#################



### Defining angular and position cuts ###
#ang_cut = " && (true_beam_Start_DirX*trackDirX + true_beam_Start_DirY*trackDirY + true_beam_Start_DirZ*trackDirZ > .93)"
#
#pos_cut = " && ( (true_beam_Start_X + -1.*true_beam_Start_Z*(true_beam_Start_DirX/true_beam_Start_DirZ) - startX) > -4. )"
#pos_cut = pos_cut + " && ( (true_beam_Start_X + -1.*true_beam_Start_Z*(true_beam_Start_DirX/true_beam_Start_DirZ) - startX) < 4. )"
#pos_cut = pos_cut + " && ( (true_beam_Start_Y + -1.*true_beam_Start_Z*(true_beam_Start_DirY/true_beam_Start_DirZ) - startY) > -7. )"
#pos_cut = pos_cut + " && ( (true_beam_Start_Y + -1.*true_beam_Start_Z*(true_beam_Start_DirY/true_beam_Start_DirZ) - startY) < 8. )"
#pos_cut = pos_cut + " && ( startZ > 16. && startZ < 20.)"

ang_cut = ang_cut_str()
pos_cut = pos_cut_str()
##########################################


#def ang_pos_test_cut(e):
#  if (e.true_beam_Start_DirX*e.trackDirX + e.true_beam_Start_DirY*e.trackDirY + e.true_beam_Start_DirZ*e.trackDirZ < .93): return 0
#
#  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) < -4. ): return 0
#
#  if ( (e.true_beam_Start_X + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirX/e.true_beam_Start_DirZ) - e.startX) > 4. ): return 0
#
#  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) < -7. ): return 0
#
#  if ( (e.true_beam_Start_Y + -1.*e.true_beam_Start_Z*(e.true_beam_Start_DirY/e.true_beam_Start_DirZ) - e.startY) > 8. ): return 0
#
#  if( e.startZ < 16. or e.startZ > 20. ): return 0
#
#  return 1



### First cut: with angular and position cuts ###
first_cut_dir = fout.mkdir( "first_cut_dir", "Cuts include start position and angular cuts")
first_cut_dir.cd()

lenhists = dict()
lenstack = RT.THStack("lenstack_ang_pos","")

endP_hists = dict()
endP_stack = RT.THStack("endPstack_ang_pos","")

endZ_hists = dict()
endZ_stack = RT.THStack("endZstack_ang_pos","")

endY_hists = dict()
endY_stack = RT.THStack("endYstack_ang_pos","")

chi2hists = dict()
chi2stack = RT.THStack("chi2stack_ang_pos","")

cnnhists = dict()
cnnstack = RT.THStack("cnnstack_ang_pos","")



for name in names:
  cut = cuts[name]
  #print name, cut

  tree.Draw( "reco_beam_len>>len_ang_pos_cut_"+name+"(40,0.,500.)", cut  + base_cut + ang_cut + pos_cut)
  lenhists[name] = RT.gDirectory.Get("len_ang_pos_cut_"+name)
  #lenhists[name].Write()

  lenhists[name].SetFillColor(colors[name]) 
  lenhists[name].SetLineColor(colors[name]) 
  lenstack.Add(lenhists[name])

  tree.Draw( "reco_beam_true_byHits_endP>>endP_ang_pos_cut_"+name+"(40,0.,1.2)", cut  + base_cut + ang_cut + pos_cut)
  endP_hists[name] = RT.gDirectory.Get("endP_ang_pos_cut_"+name)
  #endP_hists[name].Write()

  endP_hists[name].SetFillColor(colors[name]) 
  endP_hists[name].SetLineColor(colors[name]) 
  endP_stack.Add(endP_hists[name])

  tree.Draw( "reco_beam_endZ>>endZ_ang_pos_cut_"+name+"(40,0.,500.)", cut  + base_cut + ang_cut + pos_cut)
  endZ_hists[name] = RT.gDirectory.Get("endZ_ang_pos_cut_"+name)
  #endZ_hists[name].Write()

  endZ_hists[name].SetFillColor(colors[name]) 
  endZ_hists[name].SetLineColor(colors[name]) 
  endZ_stack.Add(endZ_hists[name])

  tree.Draw( "reco_beam_endY>>endY_ang_pos_cut_"+name+"(40,0.,500.)", cut  + base_cut + ang_cut + pos_cut)
  endY_hists[name] = RT.gDirectory.Get("endY_ang_pos_cut_"+name)
  #endY_hists[name].Write()

  endY_hists[name].SetFillColor(colors[name]) 
  endY_hists[name].SetLineColor(colors[name]) 
  endY_stack.Add(endY_hists[name])

  tree.Draw( "reco_beam_Chi2_proton/reco_beam_Chi2_ndof>>chi2_ang_pos_cut_"+name+"(100,0,400)", cut + base_cut + ang_cut + pos_cut)
  chi2hists[name] = RT.gDirectory.Get("chi2_ang_pos_cut_"+name)
  chi2hists[name].SetFillColor(colors[name]) 
  chi2hists[name].SetLineColor(colors[name]) 
  chi2stack.Add(chi2hists[name])

  tree.Draw( "reco_beam_PFP_trackScore>>cnn_ang_pos_cut_"+name+"(100,0,1.)", cut + base_cut + ang_cut + pos_cut)
  cnnhists[name] = RT.gDirectory.Get("cnn_ang_pos_cut_"+name)
  cnnhists[name].SetFillColor(colors[name]) 
  cnnhists[name].SetLineColor(colors[name]) 
  cnnstack.Add(cnnhists[name])


  #print

lenstack.Write()
endP_stack.Write()
endZ_stack.Write()
endY_stack.Write()
chi2stack.Write()
cnnstack.Write()
#endP_first_cut_signal = endP_hists["PrimaryPionInteract"].Clone("first_cut_signal")
#endP_first_cut_signal.Sumw2()
#endP_first_cut_signal.SetFillColor(0)
#endP_first_cut_signal.Divide(endP_total_signal)
#endP_first_cut_signal.Write("endP_first_cut_efficiency")
###################################################





## Now with length cuts ###
second_cut_dir = fout.mkdir( "second_cut_dir", "Cuts include start position and angular cuts and track length cut")
second_cut_dir.cd()

lenhists = dict()
lenstack = RT.THStack("lenstack_ang_pos_len","")

endP_hists = dict()
endP_stack = RT.THStack("endPstack_ang_pos_len","")

endZ_hists = dict()
endZ_stack = RT.THStack("endZstack_ang_pos_len","")

endY_hists = dict()
endY_stack = RT.THStack("endYstack_ang_pos_len","")

chi2hists = dict()
chi2stack = RT.THStack("chi2stack_ang_pos_len","")

cnnhists = dict()
cnnstack = RT.THStack("cnnstack_ang_pos_len","")

endZ_cut = " && reco_beam_endZ < 226."

for name in names:
  cut = cuts[name]
  #print name, cut

  tree.Draw( "reco_beam_len>>len_ang_pos_endZ_cut_"+name+"(40,0.,500.)", cut  + base_cut + ang_cut + pos_cut + endZ_cut)
  lenhists[name] = RT.gDirectory.Get("len_ang_pos_endZ_cut_"+name)
  #lenhists[name].Write()

  lenhists[name].SetFillColor(colors[name]) 
  lenhists[name].SetLineColor(colors[name]) 
  lenstack.Add(lenhists[name])

  tree.Draw( "reco_beam_true_byHits_endP>>endP_ang_pos_endZ_cut_"+name+"(40,0.,1.2)", cut  + base_cut + ang_cut + pos_cut + endZ_cut)
  endP_hists[name] = RT.gDirectory.Get("endP_ang_pos_endZ_cut_"+name)
  #endP_hists[name].Write()

  endP_hists[name].SetFillColor(colors[name]) 
  endP_hists[name].SetLineColor(colors[name]) 
  endP_stack.Add(endP_hists[name])

  tree.Draw( "reco_beam_endZ>>endZ_ang_pos_endZ_cut_"+name+"(40,0.,500.)", cut  + base_cut + ang_cut + pos_cut + endZ_cut)
  endZ_hists[name] = RT.gDirectory.Get("endZ_ang_pos_endZ_cut_"+name)
  #endZ_hists[name].Write()

  endZ_hists[name].SetFillColor(colors[name]) 
  endZ_hists[name].SetLineColor(colors[name]) 
  endZ_stack.Add(endZ_hists[name])

  tree.Draw( "reco_beam_endY>>endY_ang_pos_endZ_cut_"+name+"(40,0.,500.)", cut  + base_cut + ang_cut + pos_cut + endZ_cut)
  endY_hists[name] = RT.gDirectory.Get("endY_ang_pos_endZ_cut_"+name)
  #endY_hists[name].Write()

  endY_hists[name].SetFillColor(colors[name]) 
  endY_hists[name].SetLineColor(colors[name]) 
  endY_stack.Add(endY_hists[name])

  tree.Draw( "reco_beam_Chi2_proton/reco_beam_Chi2_ndof>>chi2_ang_pos_endZ_cut_"+name+"(100,0,400)", cut + base_cut + ang_cut + pos_cut + endZ_cut)
  chi2hists[name] = RT.gDirectory.Get("chi2_ang_pos_endZ_cut_"+name)
  chi2hists[name].SetFillColor(colors[name]) 
  chi2hists[name].SetLineColor(colors[name]) 
  chi2stack.Add(chi2hists[name])

  tree.Draw( "reco_beam_PFP_trackScore>>cnn_ang_pos_endZ_cut_"+name+"(100,0,1.)", cut + base_cut + ang_cut + pos_cut + endZ_cut)
  cnnhists[name] = RT.gDirectory.Get("cnn_ang_pos_endZ_cut_"+name)
  cnnhists[name].SetFillColor(colors[name]) 
  cnnhists[name].SetLineColor(colors[name]) 
  cnnstack.Add(cnnhists[name])


  #print

lenstack.Write()
endP_stack.Write()
endZ_stack.Write()
endY_stack.Write()
chi2stack.Write()
cnnstack.Write()


## Now with chi2 cuts ###
third_cut_dir = fout.mkdir( "third_cut_dir", "Cuts include start position and angular cuts and track length cut and chi2 cut")
third_cut_dir.cd()

lenhists = dict()
lenstack = RT.THStack("lenstack_ang_pos_chi2","")

endP_hists = dict()
endP_stack = RT.THStack("endPstack_ang_pos_chi2","")

endZ_hists = dict()
endZ_stack = RT.THStack("endZstack_ang_pos_chi2","")

chi2hists = dict()
chi2stack = RT.THStack("chi2stack_ang_pos_chi2","")

cnnhists = dict()
cnnstack = RT.THStack("cnnstack_ang_pos_chi2","")

endZ_cut = " && reco_beam_endZ < 226."

chi2_cut = " && reco_beam_Chi2_proton/reco_beam_Chi2_ndof > 140. "

for name in names:
  cut = cuts[name]
  #print name, cut

  tree.Draw( "reco_beam_len>>len_ang_pos_endZ_chi2_cut_"+name+"(40,0.,500.)", cut  + base_cut + ang_cut + pos_cut + endZ_cut + chi2_cut)
  lenhists[name] = RT.gDirectory.Get("len_ang_pos_endZ_chi2_cut_"+name)
  #lenhists[name].Write()

  lenhists[name].SetFillColor(colors[name]) 
  lenhists[name].SetLineColor(colors[name]) 
  lenstack.Add(lenhists[name])

  tree.Draw( "reco_beam_true_byHits_endP>>endP_ang_pos_endZ_chi2_cut_"+name+"(40,0.,1.2)", cut  + base_cut + ang_cut + pos_cut + endZ_cut + chi2_cut)
  endP_hists[name] = RT.gDirectory.Get("endP_ang_pos_endZ_chi2_cut_"+name)
  #endP_hists[name].Write()

  endP_hists[name].SetFillColor(colors[name]) 
  endP_hists[name].SetLineColor(colors[name]) 
  endP_stack.Add(endP_hists[name])

  tree.Draw( "reco_beam_endZ>>endZ_ang_pos_endZ_chi2_cut_"+name+"(40,0.,500.)", cut  + base_cut + ang_cut + pos_cut + endZ_cut + chi2_cut)
  endZ_hists[name] = RT.gDirectory.Get("endZ_ang_pos_endZ_chi2_cut_"+name)
  #endZ_hists[name].Write()

  endZ_hists[name].SetFillColor(colors[name]) 
  endZ_hists[name].SetLineColor(colors[name]) 
  endZ_stack.Add(endZ_hists[name])

  tree.Draw( "reco_beam_Chi2_proton/reco_beam_Chi2_ndof>>chi2_ang_pos_endZ_chi2_cut_"+name+"(100,0,400)", cut + base_cut + ang_cut + pos_cut + endZ_cut + chi2_cut)
  chi2hists[name] = RT.gDirectory.Get("chi2_ang_pos_endZ_chi2_cut_"+name)
  chi2hists[name].SetFillColor(colors[name]) 
  chi2hists[name].SetLineColor(colors[name]) 
  chi2stack.Add(chi2hists[name])

  tree.Draw( "reco_beam_PFP_trackScore>>cnn_ang_pos_endZ_chi2_cut_"+name+"(100,0,1.)", cut + base_cut + ang_cut + pos_cut + endZ_cut + chi2_cut)
  cnnhists[name] = RT.gDirectory.Get("cnn_ang_pos_endZ_chi2_cut_"+name)
  cnnhists[name].SetFillColor(colors[name]) 
  cnnhists[name].SetLineColor(colors[name]) 
  cnnstack.Add(cnnhists[name])


  #print

lenstack.Write()
endP_stack.Write()
endZ_stack.Write()
chi2stack.Write()
cnnstack.Write()


exit()


## Now with dedx cuts ###
third_cut_dir = fout.mkdir( "third_cut_dir", "Cuts include start position and angular cuts and dedx cut")
third_cut_dir.cd()

lenhists = dict()
lenstack = RT.THStack("lenstack_ang_pos_dedx","")

endP_hists = dict()
endP_stack = RT.THStack("endPstack_ang_pos_dedx","")

endZ_hists = dict()
endZ_stack = RT.THStack("endZstack_ang_pos_dedx","")

for name in names:
  #print name

  lenhists[name] = RT.TH1D("len_ang_pos_dedx_cut_"+name, "", 40, 0., 500.)
  lenhists[name].SetFillColor(colors[name]) 
  lenhists[name].SetLineColor(colors[name]) 

  endP_hists[name] = RT.TH1D("endP_ang_pos_dedx_cut_"+name, "", 40, 0., 1.2)
  endP_hists[name].SetFillColor(colors[name]) 
  endP_hists[name].SetLineColor(colors[name]) 

  endZ_hists[name] = RT.TH1D("endZ_ang_pos_dedx_cut_"+name, "", 40, 0., 1.2)
  endZ_hists[name].SetFillColor(colors[name]) 
  endZ_hists[name].SetLineColor(colors[name]) 

  #print

for e in tree:

  if len([i for i in e.reco_beam_dEdX]) < 1: continue
  if not ang_pos_test_cut(e): continue

  total_dedx = sum( [i for i in e.reco_beam_dEdX] )
  avg_dedx =  total_dedx /  len( [i for i in e.reco_beam_dEdX] ) 
  len_avg_dedx = total_dedx / e.reco_beam_len 

  if( avg_dedx > 5. ): continue

  cut = testcuts_FS(e)
  if cut == "bad" or cut == "PrimaryElectron" or cut == "PrimaryProton" or cut == "NeutronInel" or cut == "ProtonInel" or cut == "Other": continue

  lenhists[cut].Fill( e.reco_beam_len )
  endP_hists[cut].Fill( e.reco_beam_true_byE_endP )
  endZ_hists[cut].Fill(e.reco_beam_endZ)

for name in names:
  lenstack.Add( lenhists[name] )
  endP_stack.Add( endP_hists[name] )
  endZ_stack.Add( endZ_hists[name] )

  #lenhists[name].Write()
  #endP_hists[name].Write()

lenstack.Write()
endP_stack.Write()
endZ_stack.Write()
###########################

### ###
#fourth_cut_dir = fout.mkdir( "fourth_cut_dir", "Cuts include start position and angular cuts and track length cut and dedx cut")
#fourth_cut_dir.cd()
#
#lenhists = dict()
#lenstack = THStack("lenstack_ang_pos_dedx_len","")
#
#endP_hists = dict()
#endP_stack = THStack("endPstack_ang_pos_dedx_len","")
#
#for name in names:
#  #print name
#
#  lenhists[name] = TH1D("len_ang_pos_dedx_endZ_cut_"+name, "", 40, 0., 500.)
#  lenhists[name].SetFillColor(colors[name]) 
#  lenhists[name].SetLineColor(colors[name]) 
#
#  endP_hists[name] = TH1D("endP_ang_pos_dedx_endZ_cut_"+name, "", 40, 0., 1.2)
#  endP_hists[name].SetFillColor(colors[name]) 
#  endP_hists[name].SetLineColor(colors[name]) 
#
#  #print
#
#for e in tree:
#
#  if len([i for i in e.reco_beam_dEdX]) < 1: continue
#  if not ang_pos_test_cut(e): continue
#  if e.reco_beam_len > 250.: continue
#
#  total_dedx = sum( [i for i in e.reco_beam_dEdX] )
#  avg_dedx =  total_dedx /  len( [i for i in e.reco_beam_dEdX] ) 
#  len_avg_dedx = total_dedx / e.reco_beam_len 
#
#  if( avg_dedx > 5. ): continue
#
#  cut = testcuts_FS(e)
#  if cut == "bad" or cut == "PrimaryElectron" or cut == "PrimaryProton" or cut == "NeutronInel" or cut == "ProtonInel" or cut == "Other": continue
#
#  lenhists[cut].Fill( e.reco_beam_len )
#  endP_hists[cut].Fill( e.reco_beam_true_byE_endP )
#
#for name in names:
#  lenstack.Add( lenhists[name] )
#  endP_stack.Add( endP_hists[name] )
#
#  #lenhists[name].Write()
#  #endP_hists[name].Write()
#
#lenstack.Write()
#endP_stack.Write()
#
#endP_fourth_cut_signal = endP_hists["PrimaryPionInteract"].Clone("fourth_cut_signal")
#endP_fourth_cut_signal.Sumw2()
#endP_fourth_cut_signal.SetFillColor(0)
#endP_fourth_cut_signal.Divide(endP_total_signal)
#endP_fourth_cut_signal.Write("endP_fourth_cut_efficiency")
############################

