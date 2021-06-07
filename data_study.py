import ROOT as RT
import sys
from defcuts import defcuts, testcuts, testcuts_FS, ang_pos_test_cut, ang_cut_str, pos_cut_str
from defcuts import data_ang_pos_test_cut, data_ang_cut_str, data_pos_cut_str
from array import array
from set_style import *



RT.gROOT.SetBatch(1)

f = RT.TFile( sys.argv[1] )

tree = f.Get("pionana/beamana")

base_cut = " && reco_beam_type == 13"

fout = RT.TFile( sys.argv[2], "RECREATE" )
outtree = RT.TTree("tree", "")
endZ = array("d", [0])
length = array("d", [0])
startZ = array("d", [0])
startX = array("d", [0])
startY = array("d", [0])
beamZ = array("d", [0])
beamX = array("d", [0])
beamY = array("d", [0])
startDirZ = array("d", [0])
startDirX = array("d", [0])
startDirY = array("d", [0])
caloDirZ = array("d", [0])
caloDirX = array("d", [0])
caloDirY = array("d", [0])
caloZ = array("d", [0])
caloX = array("d", [0])
caloY = array("d", [0])
beamDirZ = array("d", [0])
beamDirX = array("d", [0])
beamDirY = array("d", [0])
chi2 = array("d", [0])
cnn = array("d", [0])
cnn_collection = array("d", [0])


outtree.Branch("endZ",        endZ,       "endZ/D"  )
outtree.Branch("length",      length,     "length/D")
outtree.Branch("startZ",      startZ,     "startZ/D")
outtree.Branch("startX",      startX,     "startX/D")
outtree.Branch("startY",      startY,     "startY/D")
outtree.Branch("beamZ",       beamZ,      "beamZ/D")
outtree.Branch("beamX",       beamX,      "beamX/D")
outtree.Branch("beamY",       beamY,      "beamY/D")
outtree.Branch("startDirZ",   startDirZ,  "startDirZ/D")
outtree.Branch("startDirX",   startDirX,  "startDirX/D")
outtree.Branch("startDirY",   startDirY,  "startDirY/D")
outtree.Branch("caloDirZ",   caloDirZ,  "caloDirZ/D")
outtree.Branch("caloDirX",   caloDirX,  "caloDirX/D")
outtree.Branch("caloDirY",   caloDirY,  "caloDirY/D")
outtree.Branch("caloZ",   caloZ,  "caloZ/D")
outtree.Branch("caloX",   caloX,  "caloX/D")
outtree.Branch("caloY",   caloY,  "caloY/D")
outtree.Branch("beamDirZ",    beamDirZ,   "beamDirZ/D")
outtree.Branch("beamDirX",    beamDirX,   "beamDirX/D")
outtree.Branch("beamDirY",    beamDirY,   "beamDirY/D")
outtree.Branch("chi2",    chi2,   "chi2/D")
outtree.Branch("cnn",    cnn,   "cnn/D")
outtree.Branch("cnn_collection",    cnn_collection,   "cnn_collection/D")





for e in tree:
  if not( e.reco_beam_type == 13 and ( 211 in [i for i in e.data_BI_PDG_candidates])): continue

  if not( e.data_BI_nMomenta == 1 and e.data_BI_nTracks == 1): continue

  endZ[0] = e.reco_beam_endZ
  length[0] = e.reco_beam_len
  startZ[0] = e.reco_beam_startZ
  startX[0] = e.reco_beam_startX
  startY[0] = e.reco_beam_startY

  endZ[0] = e.reco_beam_endZ

  beamZ[0] = e.data_BI_Z
  beamX[0] = e.data_BI_X
  beamY[0] = e.data_BI_Y

  beamDirZ[0] = e.data_BI_dirZ
  beamDirX[0] = e.data_BI_dirX
  beamDirY[0] = e.data_BI_dirY

  startDirZ[0] = e.reco_beam_trackDirZ
  startDirX[0] = e.reco_beam_trackDirX
  startDirY[0] = e.reco_beam_trackDirY

  if (len([i for i in e.reco_beam_calo_startDirZ])):
    caloDirZ[0] = e.reco_beam_calo_startDirZ[0]
    caloDirX[0] = e.reco_beam_calo_startDirX[0]
    caloDirY[0] = e.reco_beam_calo_startDirY[0]
  else:
    caloDirZ[0] = -1.
    caloDirX[0] = -1.
    caloDirY[0] = -1.

  caloZ[0] = e.reco_beam_calo_startZ
  caloX[0] = e.reco_beam_calo_startX
  caloY[0] = e.reco_beam_calo_startY

  chi2[0] = e.reco_beam_Chi2_proton / e.reco_beam_Chi2_ndof
  cnn[0] = e.reco_beam_PFP_trackScore
  cnn_collection[0] = e.reco_beam_PFP_trackScore_collection


  outtree.Fill()

outtree.Draw( "length>>lenhist(40,0,500.)" )
lenhist = RT.gDirectory.Get("lenhist")

outtree.Draw( "startX>>startXhist(40, -100., 100.)" )
startXhist = RT.gDirectory.Get("startXhist")

outtree.Draw( "startY>>startYhist(40, 380., 500.)" )
startYhist = RT.gDirectory.Get("startYhist")

outtree.Draw( "startZ>>startZhist(50, 0., 50.)" )
startZhist = RT.gDirectory.Get("startZhist")

outtree.Draw( "startX-beamX>>deltaXhist(50, -100., 100.)" )
deltaXhist = RT.gDirectory.Get("deltaXhist")

outtree.Draw( "startY-beamY>>deltaYhist(50, -100., 100.)" )
deltaYhist = RT.gDirectory.Get("deltaYhist")

outtree.Draw( "beamX>>beamXhist(40,-100.,100.)" )
beamXhist = RT.gDirectory.Get("beamXhist")

outtree.Draw( "beamY>>beamYhist(40, 380., 500.)" )
beamYhist = RT.gDirectory.Get("beamYhist")

outtree.Draw( "beamDirX>>beam_dirXhist(100, -.3, -.1)" )
beam_dirXhist = RT.gDirectory.Get("beam_dirXhist")

outtree.Draw( "beamDirY>>beam_dirYhist(100, -.3, -.1)" )
beam_dirYhist = RT.gDirectory.Get("beam_dirYhist")

outtree.Draw( "beamDirZ>>beam_dirZhist(100, .9, 1.)" )
beam_dirZhist = RT.gDirectory.Get("beam_dirZhist")

outtree.Draw( "startDirX>>trackDirXhist(40, -.75, .5)" )
trackDirXhist = RT.gDirectory.Get("trackDirXhist")

outtree.Draw( "startDirY>>trackDirYhist(40, -.75, .5)" )
trackDirYhist = RT.gDirectory.Get("trackDirYhist")

outtree.Draw( "startDirZ>>trackDirZhist(20, .75, 1.)" )
trackDirZhist = RT.gDirectory.Get("trackDirZhist")

outtree.Draw( "caloDirX>>caloDirXhist(40, -.75, .5)" )
caloDirXhist = RT.gDirectory.Get("caloDirXhist")

outtree.Draw( "caloDirY>>caloDirYhist(40, -.75, .5)" )
caloDirYhist = RT.gDirectory.Get("caloDirYhist")

outtree.Draw( "caloDirZ>>caloDirZhist(20, .75, 1.)" )
caloDirZhist = RT.gDirectory.Get("caloDirZhist")

outtree.Draw( "caloX>>caloXhist(40, -100., 100.)" )
caloXhist = RT.gDirectory.Get("caloXhist")

outtree.Draw( "caloY>>caloYhist(40, 380., 500.)" )
caloYhist = RT.gDirectory.Get("caloYhist")

outtree.Draw( "caloZ>>caloZhist(20, -20., 20.)" )
caloZhist = RT.gDirectory.Get("caloZhist")

outtree.Draw( "(beamDirX*startDirX + beamDirY*startDirY + beamDirZ*startDirZ)>>coshist(50, .75, 1.)")
coshist = RT.gDirectory.Get("coshist")

outtree.Draw( "(caloDirX*startDirX + caloDirY*startDirY + caloDirZ*startDirZ)>>cosSCEhist(50, .75, 1.)")
cosSCEhist = RT.gDirectory.Get("cosSCEhist")

outtree.Draw( "chi2>>chi2hist(100, 0., 400.)" )
chi2hist = RT.gDirectory.Get("chi2hist")

outtree.Draw( "cnn>>cnnhist(100, 0., 1.)" )
cnnhist = RT.gDirectory.Get("cnnhist")

outtree.Draw( "cnn_collection>>cnn_collectionhist(100, 0., 1.)" )
cnn_collectionhist = RT.gDirectory.Get("cnn_collectionhist")

outtree.Draw( "endZ>>endZhist(40, 0., 500.)" )
endZhist = RT.gDirectory.Get("endZhist")

fout.cd()

set_style(lenhist, "Track Length (cm)", "")
markers(lenhist)
lenhist.Write()

set_style(startXhist, "Track Start X (cm)", "")
markers(startXhist)
startXhist.Write()

set_style(startYhist, "Track Start Y (cm)", "")
markers(startYhist)
startYhist.Write()

set_style(startZhist, "Track Start Z (cm)", "")
markers(startZhist)
startZhist.Write()

set_style(deltaXhist, "#DeltaX (cm)", "")
markers(deltaXhist)
deltaXhist.Write()

set_style(deltaYhist, "#DeltaY (cm)", "")
markers(deltaYhist)
deltaYhist.Write()

set_style(coshist, "Cos(#theta)", "")
markers(coshist)
coshist.Write()

set_style(cosSCEhist, "Cos(#theta)", "")
markers(cosSCEhist)
cosSCEhist.Write()

set_style(chi2hist, "#chi^{2}", "")
markers(chi2hist)
chi2hist.Write()

set_style(cnnhist, "Track Score", "")
markers(cnnhist)
cnnhist.Write()

set_style(cnn_collectionhist, "Track Score", "")
markers(cnn_collectionhist)
cnn_collectionhist.Write()

set_style(endZhist, "Track End Z (cm)", "")
markers(endZhist)
endZhist.Write()

set_style(beamXhist, "Beam X (cm)", "")
markers(beamXhist)
beamXhist.Write()

set_style(beamYhist, "Beam Y (cm)", "")
markers(beamYhist)
beamYhist.Write()

set_style(beam_dirXhist, "Beam dir X", "")
markers(beam_dirXhist)
beam_dirXhist.Write()

set_style(beam_dirYhist, "Beam dir Y", "")
markers(beam_dirYhist)
beam_dirYhist.Write()

set_style(beam_dirZhist, "Beam dir Z", "")
markers(beam_dirZhist)
beam_dirZhist.Write()

set_style(trackDirXhist, "Reco dir X", "")
markers(trackDirXhist)
trackDirXhist.Write()

set_style(trackDirYhist, "Reco dir Y", "")
markers(trackDirYhist)
trackDirYhist.Write()

set_style(trackDirZhist, "Reco dir Z", "")
markers(trackDirZhist)
trackDirZhist.Write()

set_style(caloDirXhist, "Reco dir X", "")
markers(caloDirXhist)
caloDirXhist.Write()

set_style(caloDirYhist, "Reco dir Y", "")
markers(caloDirYhist)
caloDirYhist.Write()

set_style(caloDirZhist, "Reco dir Z", "")
markers(caloDirZhist)
caloDirZhist.Write()

set_style(caloXhist, "Reco X (cm)", "")
markers(caloXhist)
caloXhist.Write()

set_style(caloYhist, "Reco Y (cm)", "")
markers(caloYhist)
caloYhist.Write()

set_style(caloZhist, "Reco Z (cm)", "")
markers(caloZhist)
caloZhist.Write()

### First cut: with angular and position cuts ###
first_cut_dir = fout.mkdir( "first_cut_dir", "Cuts include start position and angular cuts")
first_cut_dir.cd()

ang_cut = " && ( ( startDirX*beamDirX + startDirY*beamDirY + startDirZ*beamDirZ ) > .93 )" 

pos_cut  = " && ( ( startX - beamX ) > 0. ) " 
pos_cut += " && ( ( startX - beamX ) < 10. ) " 
pos_cut += " && ( ( startY - beamY ) > -5. ) " 
pos_cut += " && ( ( startY - beamY ) < 10. ) " 
pos_cut += " && ( startZ > 30. ) && ( startZ < 35. )"

outtree.Draw( "length>>len_ang_pos_cut(40,0.,500.)", "1 " + ang_cut + pos_cut)
lenhist = RT.gDirectory.Get("len_ang_pos_cut")

outtree.Draw( "endZ>>endZ_ang_pos_cut(40,0.,500.)", "1 " + ang_cut + pos_cut)
endZhist = RT.gDirectory.Get("endZ_ang_pos_cut")

outtree.Draw( "startZ>>startZ_ang_pos_cut(80, 0., 80.)", "1 " + ang_cut + pos_cut)
startZhist = RT.gDirectory.Get("startZ_ang_pos_cut")

outtree.Draw( "startX>>startX_ang_pos_cut(40, -100., 100.)", "1 " + ang_cut + pos_cut)
startXhist = RT.gDirectory.Get("startX_ang_pos_cut")

outtree.Draw( "startY>>startY_ang_pos_cut(40, 380., 500.)", "1 " + ang_cut + pos_cut)
startYhist = RT.gDirectory.Get("startY_ang_pos_cut")

outtree.Draw( "chi2>>chi2_ang_pos_cut(100, 0., 400.)", "1 " + ang_cut + pos_cut)
chi2hist = RT.gDirectory.Get("chi2_ang_pos_cut")

outtree.Draw( "cnn>>cnn_ang_pos_cut(100, 0., 1.)", "1 " + ang_cut + pos_cut)
cnnhist = RT.gDirectory.Get("cnn_ang_pos_cut")


first_cut_dir.cd()
set_style(lenhist, "Track Length (cm)", "")
markers(lenhist)
lenhist.Write()

set_style(startXhist, "Track Start X (cm)", "")
markers(startXhist)
startXhist.Write()

set_style(startYhist, "Track Start Y (cm)", "")
markers(startYhist)
startYhist.Write()

set_style(startZhist, "Track Start Z (cm)", "")
markers(startZhist)
startZhist.Write()

set_style(endZhist, "Track End Z (cm)", "")
markers(endZhist)
endZhist.Write()

set_style(chi2hist, "#chi^{2}", "")
markers(chi2hist)
chi2hist.Write()

set_style(cnnhist, "cnn", "")
markers(cnnhist)
cnnhist.Write()




###################################################


## Now with length cuts ###
second_cut_dir = fout.mkdir( "second_cut_dir", "Cuts include start position and angular cuts and track length cut")
second_cut_dir.cd()


endZ_cut = " && endZ < 226. "

outtree.Draw( "length>>len_ang_pos_endZ_cut(40,0.,500.)", "1 " + ang_cut + pos_cut + endZ_cut)
lenhist = RT.gDirectory.Get("len_ang_pos_endZ_cut")

outtree.Draw( "endZ>>endZ_ang_pos_endZ_cut(40,0.,500.)", "1 " + ang_cut + pos_cut + endZ_cut)
endZhist = RT.gDirectory.Get("endZ_ang_pos_endZ_cut")

outtree.Draw( "startZ>>startZ_ang_pos_endZ_cut(80, 0., 80.)", "1 " + ang_cut + pos_cut + endZ_cut)
startZhist = RT.gDirectory.Get("startZ_ang_pos_endZ_cut")

outtree.Draw( "startX>>startX_ang_pos_endZ_cut(40, -100., 100.)", "1 " + ang_cut + pos_cut + endZ_cut)
startXhist = RT.gDirectory.Get("startX_ang_pos_endZ_cut")

outtree.Draw( "startY>>startY_ang_pos_endZ_cut(40, 380., 500.)", "1 " + ang_cut + pos_cut + endZ_cut)
startYhist = RT.gDirectory.Get("startY_ang_pos_endZ_cut")

outtree.Draw( "chi2>>chi2_ang_pos_endZ_cut(100, 0., 400.)", "1 " + ang_cut + pos_cut + endZ_cut)
chi2hist = RT.gDirectory.Get("chi2_ang_pos_endZ_cut")

outtree.Draw( "cnn>>cnn_ang_pos_endZ_cut(100, 0., 1.)", "1 " + ang_cut + pos_cut + endZ_cut)
cnnhist = RT.gDirectory.Get("cnn_ang_pos_endZ_cut")


second_cut_dir.cd()

set_style(lenhist, "Track Length (cm)", "")
markers(lenhist)
lenhist.Write()

set_style(startXhist, "Track Start X (cm)", "")
markers(startXhist)
startXhist.Write()

set_style(startYhist, "Track Start Y (cm)", "")
markers(startYhist)
startYhist.Write()

set_style(startZhist, "Track Start Z (cm)", "")
markers(startZhist)
startZhist.Write()

set_style(endZhist, "Track End Z (cm)", "")
markers(endZhist)
endZhist.Write()

set_style(chi2hist, "#chi^{2}", "")
markers(chi2hist)
chi2hist.Write()

set_style(cnnhist, "cnn", "")
markers(cnnhist)
cnnhist.Write()

###########################

## Now with chi2 cuts ###
third_cut_dir = fout.mkdir( "third_cut_dir", "Cuts include start position and angular cuts and track length cut and chi2 cut")
third_cut_dir.cd()


endZ_cut = " && endZ < 226. "
chi2_cut = " && chi2 > 140. "


outtree.Draw( "length>>len_ang_pos_endZ_chi2_cut(40,0.,500.)", "1 " + ang_cut + pos_cut + endZ_cut + chi2_cut)
lenhist = RT.gDirectory.Get("len_ang_pos_endZ_chi2_cut")

outtree.Draw( "endZ>>endZ_ang_pos_endZ_chi2_cut(40,0.,500.)", "1 " + ang_cut + pos_cut + endZ_cut + chi2_cut)
endZhist = RT.gDirectory.Get("endZ_ang_pos_endZ_chi2_cut")

outtree.Draw( "startZ>>startZ_ang_pos_endZ_chi2_cut(80, 0., 80.)", "1 " + ang_cut + pos_cut + endZ_cut + chi2_cut)
startZhist = RT.gDirectory.Get("startZ_ang_pos_endZ_chi2_cut")

outtree.Draw( "startX>>startX_ang_pos_endZ_chi2_cut(40, -100., 100.)", "1 " + ang_cut + pos_cut + endZ_cut + chi2_cut)
startXhist = RT.gDirectory.Get("startX_ang_pos_endZ_chi2_cut")

outtree.Draw( "startY>>startY_ang_pos_endZ_chi2_cut(40, 380., 500.)", "1 " + ang_cut + pos_cut + endZ_cut + chi2_cut)
startYhist = RT.gDirectory.Get("startY_ang_pos_endZ_chi2_cut")

outtree.Draw( "chi2>>chi2_ang_pos_endZ_chi2_cut(100, 0., 400.)", "1 " + ang_cut + pos_cut + endZ_cut + chi2_cut)
chi2hist = RT.gDirectory.Get("chi2_ang_pos_endZ_chi2_cut")

outtree.Draw( "cnn>>cnn_ang_pos_endZ_chi2_cut(100, 0., 1.)", "1 " + ang_cut + pos_cut + endZ_cut + chi2_cut)
cnnhist = RT.gDirectory.Get("cnn_ang_pos_endZ_chi2_cut")


third_cut_dir.cd()

set_style(lenhist, "Track Length (cm)", "")
markers(lenhist)
lenhist.Write()

set_style(startXhist, "Track Start X (cm)", "")
markers(startXhist)
startXhist.Write()

set_style(startYhist, "Track Start Y (cm)", "")
markers(startYhist)
startYhist.Write()

set_style(startZhist, "Track Start Z (cm)", "")
markers(startZhist)
startZhist.Write()

set_style(endZhist, "Track End Z (cm)", "")
markers(endZhist)
endZhist.Write()

set_style(chi2hist, "#chi^{2}", "")
markers(chi2hist)
chi2hist.Write()

set_style(cnnhist, "cnn", "")
markers(cnnhist)
cnnhist.Write()


