import ROOT as RT
from array import array
#from vertex_type import vertex_type as vt
from defcuts import ang_pos_test_cut, data_ang_pos_test_cut
from math import sqrt
#from chi2 import chi2 as do_chi2

from argparse import ArgumentParser as ap

parser = ap()

parser.add_argument( "-i", type=str, help='Input file' )
parser.add_argument( "-o", type=str, help='Output file' )
parser.add_argument( "-v", type=int, help='Do Vertex type?', default=0 )
parser.add_argument( "-c", type=float, help='CNN score for cut', default=.3 )
parser.add_argument( "--ctype", type=str, help='Which CNN to use', default='full')
parser.add_argument( "--dR", type=float, help='dR for cut', default=-1. )
parser.add_argument( "--dtype", type=int, help='Which cuts to use', default=1)
parser.add_argument( "--chi2", type=float, help='Which beam chi2 cut', default=-1.)


args = parser.parse_args()


f = RT.TFile( args.i )
tree = f.Get("pionana/beamana")

template_file = RT.TFile("dEdxrestemplates.root")
profile = template_file.Get("dedx_range_pro")

fout = RT.TFile( args.o, "RECREATE" )
outtree = RT.TTree("tree","")

dR = array("d", [0.])
shower_dR = array("d", [0.])
shower_len_per_dR = array("d", [0.])
shower_Energy = array("d", [0.])
shower_Energy_per_hit = array("d", [0.])
shower_cos = array("d", [0.])
vertex_slice = array("i", [0])
daughter_slice = array("i", [0])
is_cosmic = array("i", [0])
pdg = array("i", [0])
par_pdg = array("i", [0])
matched = array("i", [0])

chi2 = array("d", [0])
new_chi2 = array("d", [0])
cnn  = array("d", [0])
cnn_collection  = array("d", [0])
cnn_michel  = array("d", [0])
pi0_gamma = array("i", [0])
true_daughter = array("i", [0])
true_grand_daughter = array("i", [0])
true_great_grand_daughter = array("i", [0])
michel_from_pi_daughter = array("i", [0])
true_gd_ParPDG = array("i", [0])
is_self = array("i", [0])
event = array("i", [0])
subrun = array("i", [0])
run = array("i", [0])
beamTrackID = array("i", [0])
daughterTrackID = array("i", [0])
daughterPFPID = array("i", [0])

nHits = array("i", [0])

is_track = array("i", [0])
is_shower = array("i", [0])
deltaZ = array("d", [0])
endZ = array("d", [0])
startZ = array("d", [0])
shower_deltaZ = array("d", [0])
vertex_type = array("i", [0])

outtree.Branch("dR", dR, "dR/D")
outtree.Branch("vertex_slice", vertex_slice, "vertex_slice/I")
outtree.Branch("daughter_slice", daughter_slice, "daughter_slice/I")
outtree.Branch("is_cosmic", is_cosmic, "is_cosmic/I")
outtree.Branch("pdg", pdg, "pdg/I")
outtree.Branch("matched", matched, "matched/I")
outtree.Branch("par_pdg", par_pdg, "par_pdg/I")
outtree.Branch("chi2", chi2, "chi2/D")
outtree.Branch("new_chi2", new_chi2, "new_chi2/D")
outtree.Branch("cnn", cnn, "cnn/D")
outtree.Branch("cnn_collection", cnn_collection, "cnn_collection/D")
outtree.Branch("cnn_michel", cnn_michel, "cnn_michel/D")
outtree.Branch("pi0_gamma", pi0_gamma, "pi0_gamma/I")
outtree.Branch("true_daughter", true_daughter, "true_daughter/I")
outtree.Branch("true_grand_daughter", true_grand_daughter, "true_grand_daughter/I")
outtree.Branch("true_gd_ParPDG", true_gd_ParPDG, "true_gd_ParPDG/I")
outtree.Branch("true_great_grand_daughter", true_great_grand_daughter, "true_great_grand_daughter/I")
outtree.Branch("michel_from_pi_daughter", michel_from_pi_daughter, "michel_from_pi_daughter/I")
outtree.Branch("is_self", is_self, "is_self/I")
outtree.Branch("event", event, "event/I")
outtree.Branch("subrun", subrun, "subrun/I")
outtree.Branch("run", run, "run/I")
outtree.Branch("beamTrackID", beamTrackID, "beamTrackID/I")
outtree.Branch("daughterTrackID", daughterTrackID, "daughterTrackID/I")
outtree.Branch("daughterPFPID", daughterPFPID, "daughterPFPID/I")
outtree.Branch("shower_dR", shower_dR, "shower_dR/D")
outtree.Branch("shower_len_per_dR", shower_len_per_dR, "shower_len_per_dR/D")
outtree.Branch("shower_Energy", shower_Energy, "shower_Energy/D")
outtree.Branch("shower_Energy_per_hit", shower_Energy_per_hit, "shower_Energy_per_hit/D")
outtree.Branch("shower_cos", shower_cos, "shower_cos/D")
outtree.Branch("is_track", is_track, "is_track/I")
outtree.Branch("is_shower", is_shower, "is_shower/I")
outtree.Branch("deltaZ", deltaZ, "deltaZ/D")
outtree.Branch("endZ", endZ, "endZ/D")
outtree.Branch("startZ", startZ, "startZ/D")
outtree.Branch("shower_deltaZ", shower_deltaZ, "shower_deltaZ/D")

outtree.Branch("vertex_type", vertex_type, "vertex_type/I")

outtree.Branch("nHits", nHits, "nHits/I")

RT.gROOT.SetBatch(1)

#print "Entries:", tree.GetEntries()

isMC = 0
for e in tree:

  if e.MC:
    isMC = 1
    if not ( e.true_beam_PDG in [211, -13, 13] and e.reco_beam_type == 13 ): continue
    #if not ( e.reco_beam_true_byHits_matched and e.true_beam_PDG in [211, -13, 13] and e.reco_beam_type == 13 ): continue
    #if not ( ang_pos_test_cut( e, xlow=-3, xhigh=0., ylow=-1., yhigh=2., zlow=28., zhigh=32. ) ): continue
    if not ( ang_pos_test_cut( e) ): continue
  else:
    if not data_ang_pos_test_cut(e): continue 
    if not ( e.reco_beam_type == 13 and (211 in [i for i in e.data_BI_PDG_candidates])): continue
    #if not (e.reco_beam_Chi2_proton / e.reco_beam_Chi2_ndof > 50. and e.reco_beam_Chi2_proton / e.reco_beam_Chi2_ndof < 500 ): continue

  if( e.reco_beam_endZ > 226. ): continue
  if( args.chi2 > 0. ):
    if( e.reco_beam_Chi2_proton / e.reco_beam_Chi2_ndof < args.chi2): continue

  #if e.MC: vertex_type[0] = vt(e, 5., 3)
  event[0] = e.event
  subrun[0] = e.subrun
  run[0] = e.run
  beamTrackID[0] = e.reco_beam_trackID  

  matched[0] = e.reco_beam_true_byHits_matched

  vertex_slice[0] = e.reco_beam_vertex_slice

  if e.MC: pi0_decay_IDs = [i for i in e.true_beam_Pi0_decay_ID]

  pfp_ids = [i for i in e.reco_daughter_PFP_ID]
  for i in range(0, len(pfp_ids)):

    daughterPFPID[0] = pfp_ids[i]
    daughterTrackID[0] = e.reco_daughter_allTrack_ID[i]

      
    if e.MC: 
      pdg[0] = e.reco_daughter_PFP_true_byHits_PDG[i]
      is_cosmic[0] = int( e.reco_daughter_PFP_true_byHits_origin[i] == 2 )
      is_self[0] = int( e.reco_daughter_PFP_true_byHits_ID[i] == e.true_beam_ID )

    cnn[0] = e.reco_daughter_PFP_trackScore[i] 
    cnn_collection[0] = e.reco_daughter_PFP_trackScore_collection[i] 
    cnn_michel[0] = e.reco_daughter_PFP_michelScore[i]  

    dR[0] = -1.
    shower_dR[0] = -1.
    shower_len_per_dR[0] = -1.
    shower_Energy[0] = -1.
    shower_Energy_per_hit[0] = -1.
    shower_cos[0] = -1.

    nHits[0] = e.reco_daughter_PFP_nHits[i]


    bX = e.reco_beam_endX
    bY = e.reco_beam_endY
    bZ = e.reco_beam_endZ


    if( e.reco_daughter_allTrack_ID[i] != -1 and e.reco_daughter_allTrack_Chi2_ndof[i] > 0):
      #dR[0] = e.reco_daughter_allTrack_to_vertex[i]
      #daughter_slice[0] = e.reco_daughter_allTrack_slice[i]
      chi2[0] = e.reco_daughter_allTrack_Chi2_proton[i] / e.reco_daughter_allTrack_Chi2_ndof[i]
      #this_chi2 = do_chi2(e,i, profile)
      #chi2[0] = this_chi2[0] / this_chi2[1] 
      is_track[0] = 1

      start_dX = e.reco_daughter_allTrack_startX[i]
      start_dY = e.reco_daughter_allTrack_startY[i]
      start_dZ = e.reco_daughter_allTrack_startZ[i]
      end_dX = e.reco_daughter_allTrack_endX[i]
      end_dY = e.reco_daughter_allTrack_endY[i]
      end_dZ = e.reco_daughter_allTrack_endZ[i]

      startZ[0] = start_dZ
      endZ[0] = end_dZ


      start_dR = sqrt( (start_dZ - bZ)**2 + (start_dX - bX)**2 + (start_dY - bY)**2 )
      end_dR = sqrt( (end_dZ - bZ)**2 + (end_dX - bX)**2 + (end_dY - bY)**2 )

      if abs(start_dR) < abs(end_dR): 
        deltaZ[0] = start_dZ - bZ
        dR[0] = start_dR
      else: 
        deltaZ[0] =  end_dZ - bZ 
        dR[0] = end_dR
 
    else: 
      is_track[0] = 0
      startZ[0] = -999.
      endZ[0] = -999.
      chi2[0] = -999.
      dR[0] = -999.
      deltaZ[0] = -999.

    if not e.reco_daughter_allShower_ID[i] == -1:
      dX = e.reco_daughter_allShower_startX[i]
      dY = e.reco_daughter_allShower_startY[i]
      dZ = e.reco_daughter_allShower_startZ[i]

      shower_dR[0] = sqrt( (dX - bX)**2 + (dY - bY)**2 + (dZ - bZ)**2 )
      shower_len_per_dR[0] = shower_dR[0] / e.reco_daughter_allShower_len[i]
      shower_deltaZ[0] = dZ - bZ 

      shower_Energy[0] = e.reco_daughter_allShower_energy[i] 
      shower_Energy_per_hit[0] = e.reco_daughter_allShower_energy[i]/nHits[0]
      shower_cos[0] = (e.reco_daughter_allShower_dirX[i]*e.reco_beam_trackEndDirX +
                       e.reco_daughter_allShower_dirY[i]*e.reco_beam_trackEndDirY +
                       e.reco_daughter_allShower_dirZ[i]*e.reco_beam_trackEndDirZ)
      is_shower[0] = 1
    else: 
      is_shower[0] = 0
       

    pi0_gamma[0] = 0
    if e.MC: 
      #check if this is a pi0 gamma:
      if( pdg[0] == 22 ):
        if 111 in [j for j in e.true_beam_daughter_PDG]:
          if e.reco_daughter_PFP_true_byHits_ID[i] in pi0_decay_IDs:
            pi0_gamma[0] = 1

      true_daughter[0] = int( e.reco_daughter_PFP_true_byHits_parID[i] == e.true_beam_ID )
      true_grand_daughter[0] = int( e.reco_daughter_PFP_true_byHits_ID[i] in [j for j in e.true_beam_grand_daughter_ID] )
      true_great_grand_daughter[0] = int( e.reco_daughter_PFP_true_byHits_parID[i] in [j for j in e.true_beam_grand_daughter_ID] )

      michel_from_pi_daughter[0] = 0
      if true_great_grand_daughter[0] and abs(pdg[0]) == 11:   
        a = 0
        for gdID,gdProc in zip( [j for j in e.true_beam_grand_daughter_ID], [str(j) for j in e.true_beam_grand_daughter_endProcess] ):
          if e.reco_daughter_PFP_true_byHits_parID[i] == gdID and ( gdProc == "Decay" or "FastScint" in gdProc ):
            b = 0
            for dID,dProc in zip( [j for j in e.true_beam_daughter_ID], [str(j) for j in e.true_beam_daughter_endProcess] ):
              if e.true_beam_grand_daughter_parID[a] == dID and ( "FastScint" in dProc or "Decay" in dProc ) and abs(e.true_beam_daughter_PDG[b]) == 211:
                michel_from_pi_daughter[0] = 1   
                break
              b += 1
            break
          a += 1

      if true_grand_daughter[0]:
        true_gd_IDs = [j for j in e.true_beam_grand_daughter_ID]
        true_gd_ParIDs = [j for j in e.true_beam_grand_daughter_parID]

        for gd_ID, gd_ParID in zip(true_gd_IDs, true_gd_ParIDs):
          if e.reco_daughter_PFP_true_byHits_ID[i] == gd_ID:
            for d_ID, d_PDG in zip([j for j in e.true_beam_daughter_ID], [j for j in e.true_beam_daughter_PDG]):
              if gd_ParID == d_ID:
                true_gd_ParPDG[0] = d_PDG

    outtree.Fill()

if isMC: 
  if args.dtype:
    cats = [
      "Cosmic", "Nuc", "Self",
      "Gamma", "pi0Gamma",
      "Other", "GDaughter", "GGDaughter", "Michel", "Pi", "Proton", "Mu"
    ]
    
    colors = {
      "Cosmic": RT.kRed+2,
      "pi0Gamma": RT.kSpring-8,
      "Pi": RT.kTeal,
      "Mu": RT.kRed,
      "Proton": RT.kViolet-3,
      "Gamma": RT.kOrange+1,
      "Nuc": RT.kOrange+10,
      "Other": RT.kBlack,
      "GDaughter": RT.kMagenta,
      "GGDaughter": RT.kOrange,
      "Michel": RT.kPink,
      "Self": RT.kBlue
    }

    leg_names = {
      "Cosmic"    : "Cosmic",
      "pi0Gamma"  : "#pi^{0} Shower",
      "Pi"        : "#pi^{#pm}", 
      "Mu"        : "#mu",
      "Proton"    : "p",
      "Gamma"     : "#gamma",
      "Nuc"       : "Nuclear fragment",
      "Other"     : "Other",
      "GDaughter" : "Grand Daughter",
      "GGDaughter": "Great-Grand Daughter",
      "Michel"    : "Michel electron",
      "Self"      : "Self"
    }
    
    cuts = {
      "Cosmic": "is_cosmic",
      "pi0Gamma": "pi0_gamma",
      "Pi": "true_daughter && abs(pdg) == 211",
      "Mu": "true_daughter && abs(pdg) == 13",
      "Proton": "true_daughter && abs(pdg) == 2212",
      "Gamma":  "true_daughter && abs(pdg) == 22",
      "Nuc":  "true_daughter && abs(pdg) > 2212",
      "GDaughter": "true_grand_daughter && !pi0_gamma",
      "GGDaughter": "true_great_grand_daughter && !michel_from_pi_daughter",
      "Michel": "michel_from_pi_daughter || (true_daughter && abs(pdg) == 11)",
      "Other": "!(true_daughter || true_grand_daughter || true_great_grand_daughter || is_cosmic || pi0_gamma || is_self)",
      "Self": "is_self"
    }
  else:
    cats = [
      "Mu", "Nuc",
      "Gamma", 
      "Other", "e", "Pi", "Proton"
    ]
    
    colors = {
      "Gamma": RT.kSpring-8,
      "Pi": RT.kTeal,
      "Mu": RT.kRed,
      "Proton": RT.kViolet-3,
      "Nuc": RT.kOrange,
      "Other": RT.kBlack,
      "e": RT.kBlue
    }
    
    
    cuts = {
      #"pi0Gamma": "pi0_gamma",
      "Pi": "abs(pdg) == 211",
      "Mu": "abs(pdg) == 13",
      "Proton": "abs(pdg) == 2212",
      "Gamma":  "abs(pdg) == 22",
      "Nuc":  "abs(pdg) > 2212",
      "e":   "abs(pdg) == 11", 
      "Other": "!( abs(pdg) == 211 || abs(pdg) == 13 || abs(pdg) == 2212 || abs(pdg) == 22 || abs(pdg) > 2212 || abs(pdg) == 11)"
    }

    leg_names = {
      "Pi": "#pi^{#pm}",
      "Mu": "#mu",
      "Proton": "p",
      "Gamma": "#gamma",
      "Nuc": "Nuclear fragment",
      "e": "e",
      "Other": "Other"
    }
else:
  cats = ["Data"]
  colors = { "Data": 1}
  cuts = {"Data":"1"}


#Looking at dR -- distance between Closer(end/start) of daughter and closest point on beam track
dR_hists = dict()
shower_dR_hists = dict()
shower_len_per_dR_hists = dict()
shower_Energy_hists = dict()
shower_Energy_per_hit_hists = dict()
shower_cos_hists = dict()
deltaZ_hists = dict()
shower_deltaZ_hists = dict()
slice_hists = dict()
chi2_hists = dict()
cnn_hists = dict()
cnn_collection_hists = dict()
cnn_michel_hists = dict()

track_nHits_hists = dict()
shower_nHits_hists = dict()

leg = RT.TLegend(.6,.6,.85,.85)


#print outtree.GetEntries()

do_vertex_type = False
if args.v:
  "doing vertex_type"
  do_vertex_type = True

if isMC and do_vertex_type: nVT = 5
else: nVT = 1

cnn_cut = str(args.c)

if(args.c == "full"):
  cnn_str = "cnn"
else:
  cnn_str = "cnn_collection"

for i in range(0, nVT):
  for cat in cats: 

    theCut = cuts[cat] 
    if do_vertex_type:
      #print "doing cut"
      theCut = theCut + " && vertex_type == " + str(i)


    #print cuts[cat] + " && vertex_type == " + str(i) 
    outtree.Draw("dR>>h" + cat + str(i) + "(175,0,350)", theCut + " && " + cnn_str + " > " + cnn_cut + " && is_track")
    dR_hists[cat + str(i)] = RT.gDirectory.Get("h"+cat + str(i))
    #print dR_hists[cat + str(i)].Integral()
    dR_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: dR_hists[cat + str(i)].SetFillColor(colors[cat])
    #print cat, dR_hists[cat + str(i)].Integral()

    outtree.Draw("shower_dR>>h_shower_" + cat + str(i) + "(350,0,350)", theCut + " && " + cnn_str + " < " + cnn_cut + " && is_shower")
    shower_dR_hists[cat + str(i)] = RT.gDirectory.Get("h_shower_"+cat + str(i))
    shower_dR_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: shower_dR_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw("shower_len_per_dR>>h_shower_lendr" + cat + str(i) + "(500,0,20)", theCut + " && " + cnn_str + " < " + cnn_cut + " && is_shower")
    shower_len_per_dR_hists[cat + str(i)] = RT.gDirectory.Get("h_shower_lendr"+cat + str(i))
    shower_len_per_dR_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: shower_len_per_dR_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw("shower_Energy>>h_shower_Energy" + cat + str(i) + "(50,0,500)", theCut + " && " + cnn_str + " < " + cnn_cut + " && is_shower")
    shower_Energy_hists[cat + str(i)] = RT.gDirectory.Get("h_shower_Energy" + cat + str(i))
    shower_Energy_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: shower_Energy_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw("shower_Energy_per_hit>>h_shower_Energy_per_hit" + cat + str(i) + "(50,0,5)", theCut + " && " + cnn_str + " < " + cnn_cut + " && is_shower")
    shower_Energy_per_hit_hists[cat + str(i)] = RT.gDirectory.Get("h_shower_Energy_per_hit" + cat + str(i))
    shower_Energy_per_hit_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: shower_Energy_per_hit_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw("shower_cos>>h_shower_cos" + cat + str(i) + "(50, -1., 1.)", theCut + " && " + cnn_str + " < " + cnn_cut + " && is_shower")
    shower_cos_hists[cat + str(i)] = RT.gDirectory.Get("h_shower_cos" + cat + str(i))
    shower_cos_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: shower_cos_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw("deltaZ>>h_deltaZ_" + cat + str(i) + "(175,-175,175)", theCut + " && " + cnn_str + " > " + cnn_cut + " && is_track")
    deltaZ_hists[cat + str(i)] = RT.gDirectory.Get("h_deltaZ_"+cat + str(i))
    deltaZ_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: deltaZ_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw("shower_deltaZ>>h_shower_deltaZ_" + cat + str(i) + "(175,-175,175)", theCut + " && " + cnn_str + " < " + cnn_cut + " && is_shower")
    shower_deltaZ_hists[cat + str(i)] = RT.gDirectory.Get("h_shower_deltaZ_"+cat + str(i))
    shower_deltaZ_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: shower_deltaZ_hists[cat + str(i)].SetFillColor(colors[cat])

  
    outtree.Draw( "vertex_slice - daughter_slice>>hSlice" + cat + str(i) + "(600,0,600)", theCut )
    slice_hists[cat + str(i)] = RT.gDirectory.Get("hSlice"+cat + str(i))
    slice_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: slice_hists[cat + str(i)].SetFillColor(colors[cat])
  
    chi2_cut = cuts[cat]  + " && vertex_type == " + str(i)

    outtree.Draw( "chi2>>hChi2" + cat + str(i) + "(100,0,500)", theCut + " && " + cnn_str + " > " + cnn_cut + " && is_track" )
    chi2_hists[cat + str(i)] = RT.gDirectory.Get("hChi2"+cat + str(i))
    chi2_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: chi2_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw( "cnn>>hCNN" + cat + str(i) + "(100,0,1)", theCut)
    cnn_hists[cat + str(i)] = RT.gDirectory.Get("hCNN"+cat + str(i))
    cnn_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: cnn_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw( "cnn_collection>>hCNN_collection" + cat + str(i) + "(100,0,1)", theCut)
    cnn_collection_hists[cat + str(i)] = RT.gDirectory.Get("hCNN_collection"+cat + str(i))
    cnn_collection_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: cnn_collection_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw( "cnn_michel>>hCNN_michel" + cat + str(i) + "(100,0,1)", theCut  )
    cnn_michel_hists[cat + str(i)] = RT.gDirectory.Get("hCNN_michel"+cat + str(i))
    cnn_michel_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: cnn_michel_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw("nHits>>h_shower_nHits" + cat + str(i) + "(500,0,2000)", theCut + " && " + cnn_str + " < " + cnn_cut + " && is_shower")
    shower_nHits_hists[cat + str(i)] = RT.gDirectory.Get("h_shower_nHits"+cat + str(i))
    shower_nHits_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: shower_nHits_hists[cat + str(i)].SetFillColor(colors[cat])

    outtree.Draw("nHits>>h_track_nHits" + cat + str(i) + "(100,0,2000)", theCut + " && " + cnn_str + " > " + cnn_cut + " && is_track")
    track_nHits_hists[cat + str(i)] = RT.gDirectory.Get("h_track_nHits"+cat + str(i))
    track_nHits_hists[cat + str(i)].SetLineColor(colors[cat])
    if isMC: track_nHits_hists[cat + str(i)].SetFillColor(colors[cat])

  
  dR_stack = RT.THStack("dR_stack" + str(i), "")
  shower_dR_stack = RT.THStack("shower_dR_stack" + str(i), "")
  shower_len_per_dR_stack = RT.THStack("shower_len_per_dR_stack" + str(i), "")
  shower_Energy_per_hit_stack = RT.THStack("shower_Energy_per_hit_stack" + str(i), "")
  shower_Energy_stack = RT.THStack("shower_Energy_stack" + str(i), "")
  shower_cos_stack = RT.THStack("shower_cos_stack" + str(i), "")
  shower_deltaZ_stack = RT.THStack("shower_deltaZ_stack" + str(i), "")
  deltaZ_stack = RT.THStack("deltaZ_stack" + str(i), "")
  slice_stack = RT.THStack("slice_stack" + str(i), "")
  chi2_stack = RT.THStack("chi2_stack" + str(i), "")
  cnn_stack = RT.THStack("cnn_stack" + str(i), "") 
  cnn_collection_stack = RT.THStack("cnn_collection_stack" + str(i), "") 
  cnn_michel_stack = RT.THStack("cnn_michel_stack" + str(i), "") 
  shower_nHits_stack = RT.THStack("shower_nHits_stack" + str(i), "")
  track_nHits_stack  = RT.THStack("track_nHits_stack" + str(i), "")

  for cat in cats:
    dR_stack.Add(dR_hists[cat + str(i)])
    if i == 0 and isMC: 
      leg.AddEntry(dR_hists[cat + str(i)], leg_names[cat], "lf")
  
    slice_stack.Add(slice_hists[cat + str(i)])
    chi2_stack.Add(chi2_hists[cat + str(i)])
    cnn_stack.Add(cnn_hists[cat + str(i)])
    cnn_collection_stack.Add(cnn_collection_hists[cat + str(i)])
    cnn_michel_stack.Add(cnn_michel_hists[cat + str(i)])
    shower_dR_stack.Add(shower_dR_hists[cat + str(i)])
    shower_len_per_dR_stack.Add(shower_len_per_dR_hists[cat + str(i)])
    shower_Energy_per_hit_stack.Add(shower_Energy_per_hit_hists[cat + str(i)])
    shower_Energy_stack.Add(shower_Energy_hists[cat + str(i)])
    shower_cos_stack.Add(shower_cos_hists[cat + str(i)])
    deltaZ_stack.Add(deltaZ_hists[cat + str(i)])
    shower_deltaZ_stack.Add(shower_deltaZ_hists[cat + str(i)])
    shower_nHits_stack.Add(shower_nHits_hists[cat + str(i)])
    track_nHits_stack.Add(track_nHits_hists[cat + str(i)])
  
  dR_stack.SetTitle(";dR (cm);")
  dR_stack.Write("dR_stack" + str(i))
  shower_dR_stack.Write("shower_dR_stack" + str(i))
  shower_len_per_dR_stack.Write("shower_len_per_dR_stack" + str(i))
  shower_Energy_stack.Write("shower_Energy_stack" + str(i))
  shower_Energy_per_hit_stack.Write("shower_Energy_per_hit_stack" + str(i))
  shower_cos_stack.Write("shower_cos_stack" + str(i))
  deltaZ_stack.Write("deltaZ_stack" + str(i))
  shower_deltaZ_stack.Write("shower_deltaZ_stack" + str(i))
  shower_nHits_stack.Write("shower_nHits_stack" + str(i))
  track_nHits_stack.Write("track_nHits_stack" + str(i))
  slice_stack.Write("slice_stack" + str(i))
  chi2_stack.Write("chi2_stack" + str(i))
  cnn_stack.Write("cnn_stack" + str(i))
  cnn_collection_stack.Write("cnn_collection_stack" + str(i))
  cnn_michel_stack.Write("cnn_michel_stack" + str(i))
leg.Write("leg")

fout.cd()
outtree.Write()
fout.Close()
