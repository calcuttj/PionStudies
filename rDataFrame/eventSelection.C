#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "lambda.h"
#include "eventSelection.h"
#include <ROOT/RDataFrame.hxx>

#include "TGraphAsymmErrors.h"

#include "backgrounds.h"
#include "selection_defs.h"

#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <chrono>

//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;

using namespace std::chrono; 


//***********************
//Main Function
int eventSelection_primMu(const string mcFile, const string dataFile,
                           bool doBatch = false) {

  //This prevents the canvas from being draw at the end
  //Useful for when on the gpvms 
  gROOT->SetBatch(doBatch);
  //ROOT::EnableImplicitMT(4);

  gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
  
  ROOT::RDataFrame frame(inputTree, mcFile);
  ROOT::RDataFrame data_frame(inputTree, dataFile);
  //TFile *output = new TFile ("output_eventSelection.root", "RECREATE");
  //THStack *stack_cutFlow = new THStack("cutFlow", "Cut Flow MC and Data");

  //This creates columns in the frame
  //which we'll use to filter through later
  auto mc_all = frame
    .Define("primaryMuon", primaryMuon, 
            {"reco_beam_true_byHits_PDG", "reco_beam_true_byHits_process",
             "reco_beam_true_byHits_matched", "reco_beam_true_byHits_origin"})

    .Define("isCosmic", isCosmic, {"reco_beam_true_byHits_origin"})

    .Define("isDecay", isDecay, {"reco_beam_true_byHits_process",
                                 "reco_beam_true_byHits_origin"})

    .Define("isExtraBeam", isExtraBeam,
            {"reco_beam_true_byHits_process", "reco_beam_true_byHits_matched",
             "reco_beam_true_byHits_origin"})

    .Define("upstreamInt", upstreamInt,
            {"reco_beam_true_byHits_process", "reco_beam_true_byHits_origin"})

    //tag if there is a true Pion with high momentum in event, threshold
    //defined in eventSelection.h
    .Define("true_daughter_pion_momentumHigh", tagDaughterPionMomentumHigh, 
          {"true_beam_daughter_PDG", "true_beam_daughter_startP",
           "true_daughter_nPiPlus", "true_daughter_nPiMinus"})

    
    .Define("true_primPionInel", tagPrimPionInel,
            {"true_beam_PDG", "true_beam_endProcess"})

    .Define("true_primPionInel_withElastic", tagPrimPionInel_withElastic,
            {"true_beam_PDG", "true_beam_endProcess",
             "true_beam_nElasticScatters"})

    .Define("true_combinedSignal", tagAbsChEx,
            {"true_primPionInel", "true_daughter_pion_momentumHigh"})

    .Define("true_chexSignal", tagChEx, {"true_combinedSignal",
                                         "true_daughter_nPi0"})
    .Define("true_absSignal", tagAbs, {"true_combinedSignal",
                                       "true_daughter_nPi0"})
    .Define("true_nPi0Signal", tagNpi0, {"true_combinedSignal",
                                         "true_daughter_nPi0"})
    .Define("true_backGround", tagBackGround, {"true_primPionInel",
                                               "true_combinedSignal"})
    .Define("true_pion_daughter",
            "true_daughter_nPiPlus + true_daughter_nPiMinus")

    //necessary Reco Values
    //FRANCESCA: not using distance anymore.. kick out?
    /*.Define("daughter_distance3D", compute_distanceVertex,
            {"reco_beam_endX", "reco_beam_endY", "reco_beam_endZ", 
             "reco_daughter_allTrack_startX", "reco_daughter_allTrack_startY",
             "reco_daughter_allTrack_startZ", "reco_daughter_allTrack_endX",
             "reco_daughter_allTrack_endY", "reco_daughter_allTrack_endZ"})

    //Jake: should we replace the second set of shower locatons with something 
    //      like 9999. for clarity? 
    .Define("daughter_distance3D_shower", compute_distanceVertex, 
            {"reco_beam_endX", "reco_beam_endY", "reco_beam_endZ",
             "reco_daughter_allShower_startX", "reco_daughter_allShower_startY",
             "reco_daughter_allShower_startZ", "reco_daughter_allShower_startX",
             "reco_daughter_allShower_startY", "reco_daughter_allShower_startZ"})*/

    
    //dQdX values Libo truncated mean Libo keeping dQdX median +- 1sigma 
    //and then calculate mean. i.e. throwing away the upper and lower 16% wrt to median (keeping 68.2 percent the values)
    .Define("libo_low", [](){return 0.16;})
    .Define("libo_high",[](){return 0.16;})
    .Define("reco_daughter_allTrack_truncLibo_dQdX", truncatedMean, 
          {"libo_low", "libo_high", "reco_daughter_allTrack_dQdX_SCE"})
    .Define("reco_daughter_allTrack_truncLibo_dEdX",truncatedMean, 
          {"libo_low", "libo_high", "reco_daughter_allTrack_calibrated_dEdX_SCE"} )
    .Define("interaction_topology", interaction_topology,
            {"reco_beam_true_byHits_origin", "true_beam_PDG", "true_beam_ID",
             "true_beam_endZ", "true_beam_endProcess", "true_daughter_nPi0",
             "reco_beam_true_byHits_origin", "reco_beam_true_byHits_ID",
             "true_beam_daughter_PDG", "true_beam_daughter_startP"})
    //Filter for true primary Pion and Beam Muon
    .Filter("true_beam_PDG == 211 || true_beam_PDG == -13");

  std::cout << "MC: Defined and filtered all events" << std::endl;

  //Column definitions in data
  auto data_all = data_frame
    .Define("beamPID", data_beam_PID, {"beam_inst_PDG_candidates"})
    
   /*
    //FRANCESCA: not using distance anymore.. kick out?
    .Define("daughter_distance3D", compute_distanceVertex,
            {"reco_beam_endX", "reco_beam_endY", "reco_beam_endZ",
             "reco_daughter_allTrack_startX", "reco_daughter_allTrack_startY",
             "reco_daughter_allTrack_startZ", "reco_daughter_allTrack_endX",
             "reco_daughter_allTrack_endY",  "reco_daughter_allTrack_endZ"})

    //Jake: should we replace the second set of shower locatons with something 
    //      like 9999. for clarity? 
    .Define("daughter_distance3D_shower", compute_distanceVertex,
            {"reco_beam_endX", "reco_beam_endY", "reco_beam_endZ", 
             "reco_daughter_allShower_startX", "reco_daughter_allShower_startY", 
             "reco_daughter_allShower_startZ", "reco_daughter_allShower_startX", 
             "reco_daughter_allShower_startY", "reco_daughter_allShower_startZ"})
   */

       
    //dQdX values Libo truncated mean Libo keeping dQdX median +- 1sigma 
    //and then calculate mean. i.e. throwing away the upper and lower 26%
    .Define("libo_low", [](){return 0.16;})
    .Define("libo_high",[](){return 0.16;})
    .Define("reco_daughter_allTrack_truncLibo_dQdX", truncatedMean, 
          {"libo_low", "libo_high", "reco_daughter_allTrack_dQdX_SCE"})
    .Define("reco_daughter_allTrack_truncLibo_dEdX",truncatedMean, 
          {"libo_low", "libo_high", "reco_daughter_allTrack_calibrated_dEdX_SCE"} )


    .Filter("beamPID == true"); //Looks for just the events passing the beam line
                                //PID

  std::cout << "Data: Defined and filtered all events" << std::endl;



  //prepare Branches for all the cuts true/false. this allows to do easy filtering
  //all the cuts are in the eventSelection.h file
  //auto mc_all_cutValues = mc_all
  mc_all = mc_all
    .Define("primary_isBeamType", isBeamType, {"reco_beam_type"})

    .Define("passBeamCut", manual_beamPos_mc, 
            {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
             "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ",
             "true_beam_startDirX", "true_beam_startDirY", "true_beam_startDirZ",
             "true_beam_startX", "true_beam_startY", "true_beam_startZ"})
    .Define("passBeamCutBI", beam_cut_MC_BI,
            {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
             "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ",
             "beam_inst_X", "beam_inst_Y", "beam_inst_dirX", "beam_inst_dirY",
             "beam_inst_dirZ", "beam_inst_nMomenta", "beam_inst_nTracks"})
    
    .Define("isPrimaryMuonCandidate", candidate_primaryMuon, {"reco_daughter_PFP_michelScore"})
 
    .Define("primary_ends_inAPA3", endAPA3, {"reco_beam_endZ"})

    .Define("has_noPion_daughter", secondary_noPion,
            {"reco_daughter_PFP_trackScore_collection",
             "reco_daughter_allTrack_ID", 
             "reco_daughter_allTrack_truncLibo_dEdX"})
   
     .Define("has_shower_nHits_distance", has_shower_nHits,
            {"reco_daughter_PFP_trackScore_collection",
             "reco_daughter_PFP_nHits"})

     .Define("selected_incidentPion", "primary_isBeamType && passBeamCut && passBeamCutBI && !isPrimaryMuonCandidate")
     
     .Define("selected_abs", "selected_incidentPion && primary_ends_inAPA3 && has_noPion_daughter && !has_shower_nHits_distance")
     
     .Define("selected_cex", "selected_incidentPion && primary_ends_inAPA3 && has_noPion_daughter && has_shower_nHits_distance");
    
  // DATA
  data_all/*_cutValues*/ = data_all
    .Define("primary_isBeamType", isBeamType, {"reco_beam_type"})
    
    .Define("passBeamQuality", data_BI_quality,
            {"beam_inst_nMomenta", "beam_inst_nTracks"})
    
    .Define("passBeamCut", manual_beamPos_data,
            {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
             "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ",
             "beam_inst_X", "beam_inst_Y", "beam_inst_dirX", "beam_inst_dirY",
             "beam_inst_dirZ", "beam_inst_nMomenta", "beam_inst_nTracks"})

    .Define("isPrimaryMuonCandidate", candidate_primaryMuon, {"reco_daughter_PFP_michelScore"}) 
    
    .Define("primary_ends_inAPA3", endAPA3, {"reco_beam_endZ"})
 
    .Define("has_noPion_daughter", secondary_noPion,
            {"reco_daughter_PFP_trackScore_collection",
             "reco_daughter_allTrack_ID", 
             "reco_daughter_allTrack_truncLibo_dEdX"})
     
     .Define("has_shower_nHits_distance", has_shower_nHits,
            {"reco_daughter_PFP_trackScore_collection",
             "reco_daughter_PFP_nHits"})
     
     .Define("selected_incidentPion", "primary_isBeamType && passBeamCut && passBeamQuality && !isPrimaryMuonCandidate")
     
     .Define("selected_abs", "selected_incidentPion && primary_ends_inAPA3 && has_noPion_daughter && !has_shower_nHits_distance")
     
     .Define("selected_cex", "selected_incidentPion && primary_ends_inAPA3 && has_noPion_daughter && has_shower_nHits_distance");

  //Label within MC files who passed which CUT (this can help to see when what drops out)
  auto mc_output_with_label = mc_all/*_cutValues*/;

  //Create branches to be created in output file
  auto data_output_with_label = data_all/*_cutValues*/;

  //*******************
  //Start Cutting MC
  //******************
  //    Cuts are concatenated
  /* *******Beam Cut******* */

  //Filter out non track-like beam objects 
  std::cout << "Starting cuts" << std::endl;
  auto time0 = high_resolution_clock::now();
  auto mc_snap_all = mc_output_with_label.Snapshot(
      "pionana/beamana", "eventSelection_mc_all.root");

  auto mcCUT_beamType = mc_output_with_label.Filter("primary_isBeamType");
  auto time1 = high_resolution_clock::now();

  std::cout << "MC beam type: " <<
               duration_cast<seconds>(time1 - time0).count() <<
               std::endl;

  auto mc_snap_beam_type = mcCUT_beamType.Snapshot(
      "pionana/beamana", "eventSelection_mc_beamType.root");

  //Beam quality cuts (start position/direction)
  //auto mcCUT_beamCut = mcCUT_beamType.Filter("passBeamCut");
  std::cout << "MC beam cut: " /*<< *N_mcCUT_beamCut */<< std::endl;

  ////Make sure the beam track ends before APA2
  //auto mcCUT_endAPA3 = mcCUT_beamCut.Filter("primary_ends_inAPA3");
  auto mcCUT_endAPA3 = mcCUT_beamType.Filter("primary_ends_inAPA3 && passBeamCut && !(isPrimaryMuonCandidate)");
  auto time2 = high_resolution_clock::now();
  std::cout << "MC APA3 cut: " <<
               duration_cast<seconds>(time2 - time1).count() <<
               std::endl;

  //Make a file with only primary pions in it
  auto mc_snap_primPion = mcCUT_endAPA3.Snapshot(
      "pionana/beamana", "eventSelection_mc_PRIMARYPION.root");

  /* ****** COMBINED SAMPLE ******/

  //no  Pion-like daughter objects 
  auto mcCUT_noPionDaughter = mcCUT_endAPA3.Filter("has_noPion_daughter");

  //Create an output file for the (combined) selected events
  auto mc_COMBINED_Signal = mcCUT_noPionDaughter;
  auto mc_snap_combined = mc_COMBINED_Signal.Snapshot(
      "pionana/beamana", "eventSelection_mc_COMBINED.root");

  //Find pi0 showers from the combined selection
  //and make an output file for them (Cex events) 
  auto mcSIGNAL_cex = mc_COMBINED_Signal.Filter("has_shower_nHits_distance");
  auto mc_snap_cex = mcSIGNAL_cex.Snapshot(
      "pionana/beamana", "eventSelection_mc_CEX.root");

  //Find the selected events without pi0 showers
  //and make an output file (Abs events)
  auto mcSIGNAL_abs = mc_COMBINED_Signal.Filter("!(has_shower_nHits_distance)");
  auto mc_snap_abs = mcSIGNAL_abs.Snapshot(
      "pionana/beamana", "eventSelection_mc_ABS.root");

  //Make a file for events with a tagged pion daughter
  auto mcCUT_PionDaughter = mcCUT_endAPA3.Filter("!has_noPion_daughter");
  auto mc_snap_pion_daughter = mcCUT_PionDaughter.Snapshot(
      "pionana/beamana", "eventSelection_mc_rejected.root");
                        

  //Start Cutting DATA
  //******************
  //    Cuts are concatenated
  /* *******Beam Cut******* */

  //Filter out non track-like beam objects 
  auto data_snap_all = data_output_with_label.Snapshot(
      "pionana/beamana", "eventSelection_data_all.root");

  auto dataCUT_beamQuality = data_output_with_label.Filter("passBeamQuality && primary_isBeamType");
  auto data_snap_beamQuality = dataCUT_beamQuality.Snapshot(
      "pionana/beamana", "eventSelection_data_BeamQuality.root");

  //Ends before APA2
  auto dataCUT_endAPA3 = dataCUT_beamQuality.Filter("primary_ends_inAPA3 && passBeamCut && !(isPrimaryMuonCandidate)");
  
  /* ****** COMBINED SAMPLE ******/

  //no  Pion-like daughter objects 
  auto dataCUT_noPionDaughter = dataCUT_endAPA3.Filter("has_noPion_daughter");

  auto data_COMBINED_Signal = dataCUT_noPionDaughter;
  auto data_snap_combined = data_COMBINED_Signal.Snapshot(
      "pionana/beamana", "eventSelection_data_COMBINED.root");

  //Find pi0 showers from the combined selection
  //and make an output file for them (Cex events) 
  auto dataSIGNAL_cex = data_COMBINED_Signal.Filter("has_shower_nHits_distance");
  auto data_snap_cex = dataSIGNAL_cex.Snapshot(
      "pionana/beamana", "eventSelection_data_CEX.root");

  //Find the selected events without pi0 showers
  //and make an output file (Abs events)
  auto dataSIGNAL_abs = data_COMBINED_Signal.Filter(
      "!(has_shower_nHits_distance)");
  auto data_snap_abs = dataSIGNAL_abs.Snapshot(
      "pionana/beamana", "eventSelection_data_ABS.root");
  return 0;
}

