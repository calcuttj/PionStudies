#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "TMatrixDBase.h"
#include "TArray.h"
#include "../lambda.h"
#include "../betheBloch.h"
#include "eSlice.h"
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include <vector>

//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;

//Applying the eSlice Method

//***********************
//Main Function

//--------------------------------------------------------
//       Generate Eff and Pur of eventSelection up until Pandora Reco stage
//       Step needed to afterwards unfpold
//MC True --------------------> MC Reco ------------------> Selected Interaction, Selected Incident
// 
//correct for Reco Eff         Pandora Reco                 EventSelection Nj, and add BG
//                             smearing Mi -->Nj'             purity and eff of eventSelection
//                                                          in Reco bin j for int and inc sample
//
//-------------------------------------------------------
//need to go back the steps
// 1) remove BG that is not Signal, do not care about smearing, vector of purity/efficiency for the reco bin
// 2) need to adress unfolding --> different Macro
// //--------------------------------------------------------

int eSliceMethod_doEffPur_evSel(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");
   

   string output_name = "effPur_eventSelection_bin_" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //Selected Process and RecoE Int and Inc Histos
   //THIS is what I get from DATA
   TH1D* h_recoE_selAbs_int = new TH1D("h_recoE_selAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_inc_initE = new TH1D("h_recoE_selPion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_inc_interE = new TH1D("h_recoE_selPion_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_incident = new TH1D("h_recoE_selPion_incident", "", nBin_int, eEnd, eStart);
   
   TH1D* h_recoE_selAbs_trueAbs_int = new TH1D("h_recoE_selAbs_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_trueTotInel_int = new TH1D("h_recoE_selPion_trueTotInel_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_truePion_inc_initE = new TH1D("h_recoE_selPion_truePion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_truePion_inc_interE = new TH1D("h_recoE_selPion_truePion_inc_interE", "", nBin_int, eEnd, eStart);
   
   TH1D* h_recoE_postPandora_truePion_inc_initE = new TH1D("h_recoE_postPandora_truePion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_postPandora_truePion_inc_interE = new TH1D("h_recoE_postPandora_truePion_inc_interE", "", nBin_int, eEnd, eStart);

   TH1D* h_recoE_postPandora_trueAbs_int = new TH1D("h_recoE_postPandora_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_postPandora_trueTotInel_int = new TH1D("h_recoE_postPandora_trueTotInel_int", "", nBin_int, eEnd, eStart);

  
   //From evSel --> back to Reco MC Nj --> Nj'
   TH1D* h_pur_removeBG_abs_int = new TH1D("h_pur_removeBG_abs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_abs_int = new TH1D("h_eff_eventSel_abs_int", "", nBin_int, eEnd, eStart);
   
   TH1D* h_pur_removeBG_totInel_int = new TH1D("h_pur_removeBG_totInel_int", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_totInel_int = new TH1D("h_eff_eventSel_totInel_int", "", nBin_int, eEnd, eStart);

   TH1D* h_pur_removeBG_inc_initE = new TH1D("h_pur_removeBG_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_inc_initE = new TH1D("h_eff_eventSel_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_pur_removeBG_inc_interE = new TH1D("h_pur_removeBG_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_inc_interE = new TH1D("h_eff_eventSel_inc_interE", "", nBin_int, eEnd, eStart);

   auto frame = inputFrame      
      .Define("true_initKE", "true_firstEntryIncident")
      .Define("true_interKE", "true_interactingKE_fromLength")
      .Filter("true_beam_endZ > 0");



   //Build the True Process and TrueE Int and Inc Histograms that we need to compare unsmeared things to
   //
   //all available after beamCuts
   auto mc_preReco_allPions = frame.Filter("true_beam_PDG == 211");
   auto eventSel_post_pandoraReco = frame.Filter("primary_isBeamType");
   //selected incident Pions & selected absorption
   //make sure they don't have reco initE == reco InterE bc in Selection they will be rejected as I can't fill my histos with such evnts
   auto eventSel_incidentPion = frame.Filter("selected_incidentPion");
   auto eventSel_abs = frame.Filter("selected_abs");

   //------------------------------------------------------------
   //
   //Create the selected and Reco histos that need to be unsmeared back
   //
   //
   //------------------------------------------------------------
   //
   //for this already removed recoInitE == recoInterE
   eventSel_incidentPion
      .Foreach( [h_recoE_selPion_inc_initE, h_recoE_selPion_inc_interE] (double init_KE, double inter_KE) { 

            fill_initE_interE( h_recoE_selPion_inc_initE, h_recoE_selPion_inc_interE, init_KE, inter_KE);
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});


   build_incidentHist( h_recoE_selPion_inc_initE, h_recoE_selPion_inc_interE, h_recoE_selPion_incident );
   h_recoE_selPion_incident->Write();

   //for this already removed recoInitE == recoInterE
   eventSel_abs
      .Foreach( [h_recoE_selAbs_int] (double init_KE, double inter_KE) { 

               fill_interacting( h_recoE_selAbs_int, init_KE, inter_KE);
            }
            ,{"reco_firstEntryIncident","reco_interactingKE"});

   h_recoE_selPion_inc_interE->Write();
   h_recoE_selPion_inc_initE->Write();
   h_recoE_selAbs_int->Write();

   //------------------------------------------------------------
   //
   //Start creating the efficiencies and purities to subtract the BG and account 
   //for lost events from eventSelection
   //
   //------------------------------------------------------------

   //True Pions available after beamCuts
   eventSel_post_pandoraReco
      //Events where trueInitE == trueInterE should not be available for this selection they are not a signal
      .Filter("true_beam_PDG == 211")
      .Foreach( [h_recoE_postPandora_truePion_inc_initE, h_recoE_postPandora_truePion_inc_interE] (double init_KE, double inter_KE) { 

            fill_initE_interE( h_recoE_postPandora_truePion_inc_initE, h_recoE_postPandora_truePion_inc_interE, init_KE, inter_KE);
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});


   //True Pions in selected incident
   eventSel_incidentPion
      //Events where trueInitE == trueInterE should not be available for this selection they are not a signal
      .Filter("true_beam_PDG == 211")
      .Foreach( [h_recoE_selPion_truePion_inc_initE, h_recoE_selPion_truePion_inc_interE] (double init_KE, double inter_KE) { 

            fill_initE_interE( h_recoE_selPion_truePion_inc_initE, h_recoE_selPion_truePion_inc_interE, init_KE, inter_KE);
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});

   //True Abs available after beamCuts
   eventSel_post_pandoraReco
      //Events where trueInitE == trueInterE should not be available for this selection they are not a signal
      .Filter("true_absSignal")
      .Foreach( [h_recoE_postPandora_trueAbs_int] (double init_KE, double inter_KE) { 

               fill_interacting( h_recoE_postPandora_trueAbs_int, init_KE, inter_KE);
            }
            ,{"reco_firstEntryIncident","reco_interactingKE"});
   //True Abs available after Selection
   eventSel_abs
      //Events where trueInitE == trueInterE should not be available for this selection they are not a signal
      .Filter("true_absSignal")
      .Foreach( [h_recoE_selAbs_trueAbs_int] (double init_KE, double inter_KE) { 

               fill_interacting( h_recoE_selAbs_trueAbs_int, init_KE, inter_KE);
            }
            ,{"reco_firstEntryIncident","reco_interactingKE"});

   //True TotInel available after beamCuts
   eventSel_post_pandoraReco
      //Events where trueInitE == trueInterE should not be available for this selection they are not a signal
      .Filter("true_primPionInel")
      .Foreach( [h_recoE_postPandora_trueTotInel_int] (double init_KE, double inter_KE) { 

               fill_interacting( h_recoE_postPandora_trueTotInel_int, init_KE, inter_KE);
            }
            ,{"reco_firstEntryIncident","reco_interactingKE"});
   
   //True TotInel available after Selection
   eventSel_incidentPion
      //Events where trueInitE == trueInterE should not be available for this selection they are not a signal
      .Filter("true_primPionInel")
      .Foreach( [h_recoE_selPion_trueTotInel_int] (double init_KE, double inter_KE) { 

               fill_interacting( h_recoE_selPion_trueTotInel_int, init_KE, inter_KE);
            }
            ,{"reco_firstEntryIncident","reco_interactingKE"});


   //------------------------------------------------------------
   //
   //  Compute eff and Purities for first unsmearing Step
   //  pur = true / all selected
   //  eff = true / all available at beginning = pur*all selected / all available at beginning
   //
   //
   //------------------------------------------------------------
   //
   h_pur_removeBG_inc_initE->Divide( h_recoE_selPion_truePion_inc_initE, h_recoE_selPion_inc_initE );
   h_eff_eventSel_inc_initE->Multiply( h_pur_removeBG_inc_initE, h_recoE_selPion_inc_initE );
   h_eff_eventSel_inc_initE->Divide( h_recoE_postPandora_truePion_inc_initE );

   h_pur_removeBG_inc_interE->Divide( h_recoE_selPion_truePion_inc_interE, h_recoE_selPion_inc_interE );
   h_eff_eventSel_inc_interE->Multiply( h_pur_removeBG_inc_interE, h_recoE_selPion_inc_interE );
   h_eff_eventSel_inc_interE->Divide( h_recoE_postPandora_truePion_inc_interE );

   h_pur_removeBG_abs_int->Divide( h_recoE_selAbs_trueAbs_int, h_recoE_selAbs_int );
   h_eff_eventSel_abs_int->Multiply( h_pur_removeBG_abs_int, h_recoE_selAbs_int );
   h_eff_eventSel_abs_int->Divide( h_recoE_postPandora_trueAbs_int );

   h_pur_removeBG_totInel_int->Divide( h_recoE_selPion_trueTotInel_int, h_recoE_selPion_inc_interE );
   h_eff_eventSel_totInel_int->Multiply( h_pur_removeBG_totInel_int, h_recoE_selPion_inc_interE );
   h_eff_eventSel_totInel_int->Divide( h_recoE_postPandora_trueTotInel_int );

   h_pur_removeBG_inc_initE->Write();
   h_eff_eventSel_inc_initE->Write();
   h_pur_removeBG_inc_interE->Write();
   h_eff_eventSel_inc_interE->Write();

   h_pur_removeBG_abs_int->Write();
   h_eff_eventSel_abs_int->Write();

   h_pur_removeBG_totInel_int->Write();
   h_eff_eventSel_totInel_int->Write();
   

   return 0;
}

