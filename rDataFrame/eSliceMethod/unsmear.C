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
//
//macro to test the unsmearing operations
//
//Try to unsmear interacting and incident histogram from RecoE and Selection into trueE and True Process
//
//MC True --------------------> MC Reco ------------------> Selected Interaction, Selected Incident
// 
//correct for Reco Eff         Pandora Reco                 EventSelection Nj, and add BG
//                             smearing Mi -->Nj'             purity and eff of eventSelection
//                                                          in Reco bin j for int and inc sample
//
//-------------------------------------------------------
//need to go back the steps
// 1) remove BG that is not Signal, do not care about smearing, vector of purity/efficiency for the reco bin
// 2) need to apply inverse of the smearing matrix (true to reco)
// 3) need to correct for Reconstruction inefficiency
//
// Validate by having the true Process Int and Inc in trueE all with true_beam_endZ
// Try to unsmear the selected Events in recoE to the above mentioned
//--------------------------------------------------------

int unsmear(const string mcFilepath, bool doEffPur, bool doSmearing, bool doRecoEff){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");

   //Test getting Gauss Inverse
   //TFile f1("smearMatrixGauss_20MeV.root");
   //TH2D *invGauss = (TH2D*)f1.Get("inverseMatrix_gaussFit_interacting");
   //f1.Close();
   

   //string output_name = "unsmear_halfMC_" + std::to_string((int) bin_size_int) + "MeV.root";
   string output_name = "unsmear_" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //TrueProcess and TrueE Int and Inc Histos before Pandora Reco
   TH1D* h_trueE_prePandora_trueAbs_int = new TH1D("h_trueE_prePandora_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_prePandora_trueTotInel_int = new TH1D("h_trueE_prePandora_trueTotInel_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_prePandora_truePion_inc_initE = new TH1D("h_trueE_prePandora_truePion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_prePandora_truePion_inc_interE = new TH1D("h_trueE_prePandora_truePion_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_prePandora_truePion_incident = new TH1D("h_trueE_prePandora_truePion_incident", "", nBin_int, eEnd, eStart);
   
   //TrueProcess and TrueE Int and Inc Histos after Pandora Reco
   TH1D* h_trueE_postPandora_trueAbs_int = new TH1D("h_trueE_postPandora_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_postPandora_trueTotInel_int = new TH1D("h_trueE_postPandora_trueTotInel_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_postPandora_truePion_inc_initE = new TH1D("h_trueE_postPandora_truePion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_postPandora_truePion_inc_interE = new TH1D("h_trueE_postPandora_truePion_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_postPandora_truePion_incident = new TH1D("h_trueE_postPandora_truePion_incident", "", nBin_int, eEnd, eStart);

   //Selected Process and trueE Int and Inc Histos
   TH1D* h_trueE_selAbs_int = new TH1D("h_trueE_selAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_selPion_inc_initE = new TH1D("h_trueE_selPion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_selPion_inc_interE = new TH1D("h_trueE_selPion_inc_interE", "", nBin_int, eEnd, eStart);
   //TH1D* h_trueE_selPion_incident = new TH1D("h_trueE_selPion_incident", "", nBin_int, eEnd, eStart);

   //Selected Process and trueE Int and Inc Histos
   TH1D* h_trueE_trueAbs_selAbs_int = new TH1D("h_trueE_trueAbs_selAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_trueTotInel_selPion_int = new TH1D("h_trueE_trueTotInel_selPion_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_selPion_inc_initE = new TH1D("h_trueE_truePion_selPion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_selPion_inc_interE = new TH1D("h_trueE_truePion_selPion_inc_interE", "", nBin_int, eEnd, eStart);
   //TH1D* h_trueE_truePion_selPion_incident = new TH1D("h_trueE_truePion_selPion_incident", "", nBin_int, eEnd, eStart);

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

   //Unsmeared Histo, compare to h_trueE_trueProc
   TH1D* h_help_unsmear_abs_int = new TH1D("h_help_unsmear_abs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_help_unsmear_totInel_int = new TH1D("h_help_unsmear_totInel_int", "", nBin_int, eEnd, eStart);
   TH1D* h_help_unsmear_inc_initE = new TH1D("h_help_unsmear_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_help_unsmear_inc_interE = new TH1D("h_help_unsmear_inc_interE", "", nBin_int, eEnd, eStart);

   TH1D* h_unsmeared_abs_int = new TH1D("h_unsmeared_abs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_unsmeared_totInel_int = new TH1D("h_unsmeared_totInel_int", "", nBin_int, eEnd, eStart);
   TH1D* h_unsmeared_inc_initE = new TH1D("h_unsmeared_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_unsmeared_inc_interE = new TH1D("h_unsmeared_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_unsmeared_incident = new TH1D("h_unsmeared_incident", "", nBin_int, eEnd, eStart);
   
   TH1D* h_corrPurEffEvSel_abs_int = new TH1D("h_corrPurEffEvSel_abs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_corrPurEffEvSel_totInel_int = new TH1D("h_corrPurEffEvSel_totInel_int", "", nBin_int, eEnd, eStart);
   TH1D* h_corrPurEffEvSel_inc_initE = new TH1D("h_corrPurEffEvSel_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_corrPurEffEvSel_inc_interE = new TH1D("h_corrPurEffEvSel_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_corrPurEffEvSel_incident = new TH1D("h_corrPurEffEvSel_incident", "", nBin_int, eEnd, eStart);
   
   //From evSel --> back to Reco MC Nj --> Nj'
   TH1D* h_pur_removeBG_abs_int = new TH1D("h_pur_removeBG_abs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_abs_int = new TH1D("h_eff_eventSel_abs_int", "", nBin_int, eEnd, eStart);
   
   TH1D* h_pur_removeBG_totInel_int = new TH1D("h_pur_removeBG_totInel_int", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_totInel_int = new TH1D("h_eff_eventSel_totInel_int", "", nBin_int, eEnd, eStart);

   TH1D* h_pur_removeBG_inc_initE = new TH1D("h_pur_removeBG_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_inc_initE = new TH1D("h_eff_eventSel_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_pur_removeBG_inc_interE = new TH1D("h_pur_removeBG_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_inc_interE = new TH1D("h_eff_eventSel_inc_interE", "", nBin_int, eEnd, eStart);

   TH1D* h_trueE_pur_removeBG_abs_int = new TH1D("h_trueE_pur_removeBG_abs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_eff_eventSel_abs_int = new TH1D("h_trueE_eff_eventSel_abs_int", "", nBin_int, eEnd, eStart);
   
   TH1D* h_trueE_pur_removeBG_totInel_int = new TH1D("h_trueE_pur_removeBG_totInel_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_eff_eventSel_totInel_int = new TH1D("h_trueE_eff_eventSel_totInel_int", "", nBin_int, eEnd, eStart);

   TH1D* h_trueE_pur_removeBG_inc_initE = new TH1D("h_trueE_pur_removeBG_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_eff_eventSel_inc_initE = new TH1D("h_trueE_eff_eventSel_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_pur_removeBG_inc_interE = new TH1D("h_trueE_pur_removeBG_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_eff_eventSel_inc_interE = new TH1D("h_trueE_eff_eventSel_inc_interE", "", nBin_int, eEnd, eStart);

   TH1D* h_trueE_effPurCorr_abs_int = new TH1D("h_trueE_effPurCorr_abs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_effPurCorr_totInel_int = new TH1D("h_trueE_effPurCorr_totInel_int", "", nBin_int, eEnd, eStart);
   
   TH1D* h_trueE_effPurCorr_inc_initE = new TH1D("h_trueE_effPurCorr_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_effPurCorr_inc_interE = new TH1D("h_trueE_effPurCorr_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_effPurCorr_incident = new TH1D("h_trueE_effPurCorr_incident", "", nBin_int, eEnd, eStart);

   TH1D* h_trueE_eff_reconstruction_inc_initE = new TH1D("h_trueE_eff_reconstruction_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_eff_reconstruction_inc_interE = new TH1D("h_trueE_eff_reconstruction_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_eff_reconstruction_abs_int = new TH1D("h_trueE_eff_reconstruction_abs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_eff_reconstruction_totInel_int = new TH1D("h_trueE_eff_reconstruction_totInel_int", "", nBin_int, eEnd, eStart);

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

   
   //----------------------------------
   //Available Pions prior to Reconstruction --> will be used to correct for reco efficiency
   //------------------------------------------------------------

   mc_preReco_allPions //already filtered for pions
      .Foreach( [h_trueE_prePandora_truePion_inc_initE, h_trueE_prePandora_truePion_inc_interE] (double init_KE, double inter_KE) { 

            fill_initE_interE( h_trueE_prePandora_truePion_inc_initE, h_trueE_prePandora_truePion_inc_interE, init_KE, inter_KE);

            }
            ,{"true_initKE", "true_interKE"});

   build_incidentHist( h_trueE_prePandora_truePion_inc_initE, h_trueE_prePandora_truePion_inc_interE, h_trueE_prePandora_truePion_incident );
   h_trueE_prePandora_truePion_inc_initE->Write();
   h_trueE_prePandora_truePion_inc_interE->Write();
   h_trueE_prePandora_truePion_incident->Write();

   mc_preReco_allPions
      .Filter("true_absSignal") //should also take into account pions that decay

      .Foreach( [h_trueE_prePandora_trueAbs_int] (double init_KE, double inter_KE) { 

            fill_interacting( h_trueE_prePandora_trueAbs_int, init_KE, inter_KE);

            }
            ,{"true_initKE", "true_interKE"});

   h_trueE_prePandora_trueAbs_int->Write();

   mc_preReco_allPions
      .Filter("true_primPionInel") 
      .Foreach( [h_trueE_prePandora_trueTotInel_int] (double init_KE, double inter_KE) { 

            fill_interacting( h_trueE_prePandora_trueTotInel_int, init_KE, inter_KE);
            }
            ,{"true_initKE", "true_interKE"});

   h_trueE_prePandora_trueTotInel_int->Write();

   //----------------------------------
   //Fill the Histos with true Process and trueEbin Incident and interacting.
   //------------------------------------------------------------
   //Incident Histo InitE
   eventSel_post_pandoraReco
      .Filter("true_beam_PDG == 211")
      .Foreach( [h_trueE_postPandora_truePion_inc_initE, h_trueE_postPandora_truePion_inc_interE] (double init_KE, double inter_KE) { 
            
            fill_initE_interE(h_trueE_postPandora_truePion_inc_initE, h_trueE_postPandora_truePion_inc_interE, init_KE, inter_KE);
            
            }
            ,{"true_initKE", "true_interKE"});

   build_incidentHist( h_trueE_postPandora_truePion_inc_initE, h_trueE_postPandora_truePion_inc_interE, h_trueE_postPandora_truePion_incident );
   h_trueE_postPandora_truePion_inc_initE->Write();
   h_trueE_postPandora_truePion_inc_interE->Write();
   h_trueE_postPandora_truePion_incident->Write();

   //Interacting Histo ABS
   eventSel_post_pandoraReco
      .Filter("true_absSignal") //should also take into account pions that decay

      .Foreach( [h_trueE_postPandora_trueAbs_int] (double init_KE, double inter_KE) { 

               fill_interacting( h_trueE_postPandora_trueAbs_int, init_KE, inter_KE);
            }
            ,{"true_initKE", "true_interKE"});


   h_trueE_postPandora_trueAbs_int->Write();

   //Interacting Histo TOTINEL
   eventSel_post_pandoraReco
      .Filter("true_primPionInel") //should also take into account pions that decay

      .Foreach( [h_trueE_postPandora_trueTotInel_int] (double init_KE, double inter_KE) { 

               fill_interacting( h_trueE_postPandora_trueTotInel_int, init_KE, inter_KE);
            }
            ,{"true_initKE", "true_interKE"});


   h_trueE_postPandora_trueTotInel_int->Write();
   //------------------------------------------------------------
   //
   //Create the Reconstruction efficiencies in trueE
   //
   //------------------------------------------------------------
   h_trueE_eff_reconstruction_inc_initE->Divide( h_trueE_postPandora_truePion_inc_initE , h_trueE_prePandora_truePion_inc_initE );
   h_trueE_eff_reconstruction_inc_interE->Divide( h_trueE_postPandora_truePion_inc_interE , h_trueE_prePandora_truePion_inc_interE );
   
   h_trueE_eff_reconstruction_abs_int->Divide( h_trueE_postPandora_trueAbs_int , h_trueE_prePandora_trueAbs_int );
   h_trueE_eff_reconstruction_totInel_int->Divide( h_trueE_postPandora_trueTotInel_int , h_trueE_prePandora_trueTotInel_int );

   h_trueE_eff_reconstruction_inc_initE->Write();
   h_trueE_eff_reconstruction_inc_interE->Write();
   h_trueE_eff_reconstruction_abs_int->Write();
   h_trueE_eff_reconstruction_totInel_int->Write();
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


   //True Pions trueE in Pion Selection
   //
   eventSel_incidentPion
      .Filter("true_beam_PDG == 211")
      .Foreach( [h_trueE_truePion_selPion_inc_initE, h_trueE_truePion_selPion_inc_interE] (double init_KE, double inter_KE) { 

            fill_initE_interE( h_trueE_truePion_selPion_inc_initE, h_trueE_truePion_selPion_inc_interE, init_KE, inter_KE);
            
            }
            ,{"true_initKE", "true_interKE"});

   //Selected Pions trueE in Pion Selection
   //
   eventSel_incidentPion
      .Foreach( [h_trueE_selPion_inc_initE, h_trueE_selPion_inc_interE] (double init_KE, double inter_KE) { 
         
            fill_initE_interE( h_trueE_selPion_inc_initE, h_trueE_selPion_inc_interE, init_KE, inter_KE);
            }
            ,{"true_initKE", "true_interKE"});

   //True Abs trueE in Pion Selection
   //
   eventSel_abs
      .Filter("true_absSignal")      
      .Foreach( [h_trueE_trueAbs_selAbs_int] (double init_KE, double inter_KE) { 

               fill_interacting( h_trueE_trueAbs_selAbs_int, init_KE, inter_KE);
            }
            ,{"true_initKE", "true_interKE"});

   //Selected Abs trueE in Abs Selection
   //
   eventSel_abs
      .Foreach( [h_trueE_selAbs_int] (double init_KE, double inter_KE) { 

               fill_interacting( h_trueE_selAbs_int, init_KE, inter_KE);
            }
            ,{"true_initKE", "true_interKE"});

   //True TotINEL trueE in Pion Selection
   //
   eventSel_incidentPion
      .Filter("true_primPionInel")      
      .Foreach( [h_trueE_trueTotInel_selPion_int] (double init_KE, double inter_KE) { 

               fill_interacting( h_trueE_trueTotInel_selPion_int, init_KE, inter_KE);
            }
            ,{"true_initKE", "true_interKE"});


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

   h_trueE_pur_removeBG_inc_initE->Divide( h_trueE_truePion_selPion_inc_initE, h_trueE_selPion_inc_initE);
   h_trueE_eff_eventSel_inc_initE->Multiply( h_trueE_pur_removeBG_inc_initE, h_trueE_selPion_inc_initE);
   h_trueE_eff_eventSel_inc_initE->Divide( h_trueE_postPandora_truePion_inc_initE );

   h_trueE_pur_removeBG_inc_interE->Divide( h_trueE_truePion_selPion_inc_interE, h_trueE_selPion_inc_interE);
   h_trueE_eff_eventSel_inc_interE->Multiply( h_trueE_pur_removeBG_inc_interE, h_trueE_selPion_inc_interE);
   h_trueE_eff_eventSel_inc_interE->Divide( h_trueE_postPandora_truePion_inc_interE );

   h_trueE_pur_removeBG_abs_int->Divide( h_trueE_trueAbs_selAbs_int, h_trueE_selAbs_int);
   h_trueE_eff_eventSel_abs_int->Multiply( h_trueE_pur_removeBG_abs_int, h_trueE_selAbs_int);
   h_trueE_eff_eventSel_abs_int->Divide( h_trueE_postPandora_trueAbs_int );

   h_trueE_pur_removeBG_totInel_int->Divide( h_trueE_trueTotInel_selPion_int, h_trueE_selPion_inc_interE);
   h_trueE_eff_eventSel_totInel_int->Multiply( h_trueE_pur_removeBG_totInel_int, h_trueE_selPion_inc_interE);
   h_trueE_eff_eventSel_totInel_int->Divide( h_trueE_postPandora_trueTotInel_int );
   
   h_trueE_pur_removeBG_inc_initE->Write();
   h_trueE_eff_eventSel_inc_initE->Write();
   h_trueE_pur_removeBG_inc_interE->Write();
   h_trueE_eff_eventSel_inc_interE->Write();

   h_trueE_pur_removeBG_abs_int->Write();
   h_trueE_eff_eventSel_abs_int->Write();

   h_trueE_pur_removeBG_totInel_int->Write();
   h_trueE_eff_eventSel_totInel_int->Write();

   //------------------------------------------------------------
   //
   //  Now build Inverse Smearing Matrix to translate back to recoE j --> trueE i
   //
   //  matrix Tij*Nj=Mi
   //  Tij has true on y and reco on x axis, normalised on reco
   //------------------------------------------------------------
   //

   TH2D* h2_smearing_abs_int = new TH2D("h2_smearing_abs_int", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_smearing_totInel_int = new TH2D("h2_smearing_totInel_int", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_smearing_incident_initE = new TH2D("h2_smearing_incident_initE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_smearing_incident_interE = new TH2D("h2_smearing_incident_interE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   h2_smearing_abs_int->SetTitle("Smearing Matrix for trueAbs in Sample after Pandora Reco; reco Energy [MeV]; true Energy [MeV]");
   h2_smearing_totInel_int->SetTitle("Smearing Matrix for trueTotInel in Sample after Pandora Reco; reco Energy [MeV]; true Energy [MeV]");
   h2_smearing_incident_initE->SetTitle("Smearing Matrix for truePi initialE in Sample after Pandora Reco; reco Energy [MeV]; true Energy [MeV]");
   h2_smearing_incident_interE->SetTitle("Smearing Matrix for truePi interE in Sample after Pandora Reco; reco Energy [MeV]; true Energy [MeV]");

   TH2D* h2_inverse_smearing_abs_int = new TH2D("h2_inverse_smearing_abs_int", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_inverse_smearing_totInel_int = new TH2D("h2_inverse_smearing_totInel_int", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_inverse_smearing_incident_initE = new TH2D("h2_inverse_smearing_incident_initE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_inverse_smearing_incident_interE = new TH2D("h2_inverse_smearing_incident_interE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   h2_inverse_smearing_abs_int->SetTitle("Inverse smearing Matrix for trueAbs in Sample after Pandora Reco; true Energy [MeV]; reco Energy [MeV]");
   h2_inverse_smearing_totInel_int->SetTitle("Inverse smearing Matrix for trueTotInel in Sample after Pandora Reco; true Energy [MeV]; reco Energy [MeV]");
   h2_inverse_smearing_incident_initE->SetTitle("Inverse smearing Matrix for truePi initialE in Sample after Pandora Reco; true Energy [MeV]; reco Energy [MeV]");
   h2_inverse_smearing_incident_interE->SetTitle("Inverse smearing Matrix for truePi interE in Sample after Pandora Reco; true Energy [MeV]; reco Energy [MeV]");

   //Fill smearing matrix with entry 1 on the diagonal to avoid them being "uninvertible"
   for(int i = 1; i <= nBin_int; i++){
      h2_smearing_incident_initE->SetBinContent(i,i,1);
      h2_smearing_incident_interE->SetBinContent(i,i,1);
      h2_smearing_abs_int->SetBinContent(i,i,1);
      h2_smearing_totInel_int->SetBinContent(i,i,1);
   };
   //========================================================
   //Build the smearing Incident Histogram initE
   //IGNORE events that have true initE-bin == interE-bin
   //---------

   eventSel_post_pandoraReco 
      .Filter("true_beam_PDG == 211")
      .Foreach( [h2_smearing_incident_initE, h2_smearing_incident_interE] (double true_initE, double true_interE, double reco_initE, double reco_interE ){

            fill_smearing_matrix(h2_smearing_incident_initE, true_initE, true_interE, reco_initE, reco_interE, false);
            
            fill_smearing_matrix(h2_smearing_incident_interE, true_initE, true_interE, reco_initE, reco_interE, true);

            }
            ,{"true_initKE", "true_interKE", "reco_firstEntryIncident", "reco_interactingKE"}); 


   //NORMALISATION of initE and interE
   ////
   
   normalise_smearing(h2_smearing_incident_initE);
   normalise_smearing(h2_smearing_incident_interE);

   h2_smearing_incident_initE->Sumw2(0);
   h2_smearing_incident_initE->Write();

   h2_smearing_incident_interE->Sumw2(0);
   h2_smearing_incident_interE->Write();

   invert_smearing(h2_smearing_incident_initE, h2_inverse_smearing_incident_initE );
   h2_inverse_smearing_incident_initE->Write();
   
   invert_smearing(h2_smearing_incident_interE, h2_inverse_smearing_incident_interE );
   h2_inverse_smearing_incident_interE->Write();



   //========================================================
   //Build the smearing Interacting Histogram
   //---------
   //eventSel_abs
   eventSel_post_pandoraReco
      .Filter("true_absSignal")
      //.Filter("true_absSignal && true_reco_initE_eq_interE")
      .Foreach( [h2_smearing_abs_int] (double true_initE, double true_interE, double reco_initE, double reco_interE ){

               fill_smearing_matrix( h2_smearing_abs_int, true_initE, true_interE, reco_initE, reco_interE, true );
            }
            ,{"true_initKE", "true_interKE", "reco_firstEntryIncident", "reco_interactingKE"}); 

   //NORMALISATION Normalise to recoE column
   //Go through true columns and normalise the entries 

   normalise_smearing(h2_smearing_abs_int);

   h2_smearing_abs_int->Sumw2(0);
   h2_smearing_abs_int->Write();

   invert_smearing(h2_smearing_abs_int, h2_inverse_smearing_abs_int );
   h2_inverse_smearing_abs_int->Write();

   eventSel_post_pandoraReco
      .Filter("true_primPionInel")
      //.Filter("true_absSignal && true_reco_initE_eq_interE")
      .Foreach( [h2_smearing_totInel_int] (double true_initE, double true_interE, double reco_initE, double reco_interE ){

               fill_smearing_matrix( h2_smearing_totInel_int, true_initE, true_interE, reco_initE, reco_interE, true );
            }
            ,{"true_initKE", "true_interKE", "reco_firstEntryIncident", "reco_interactingKE"}); 

   //NORMALISATION Normalise to recoE column
   //Go through true columns and normalise the entries 

   normalise_smearing(h2_smearing_totInel_int);

   h2_smearing_totInel_int->Sumw2(0);
   h2_smearing_totInel_int->Write();

   invert_smearing(h2_smearing_totInel_int, h2_inverse_smearing_totInel_int );
   h2_inverse_smearing_totInel_int->Write();
 
   //=====================================================
   //------------------------------------------------------
   //Start Unsmearing Incident Init E and InterE
   //------------------------------------------------------
   //
   //    Remove the BGs and eventSelection inefficiency
   //    multiply by purity and divide by efficency
   if(doEffPur){
   h_help_unsmear_inc_initE->Multiply(h_recoE_selPion_inc_initE, h_pur_removeBG_inc_initE );
   h_help_unsmear_inc_initE->Divide( h_eff_eventSel_inc_initE );

   h_help_unsmear_inc_interE->Multiply(h_recoE_selPion_inc_interE, h_pur_removeBG_inc_interE );
   h_help_unsmear_inc_interE->Divide( h_eff_eventSel_inc_interE );

   }
   else{
      h_help_unsmear_inc_initE = h_recoE_selPion_inc_initE;
      h_help_unsmear_inc_interE = h_recoE_selPion_inc_interE;
   }
   //   Undo the smearing by multiplying Tij*Nj, so every row of true percentage with every entry of help_unsmeared
   //   Row's are i(y), columns are j(x) 
   //
   //Loop through smearing matrix Incident
   

   if(doSmearing){
   
   for(int i = 1; i <= nBin_int; i++){
      double help_sum_initE = 0, help_sum_interE = 0;
      //rows
      for(int j = 1; j <= nBin_int; j++){
         help_sum_initE += h_help_unsmear_inc_initE->GetBinContent(j)*h2_inverse_smearing_incident_initE->GetBinContent( i, j);
         help_sum_interE += h_help_unsmear_inc_interE->GetBinContent(j)*h2_inverse_smearing_incident_interE->GetBinContent( i, j);
      }; 
      h_unsmeared_inc_initE->SetBinContent( i , help_sum_initE);
      h_unsmeared_inc_interE->SetBinContent( i , help_sum_interE);
   };
   
   /*for(int i = 1; i <= nBin_int; i++){
      double help_sum_initE = 0, help_sum_interE = 0;
      //rows
      for(int j = 1; j <= nBin_int; j++){
         help_sum_initE += h_help_unsmear_inc_initE->GetBinContent(j)*h2_inverse_smearing_incident_initE->GetBinContent( i, j);
         help_sum_interE += h_help_unsmear_inc_interE->GetBinContent(j)*invGauss->GetBinContent( i, j);
      }; 
      h_unsmeared_inc_initE->SetBinContent( i , help_sum_initE);
      h_unsmeared_inc_interE->SetBinContent( i , help_sum_interE);
   };*/
   }
   else{
      h_unsmeared_inc_initE = h_help_unsmear_inc_initE;
      h_unsmeared_inc_interE = h_help_unsmear_inc_interE;

   }

   h_unsmeared_inc_initE->Write();
   h_unsmeared_inc_interE->Write();

   if(doRecoEff){
      h_corrPurEffEvSel_inc_initE->Divide( h_unsmeared_inc_initE , h_trueE_eff_reconstruction_inc_initE );
      h_corrPurEffEvSel_inc_interE->Divide( h_unsmeared_inc_interE , h_trueE_eff_reconstruction_inc_interE );
   }
   else{
      h_corrPurEffEvSel_inc_initE = h_unsmeared_inc_initE;
      h_corrPurEffEvSel_inc_interE = h_unsmeared_inc_interE;
   }

   h_corrPurEffEvSel_inc_initE->Write();
   h_corrPurEffEvSel_inc_interE->Write();

   //=====================================================
   //    Rebuild from initial and interacting Distribution the full incident unsmeared histogram
   //=====================================================

   build_incidentHist( h_unsmeared_inc_initE, h_help_unsmear_inc_interE, h_unsmeared_incident );
   h_unsmeared_incident->Write();

   build_incidentHist( h_corrPurEffEvSel_inc_initE, h_help_unsmear_inc_interE, h_corrPurEffEvSel_incident );
   h_corrPurEffEvSel_incident->Write();

   //------------------------------------------------------
   //Start Unsmearing Interacting
   //------------------------------------------------------
   //
   //    Remove the BGs and eventSelection inefficiency
   //    multiply by purity and divide by efficency
   if(doEffPur){ 
   h_help_unsmear_abs_int->Multiply(h_recoE_selAbs_int, h_pur_removeBG_abs_int );
   h_help_unsmear_abs_int->Divide( h_eff_eventSel_abs_int );
   
   h_help_unsmear_totInel_int->Multiply(h_recoE_selPion_inc_interE, h_pur_removeBG_totInel_int );
   h_help_unsmear_totInel_int->Divide( h_eff_eventSel_totInel_int );
   }
   else{
      h_help_unsmear_abs_int = h_recoE_selAbs_int;
      h_help_unsmear_totInel_int = h_recoE_selPion_inc_interE;
   }
   //   Undo the smearing by multiplying Tij*Nj, so every row of true percentage with every entry of help_unsmeared
   //   Row's are i(y), columns are j(x) 
   //
   //Loop through smearing matrix Interacting
   if(doSmearing){
   for(int i = 1; i <= nBin_int; i++){
      double help_sum_abs = 0;
      double help_sum_totInel = 0;
      //rows
      for(int j = 1; j <= nBin_int; j++){
         help_sum_abs += h_help_unsmear_abs_int->GetBinContent(j)*h2_inverse_smearing_abs_int->GetBinContent( i, j);
         help_sum_totInel += h_help_unsmear_totInel_int->GetBinContent(j)*h2_inverse_smearing_totInel_int->GetBinContent( i, j);
      }; 
      h_unsmeared_abs_int->SetBinContent( i , help_sum_abs);
      h_unsmeared_totInel_int->SetBinContent( i , help_sum_totInel);
   };
   }
   else {
      h_unsmeared_abs_int = h_help_unsmear_abs_int;
      h_unsmeared_totInel_int = h_help_unsmear_totInel_int;
   }
   h_unsmeared_abs_int->Write();
   h_unsmeared_totInel_int->Write();

   if(doRecoEff){
      h_corrPurEffEvSel_abs_int->Divide( h_unsmeared_abs_int, h_trueE_eff_reconstruction_abs_int );
      h_corrPurEffEvSel_totInel_int->Divide( h_unsmeared_totInel_int, h_trueE_eff_reconstruction_totInel_int );
   }
   else{
      h_corrPurEffEvSel_abs_int = h_unsmeared_abs_int;
      h_corrPurEffEvSel_totInel_int = h_unsmeared_totInel_int;
   }

   h_corrPurEffEvSel_abs_int->Write();
   h_corrPurEffEvSel_totInel_int->Write();

   //------------------------------------------------------
   //          CANVAS
   //------------------------------------------------------
   //

   TCanvas *c_comp_abs_int = new TCanvas("c_comp_abs_int", "");

   h_trueE_prePandora_trueAbs_int->SetLineColor(kBlue);
   h_trueE_prePandora_trueAbs_int->SetLineWidth(2);
   h_trueE_prePandora_trueAbs_int->SetMarkerColor(kBlue);
   h_trueE_prePandora_trueAbs_int->SetBarOffset(0.5);
   h_trueE_prePandora_trueAbs_int->GetXaxis()->SetTitle("KE (MeV)");
   
  // h_trueE_postPandora_trueAbs_int->SetLineColor(kBlue);
  // h_trueE_postPandora_trueAbs_int->SetLineWidth(2);
  // h_trueE_postPandora_trueAbs_int->SetMarkerColor(kBlue);
  // h_trueE_postPandora_trueAbs_int->SetBarOffset(0.5);
  // h_trueE_postPandora_trueAbs_int->GetXaxis()->SetTitle("KE (MeV)");

   //h_unsmeared_abs_int->SetLineColor(kRed);
   //h_unsmeared_abs_int->SetLineWidth(2);
   //h_unsmeared_abs_int->SetFillColorAlpha( kRed, 0.2);
   //h_unsmeared_abs_int->SetMarkerColor(kRed);
   //h_unsmeared_abs_int->SetBarOffset(-0.5);
   
   h_corrPurEffEvSel_abs_int->SetLineColor(kRed);
   h_corrPurEffEvSel_abs_int->SetLineWidth(2);
   h_corrPurEffEvSel_abs_int->SetFillColorAlpha( kRed, 0.2);
   h_corrPurEffEvSel_abs_int->SetMarkerColor(kRed);
   h_corrPurEffEvSel_abs_int->SetBarOffset(-0.5);

   h_recoE_selAbs_int->SetLineStyle(2);
   h_recoE_selAbs_int->SetLineWidth(3);
   h_recoE_selAbs_int->SetFillColorAlpha(kBlack, 0.1);

   h_trueE_prePandora_trueAbs_int->Draw("HIST TEXT00");
   h_corrPurEffEvSel_abs_int->Draw("HIST SAME");
   h_recoE_selAbs_int->Draw("HIST SAME");

   auto legend_abs_int = new TLegend(0.1,0.7,0.48,0.9);
   legend_abs_int->AddEntry(h_trueE_prePandora_trueAbs_int,"True Absorption, trueE, before PandoraReco");
   legend_abs_int->AddEntry(h_corrPurEffEvSel_abs_int,"Unsmeared and Corrected for RecoEff Interacting");
   legend_abs_int->AddEntry(h_recoE_selAbs_int,"Selected Interacting Absorption");
   legend_abs_int->Draw();
   c_comp_abs_int->Write();

   TCanvas *c_comp_totInel_int = new TCanvas("c_comp_totInel_int", "");

   h_trueE_prePandora_trueTotInel_int->SetLineColor(kBlue);
   h_trueE_prePandora_trueTotInel_int->SetLineWidth(2);
   h_trueE_prePandora_trueTotInel_int->SetMarkerColor(kBlue);
   h_trueE_prePandora_trueTotInel_int->SetBarOffset(0.5);
   h_trueE_prePandora_trueTotInel_int->GetXaxis()->SetTitle("KE (MeV)");
   
  // h_trueE_postPandora_trueTotInel_int->SetLineColor(kBlue);
  // h_trueE_postPandora_trueTotInel_int->SetLineWidth(2);
  // h_trueE_postPandora_trueTotInel_int->SetMarkerColor(kBlue);
  // h_trueE_postPandora_trueTotInel_int->SetBarOffset(0.5);
  // h_trueE_postPandora_trueTotInel_int->GetXaxis()->SetTitle("KE (MeV)");

   //h_unsmeared_totInel_int->SetLineColor(kRed);
   //h_unsmeared_totInel_int->SetLineWidth(2);
   //h_unsmeared_totInel_int->SetFillColorAlpha( kRed, 0.2);
   //h_unsmeared_totInel_int->SetMarkerColor(kRed);
   //h_unsmeared_totInel_int->SetBarOffset(-0.5);
   
   h_corrPurEffEvSel_totInel_int->SetLineColor(kRed);
   h_corrPurEffEvSel_totInel_int->SetLineWidth(2);
   h_corrPurEffEvSel_totInel_int->SetFillColorAlpha( kRed, 0.2);
   h_corrPurEffEvSel_totInel_int->SetMarkerColor(kRed);
   h_corrPurEffEvSel_totInel_int->SetBarOffset(-0.5);

   h_recoE_selPion_inc_interE->SetLineStyle(2);
   h_recoE_selPion_inc_interE->SetLineWidth(3);
   h_recoE_selPion_inc_interE->SetFillColorAlpha(kBlack, 0.1);

   h_trueE_prePandora_trueTotInel_int->Draw("HIST TEXT00");
   h_corrPurEffEvSel_totInel_int->Draw("HIST SAME");
   h_recoE_selPion_inc_interE->Draw("HIST SAME");

   auto legend_totInel_int = new TLegend(0.1,0.7,0.48,0.9);
   legend_totInel_int->AddEntry(h_trueE_prePandora_trueTotInel_int,"True Total Inelastic, trueE, before PandoraReco");
   legend_totInel_int->AddEntry(h_corrPurEffEvSel_totInel_int,"Unsmeared and Corrected for RecoEff Interacting");
   legend_totInel_int->AddEntry(h_recoE_selPion_inc_interE,"Selected Interacting Total Inelastic");
   legend_totInel_int->Draw();
   c_comp_totInel_int->Write();

   TCanvas *c_comp_inc = new TCanvas("c_comp_inc", "");

   h_trueE_postPandora_truePion_incident->SetLineColor(kBlue);
   h_trueE_postPandora_truePion_incident->SetLineWidth(2);
   h_trueE_postPandora_truePion_incident->SetMarkerColor(kBlue);
   h_trueE_postPandora_truePion_incident->SetMarkerSize(0.8);
   h_trueE_postPandora_truePion_incident->SetBarOffset(0.5);
   h_trueE_postPandora_truePion_incident->GetXaxis()->SetTitle("KE (MeV)");

   h_unsmeared_incident->SetLineColor(kRed);
   h_unsmeared_incident->SetLineWidth(2);
   h_unsmeared_incident->SetFillColorAlpha( kRed, 0.2);
   h_unsmeared_incident->SetMarkerColor(kRed);
   h_unsmeared_incident->SetBarOffset(-0.5);

   h_recoE_selPion_incident->SetLineStyle(2);
   h_recoE_selPion_incident->SetLineWidth(3);
   h_recoE_selPion_incident->SetFillColorAlpha(kBlack, 0.1);

   h_trueE_postPandora_truePion_incident->Draw("HIST TEXT00");
   h_unsmeared_incident->Draw("HIST SAME");
   h_recoE_selPion_incident->Draw("HIST SAME");

   auto legend_inc = new TLegend(0.1,0.7,0.48,0.9);
   legend_inc->AddEntry(h_trueE_postPandora_truePion_incident,"True Pions, trueE, after PandoraReco");
   legend_inc->AddEntry(h_unsmeared_incident,"Unsmeared Incident");
   legend_inc->AddEntry(h_recoE_selPion_incident,"Selected Incident");
   legend_inc->Draw();
   c_comp_inc->Write();

   TCanvas *c_comp_initE = new TCanvas("c_comp_initE", "");

   h_trueE_postPandora_truePion_inc_initE->SetLineColor(kBlue);
   h_trueE_postPandora_truePion_inc_initE->SetLineWidth(2);
   h_trueE_postPandora_truePion_inc_initE->SetMarkerColor(kBlue);
   h_trueE_postPandora_truePion_inc_initE->SetMarkerSize(0.8);
   h_trueE_postPandora_truePion_inc_initE->SetBarOffset(0.5);
   h_trueE_postPandora_truePion_inc_initE->GetXaxis()->SetTitle("KE (MeV)");

   h_unsmeared_inc_initE->SetLineColor(kRed);
   h_unsmeared_inc_initE->SetLineWidth(2);
   h_unsmeared_inc_initE->SetFillColorAlpha( kRed, 0.2);
   h_unsmeared_inc_initE->SetMarkerColor(kRed);
   h_unsmeared_inc_initE->SetBarOffset(-0.5);

   h_recoE_selPion_inc_initE->SetLineStyle(2);
   h_recoE_selPion_inc_initE->SetLineWidth(3);
   h_recoE_selPion_inc_initE->SetFillColorAlpha(kBlack, 0.1);

   h_trueE_postPandora_truePion_inc_initE->Draw("HIST TEXT00");
   h_unsmeared_inc_initE->Draw("HIST SAME");
   h_recoE_selPion_inc_initE->Draw("HIST SAME");

   auto legend_initE = new TLegend(0.1,0.7,0.48,0.9);
   legend_initE->AddEntry(h_trueE_postPandora_truePion_inc_initE,"True Pions, trueE, after PandoraReco");
   legend_initE->AddEntry(h_unsmeared_inc_initE,"Unsmeared inc_initE");
   legend_initE->AddEntry(h_recoE_selPion_inc_initE,"Selected inc_initE");
   legend_initE->Draw();

   c_comp_initE->Write();

   TCanvas *c_comp_interE = new TCanvas("c_comp_interE", "");

   h_trueE_postPandora_truePion_inc_interE->SetLineColor(kBlue);
   h_trueE_postPandora_truePion_inc_interE->SetLineWidth(2);
   h_trueE_postPandora_truePion_inc_interE->SetMarkerColor(kBlue);
   h_trueE_postPandora_truePion_inc_interE->SetMarkerSize(0.8);
   h_trueE_postPandora_truePion_inc_interE->SetBarOffset(0.5);
   h_trueE_postPandora_truePion_inc_interE->GetXaxis()->SetTitle("KE (MeV)");

   h_unsmeared_inc_interE->SetLineColor(kRed);
   h_unsmeared_inc_interE->SetLineWidth(2);
   h_unsmeared_inc_interE->SetFillColorAlpha( kRed, 0.2);
   h_unsmeared_inc_interE->SetMarkerColor(kRed);
   h_unsmeared_inc_interE->SetBarOffset(-0.5);

   h_recoE_selPion_inc_interE->SetLineStyle(2);
   h_recoE_selPion_inc_interE->SetLineWidth(3);
   h_recoE_selPion_inc_interE->SetFillColorAlpha(kBlack, 0.1);

   h_trueE_postPandora_truePion_inc_interE->Draw("HIST TEXT00");
   h_unsmeared_inc_interE->Draw("HIST SAME");
   h_recoE_selPion_inc_interE->Draw("HIST SAME");

   auto legend_interE = new TLegend(0.1,0.7,0.48,0.9);
   legend_interE->AddEntry(h_trueE_postPandora_truePion_inc_interE,"True Pions, trueE, after PandoraReco");
   legend_interE->AddEntry(h_unsmeared_inc_interE,"Unsmeared inc_interE");
   legend_interE->AddEntry(h_recoE_selPion_inc_interE,"Selected inc_interE");
   legend_interE->Draw();

   c_comp_interE->Write();


   //------------------------------------------------------
   // Check that in truE eff and pur corrections work to get back to original hist.
   // This is in order to understand whether problems arise from the unsmearing matrices..
   //
   //------------------------------------------------------
   //    Remove the BGs and eventSelection inefficiency
   //    multiply by purity and divide by efficency
   
   

   h_trueE_effPurCorr_inc_initE->Multiply(h_trueE_selPion_inc_initE, h_trueE_pur_removeBG_inc_initE );
   h_trueE_effPurCorr_inc_initE->Divide( h_trueE_eff_eventSel_inc_initE );
   //reco ineff
   h_trueE_effPurCorr_inc_initE->Divide( h_trueE_eff_reconstruction_inc_initE );

   h_trueE_effPurCorr_inc_interE->Multiply(h_trueE_selPion_inc_interE, h_trueE_pur_removeBG_inc_interE );
   h_trueE_effPurCorr_inc_interE->Divide( h_trueE_eff_eventSel_inc_interE );
   //reco ineff
   h_trueE_effPurCorr_inc_interE->Divide( h_trueE_eff_reconstruction_inc_interE );

   //    Rebuild from initial and interacting Distribution the full incident unsmeared histogram
   //=====================================================

   build_incidentHist( h_trueE_effPurCorr_inc_initE, h_trueE_effPurCorr_inc_interE, h_trueE_effPurCorr_incident );
   h_trueE_effPurCorr_incident->Write();

   h_trueE_effPurCorr_abs_int->Multiply(h_trueE_selAbs_int, h_trueE_pur_removeBG_abs_int );
   h_trueE_effPurCorr_abs_int->Divide( h_trueE_eff_eventSel_abs_int );
   //Reco Efficiency
   h_trueE_effPurCorr_abs_int->Divide( h_trueE_eff_reconstruction_abs_int );

   h_trueE_effPurCorr_abs_int->Write();

   h_trueE_effPurCorr_totInel_int->Multiply(h_trueE_selPion_inc_interE, h_trueE_pur_removeBG_totInel_int );
   h_trueE_effPurCorr_totInel_int->Divide( h_trueE_eff_eventSel_totInel_int );
   //Reco Efficiency
   h_trueE_effPurCorr_totInel_int->Divide( h_trueE_eff_reconstruction_totInel_int );

   h_trueE_effPurCorr_totInel_int->Write();

   TCanvas *c_trueE_effPur_compare_abs_int = new TCanvas("c_trueE_effPur_compare_abs_int", "");

   h_trueE_prePandora_trueAbs_int->SetLineColor(kBlue);
   h_trueE_prePandora_trueAbs_int->SetLineWidth(2);
   h_trueE_prePandora_trueAbs_int->SetMarkerColor(kBlue);
   h_trueE_prePandora_trueAbs_int->SetBarOffset(0.5);
   h_trueE_prePandora_trueAbs_int->GetXaxis()->SetTitle("KE (MeV)");

   h_trueE_effPurCorr_abs_int->SetLineColor(kGreen);
   h_trueE_effPurCorr_abs_int->SetLineWidth(2);
   h_trueE_effPurCorr_abs_int->SetFillColorAlpha( kGreen, 0.2);
   h_trueE_effPurCorr_abs_int->SetMarkerColor(kGreen);
   h_trueE_effPurCorr_abs_int->SetBarOffset(-0.5);

   h_trueE_prePandora_trueAbs_int->Draw("HIST TEXT00");
   h_trueE_effPurCorr_abs_int->Draw("HIST SAME");

   auto leg_effPurCorr_abs_int = new TLegend(0.1,0.7,0.48,0.9);
   leg_effPurCorr_abs_int->AddEntry(h_trueE_prePandora_trueAbs_int,"True Absorption, trueE, after BeamCuts");
   leg_effPurCorr_abs_int->AddEntry(h_trueE_effPurCorr_abs_int,"Corrected Selected Abs in trueE");
   leg_effPurCorr_abs_int->Draw();
   c_trueE_effPur_compare_abs_int->Write();

   TCanvas *c_trueE_effPur_compare_totInel_int = new TCanvas("c_trueE_effPur_compare_totInel_int", "");

   h_trueE_prePandora_trueTotInel_int->SetLineColor(kBlue);
   h_trueE_prePandora_trueTotInel_int->SetLineWidth(2);
   h_trueE_prePandora_trueTotInel_int->SetMarkerColor(kBlue);
   h_trueE_prePandora_trueTotInel_int->SetBarOffset(0.5);
   h_trueE_prePandora_trueTotInel_int->GetXaxis()->SetTitle("KE (MeV)");

   h_trueE_effPurCorr_totInel_int->SetLineColor(kGreen);
   h_trueE_effPurCorr_totInel_int->SetLineWidth(2);
   h_trueE_effPurCorr_totInel_int->SetFillColorAlpha( kGreen, 0.2);
   h_trueE_effPurCorr_totInel_int->SetMarkerColor(kGreen);
   h_trueE_effPurCorr_totInel_int->SetBarOffset(-0.5);

   h_trueE_prePandora_trueTotInel_int->Draw("HIST TEXT00");
   h_trueE_effPurCorr_totInel_int->Draw("HIST SAME");

   auto leg_effPurCorr_totInel_int = new TLegend(0.1,0.7,0.48,0.9);
   leg_effPurCorr_totInel_int->AddEntry(h_trueE_prePandora_trueTotInel_int,"True Total Inelastic, trueE, after BeamCuts");
   leg_effPurCorr_totInel_int->AddEntry(h_trueE_effPurCorr_totInel_int,"Corrected Selected TotalInelastic in trueE");
   leg_effPurCorr_totInel_int->Draw();
   c_trueE_effPur_compare_totInel_int->Write();

   TCanvas *c_trueE_effPur_compare_trueE_inc = new TCanvas("c_trueE_effPur_compare_trueE_inc", "");

   h_trueE_prePandora_truePion_incident->SetLineColor(kBlue);
   h_trueE_prePandora_truePion_incident->SetLineWidth(2);
   h_trueE_prePandora_truePion_incident->SetMarkerColor(kBlue);
   h_trueE_prePandora_truePion_incident->SetMarkerSize(0.8);
   h_trueE_prePandora_truePion_incident->SetBarOffset(0.5);
   h_trueE_prePandora_truePion_incident->GetXaxis()->SetTitle("KE (MeV)");

   h_trueE_effPurCorr_incident->SetLineColor(kGreen);
   h_trueE_effPurCorr_incident->SetLineWidth(2);
   h_trueE_effPurCorr_incident->SetFillColorAlpha( kGreen, 0.2);
   h_trueE_effPurCorr_incident->SetMarkerColor(kGreen);
   h_trueE_effPurCorr_incident->SetBarOffset(-0.5);

   h_trueE_prePandora_truePion_incident->Draw("HIST TEXT00");
   h_trueE_effPurCorr_incident->Draw("HIST SAME");

   auto leg_effPurCorr_inc = new TLegend(0.1,0.7,0.48,0.9);
   leg_effPurCorr_inc->AddEntry(h_trueE_prePandora_truePion_incident,"True Pions, trueE, before Reco");
   leg_effPurCorr_inc->AddEntry(h_trueE_effPurCorr_incident,"Corrected Selected Pion Incident in trueE");
   leg_effPurCorr_inc->Draw();
   c_trueE_effPur_compare_trueE_inc->Write();

   TCanvas *c_trueE_effPur_compare_trueE_initE = new TCanvas("c_trueE_effPur_compare_trueE_initE", "");

   h_trueE_prePandora_truePion_inc_initE->SetLineColor(kBlue);
   h_trueE_prePandora_truePion_inc_initE->SetLineWidth(2);
   h_trueE_prePandora_truePion_inc_initE->SetMarkerColor(kBlue);
   h_trueE_prePandora_truePion_inc_initE->SetMarkerSize(0.8);
   h_trueE_prePandora_truePion_inc_initE->SetBarOffset(0.5);
   h_trueE_prePandora_truePion_inc_initE->GetXaxis()->SetTitle("KE (MeV)");

   h_trueE_effPurCorr_inc_initE->SetLineColor(kGreen);
   h_trueE_effPurCorr_inc_initE->SetLineWidth(2);
   h_trueE_effPurCorr_inc_initE->SetFillColorAlpha( kGreen, 0.2);
   h_trueE_effPurCorr_inc_initE->SetMarkerColor(kGreen);
   h_trueE_effPurCorr_inc_initE->SetBarOffset(-0.5);

   h_trueE_prePandora_truePion_inc_initE->Draw("HIST TEXT00");
   h_trueE_effPurCorr_inc_initE->Draw("HIST SAME");

   auto leg_effPurCorr_initE = new TLegend(0.1,0.7,0.48,0.9);
   leg_effPurCorr_initE->AddEntry(h_trueE_prePandora_truePion_inc_initE,"True Pions, trueE, before Reco");
   leg_effPurCorr_initE->AddEntry(h_trueE_effPurCorr_inc_initE,"Correcte Selected Pion initE in trueE");
   leg_effPurCorr_initE->Draw();

   c_trueE_effPur_compare_trueE_initE->Write();

   TCanvas *c_trueE_effPur_compare_trueE_interE = new TCanvas("c_trueE_effPur_compare_trueE_interE", "");

   h_trueE_prePandora_truePion_inc_interE->SetLineColor(kBlue);
   h_trueE_prePandora_truePion_inc_interE->SetLineWidth(2);
   h_trueE_prePandora_truePion_inc_interE->SetMarkerColor(kBlue);
   h_trueE_prePandora_truePion_inc_interE->SetMarkerSize(0.8);
   h_trueE_prePandora_truePion_inc_interE->SetBarOffset(0.5);
   h_trueE_prePandora_truePion_inc_interE->GetXaxis()->SetTitle("KE (MeV)");

   h_trueE_effPurCorr_inc_interE->SetLineColor(kGreen);
   h_trueE_effPurCorr_inc_interE->SetLineWidth(2);
   h_trueE_effPurCorr_inc_interE->SetFillColorAlpha( kGreen, 0.2);
   h_trueE_effPurCorr_inc_interE->SetMarkerColor(kGreen);
   h_trueE_effPurCorr_inc_interE->SetBarOffset(-0.5);

   h_trueE_prePandora_truePion_inc_interE->Draw("HIST TEXT00");
   h_trueE_effPurCorr_inc_interE->Draw("HIST SAME");

   auto leg_effPurCorr_interE = new TLegend(0.1,0.7,0.48,0.9);
   leg_effPurCorr_interE->AddEntry(h_trueE_prePandora_truePion_inc_interE,"True Pions, trueE, before Reco");
   leg_effPurCorr_interE->AddEntry(h_trueE_effPurCorr_inc_interE,"Correcte Selected Pion interE in trueE");
   leg_effPurCorr_interE->Draw();

   c_trueE_effPur_compare_trueE_interE->Write();
   

   return 0;
}

