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
//       Pandora Reco                 EventSelection Nj, and add BG
//       smearing Mi -->Nj'             purity and eff of eventSelection
//                                     in Reco bin j for int and inc sample
//
//-------------------------------------------------------
//need to go back the steps
// 1) remove BG that is not Signal, do not care about smearing, vector of purity/efficiency for the reco bin
// 2) need to apply inverse of the smearing matrix (true to reco)
//
// Validate by having the true Process Int and Inc in trueE after beamCuts
// Try to unsmear the selected Events in recoE to the above mentioned
//--------------------------------------------------------

int unsmear(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");

   //string output_name = "unsmear_halfMC_" + std::to_string((int) bin_size_int) + "MeV.root";
   string output_name = "unsmear_" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   //TrueProcess and TrueE Int and Inc Histos
   TH1D* h_trueE_trueAbs_int = new TH1D("h_trueE_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_inc_initE = new TH1D("h_trueE_truePion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_inc_interE = new TH1D("h_trueE_truePion_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_incident = new TH1D("h_trueE_truePion_incident", "", nBin_int, eEnd, eStart);

   //Selected Process and trueE Int and Inc Histos
   TH1D* h_trueE_selAbs_int = new TH1D("h_trueE_selAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_selPion_inc_initE = new TH1D("h_trueE_selPion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_selPion_inc_interE = new TH1D("h_trueE_selPion_inc_interE", "", nBin_int, eEnd, eStart);
   //TH1D* h_trueE_selPion_incident = new TH1D("h_trueE_selPion_incident", "", nBin_int, eEnd, eStart);

   //Selected Process and trueE Int and Inc Histos
   TH1D* h_trueE_trueAbs_selAbs_int = new TH1D("h_trueE_trueAbs_selAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_selPion_inc_initE = new TH1D("h_trueE_truePion_selPion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_selPion_inc_interE = new TH1D("h_trueE_truePion_selPion_inc_interE", "", nBin_int, eEnd, eStart);
   //TH1D* h_trueE_truePion_selPion_incident = new TH1D("h_trueE_truePion_selPion_incident", "", nBin_int, eEnd, eStart);

   //Selected Process and RecoE Int and Inc Histos
   TH1D* h_recoE_selAbs_int = new TH1D("h_recoE_selAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_inc_initE = new TH1D("h_recoE_selPion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_inc_interE = new TH1D("h_recoE_selPion_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_incident = new TH1D("h_recoE_selPion_incident", "", nBin_int, eEnd, eStart);

   TH1D* h_recoE_incidentPion_truePion_inc_initE = new TH1D("h_recoE_incidentPion_truePion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_beamCut_truePion_inc_initE = new TH1D("h_recoE_beamCut_truePion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_incidentPion_truePion_inc_interE = new TH1D("h_recoE_incidentPion_truePion_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_beamCut_truePion_inc_interE = new TH1D("h_recoE_beamCut_truePion_inc_interE", "", nBin_int, eEnd, eStart);

   TH1D* h_recoE_selAbs_trueAbs_int = new TH1D("h_recoE_selAbs_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_beamCut_trueAbs_int = new TH1D("h_recoE_beamCut_trueAbs_int", "", nBin_int, eEnd, eStart);

   //Unsmeared Histo, compare to h_trueE_trueProc
   TH1D* h_help_unsmear_int = new TH1D("h_help_unsmear_int", "", nBin_int, eEnd, eStart);
   TH1D* h_help_unsmear_inc_initE = new TH1D("h_help_unsmear_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_help_unsmear_inc_interE = new TH1D("h_help_unsmear_inc_interE", "", nBin_int, eEnd, eStart);

   TH1D* h_unsmeared_int = new TH1D("h_unsmeared_int", "", nBin_int, eEnd, eStart);
   TH1D* h_unsmeared_inc_initE = new TH1D("h_unsmeared_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_unsmeared_inc_interE = new TH1D("h_unsmeared_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_unsmeared_incident = new TH1D("h_unsmeared_incident", "", nBin_int, eEnd, eStart);
   //From evSel --> back to Reco MC Nj --> Nj'
   TH1D* h_pur_removeBG_int = new TH1D("h_pur_removeBG_int", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_int = new TH1D("h_eff_eventSel_int", "", nBin_int, eEnd, eStart);

   TH1D* h_pur_removeBG_inc_initE = new TH1D("h_pur_removeBG_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_inc_initE = new TH1D("h_eff_eventSel_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_pur_removeBG_inc_interE = new TH1D("h_pur_removeBG_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_inc_interE = new TH1D("h_eff_eventSel_inc_interE", "", nBin_int, eEnd, eStart);

   TH1D* h_trueE_pur_removeBG_int = new TH1D("h_trueE_pur_removeBG_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_eff_eventSel_int = new TH1D("h_trueE_eff_eventSel_int", "", nBin_int, eEnd, eStart);

   TH1D* h_trueE_pur_removeBG_inc_initE = new TH1D("h_trueE_pur_removeBG_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_eff_eventSel_inc_initE = new TH1D("h_trueE_eff_eventSel_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_pur_removeBG_inc_interE = new TH1D("h_trueE_pur_removeBG_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_eff_eventSel_inc_interE = new TH1D("h_trueE_eff_eventSel_inc_interE", "", nBin_int, eEnd, eStart);

   TH1D* h_trueE_effPurCorr_int = new TH1D("h_trueE_effPurCorr_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_effPurCorr_inc_initE = new TH1D("h_trueE_effPurCorr_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_effPurCorr_inc_interE = new TH1D("h_trueE_effPurCorr_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_effPurCorr_incident = new TH1D("h_trueE_effPurCorr_incident", "", nBin_int, eEnd, eStart);
   //Remove events where reco_interE == reco_incidentE bin && true_interE == true_initE bin
   //
   //events where true initE == true inter E but not in reco are considered as BG
   //events where reco initE == reco inter E but not the case in true, are considered as inefficiency
   //
   auto frame = inputFrame
      .Define("true_reco_initE_eq_interE", keepEv_notSame_eBin ,
         {"true_firstEntryIncident", "true_KEint_fromEndP", "reco_firstEntryIncident", "reco_interactingKE"})
   
   .Define("reco_initE_eq_interE", [](double reco_initE, double reco_interE){
            if( initE_sameBin_interE(reco_initE, reco_interE) ) return true;
            else return false;
         }
         ,{"reco_firstEntryIncident", "reco_interactingKE"})

   .Define("true_initE_eq_interE", [](double true_initE, double true_interE){
            if( initE_sameBin_interE(true_initE, true_interE) ) return true;
            else return false;
         }
         ,{"true_firstEntryIncident", "true_KEint_fromEndP"});
   //.Filter("true_reco_initE_eq_interE"); //false if trueEinit and recoEinit == their interacting
   //.Filter("!reco_initE_eq_interE")
   //.Filter("!true_initE_eq_interE");



   //Build the True Process and TrueE Int and Inc Histograms that we need to compare unsmeared things to
   //
   //all available after beamCuts
   auto eventSel_post_beamCut = frame.Filter("primary_isBeamType && passBeamCut && passBeamCutBI");
   //selected incident Pions & selected absorption
   //make sure they don't have reco initE == reco InterE bc in Selection they will be rejected as I can't fill my histos with such evnts
   auto eventSel_incidentPion = frame.Filter("selected_incidentPion && !reco_initE_eq_interE");
   auto eventSel_abs = frame.Filter("selected_abs && !reco_initE_eq_interE");

   //----------------------------------
   //Fill the Histos with true Process and trueEbin Incident and interacting.
   //
   //We should get the amount of events that passed the beamCuts i.e. was available initially for analysis
   //!!!! This does not take into account the reconstruction efficiency
   // Maybe will have to go a step further
   //
   //this is what should be achieved after unsmearing the selected&reco histos
   //
   //------------------------------------------------------------
   //Incident Histo InitE
   eventSel_post_beamCut
      .Filter("!true_initE_eq_interE")
      .Filter("true_primPionInel") //should also take into account pions that decay

      .Foreach( [h_trueE_truePion_inc_initE, h_trueE_truePion_inc_interE] (double true_firstEntryIncident, double true_beam_interactingEnergy) { 
            h_trueE_truePion_inc_initE->Fill( true_firstEntryIncident ); 
            h_trueE_truePion_inc_interE->Fill( true_beam_interactingEnergy ); 

            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP"});


   //Build trueIncident for comparison
   TH1* h_cum_trueE_initE = h_trueE_truePion_inc_initE->GetCumulative(false);
   TH1* h_cum_trueE_interE = h_trueE_truePion_inc_interE->GetCumulative(false);
   //h_cum_trueE_initE->Write();
   //h_cum_trueE_interE->Write();
   //
   //FILLING OF INCIDENT
   //--> if Pion is "born" in bin 10 and interacts in bin 4 then incident is filled from bin 9 to 4
   //cumulative sum of initialE: the difference between higher bin and the lower next shows how many pions were born 
   //cumulative sum of interE: the difference between the higher bin and lower shows how many pions died
   //incident histo is built such that it takes the difference of what was born in previous bin compared to what died in previous bin
   //i.e. incident histo cannot fill the first energy bin of KE 1200MeV
   for(int i=nBin_int-1; i >= 1; i--){

      double diff = h_cum_trueE_initE->GetBinContent(i+1) - h_cum_trueE_interE->GetBinContent(i+1);

      h_trueE_truePion_incident->SetBinContent( i, diff);
   };

   h_trueE_truePion_incident->Write();

   //Interacting Histo
   eventSel_post_beamCut
      .Filter("!true_initE_eq_interE")
      .Filter("true_absSignal") //should also take into account pions that decay

      .Foreach( [h_trueE_trueAbs_int] (double true_firstEntryIncident, double true_beam_interactingEnergy){

            h_trueE_trueAbs_int->Fill(true_beam_interactingEnergy);
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP"});


   h_trueE_trueAbs_int->Write();
   h_trueE_truePion_inc_initE->Write();
   h_trueE_truePion_inc_interE->Write();
   //------------------------------------------------------------
   //
   //Create the selected and Reco histos that need to be unsmeared back
   //
   //
   //------------------------------------------------------------
   //
   //for this already removed recoInitE == recoInterE
   eventSel_incidentPion
      .Foreach( [h_recoE_selPion_inc_initE, h_recoE_selPion_inc_interE] (double reco_firstEntryIncident, double reco_beam_interactingEnergy ) {

            h_recoE_selPion_inc_initE->Fill( reco_firstEntryIncident ); 
            h_recoE_selPion_inc_interE->Fill( reco_beam_interactingEnergy );     
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});

   //Build recoIncident for comparison
   TH1* h_cum_recoE_initE = h_recoE_selPion_inc_initE->GetCumulative(false);
   TH1* h_cum_recoE_interE = h_recoE_selPion_inc_interE->GetCumulative(false);
   //h_cum_recoE_initE->Write();
   //h_cum_recoE_interE->Write();
   //h_recoE_selPion_incident->SetBinContent(nBin_int, h_cum_recoE_initE->GetBinContent(nBin_int));
   for(int i=nBin_int-1; i >= 1; i--){
      double diff = h_cum_recoE_initE->GetBinContent(i+1) - h_cum_recoE_interE->GetBinContent(i+1);

      h_recoE_selPion_incident->SetBinContent( i, diff);
   };

   h_recoE_selPion_incident->Write();

   //for this already removed recoInitE == recoInterE
   eventSel_abs
      .Foreach( [h_recoE_selAbs_int] ( double reco_firstEntryIncident, double reco_beam_interactingEnergy ){

            h_recoE_selAbs_int->Fill(reco_beam_interactingEnergy);
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
   eventSel_post_beamCut
      //Events where trueInitE == trueInterE should not be available for this selection they are not a signal
      .Filter("true_primPionInel && !true_initE_eq_interE")
      .Foreach( [h_recoE_beamCut_truePion_inc_initE, h_recoE_beamCut_truePion_inc_interE] (double reco_firstEntryIncident, double reco_beam_interactingEnergy ) {

            h_recoE_beamCut_truePion_inc_initE->Fill( reco_firstEntryIncident ); 
            h_recoE_beamCut_truePion_inc_interE->Fill( reco_beam_interactingEnergy ); 
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});


   //True Pions in selected incident
   eventSel_incidentPion
      //Events where trueInitE == trueInterE should not be available for this selection they are not a signal
      .Filter("true_primPionInel && !true_initE_eq_interE")
      .Foreach( [h_recoE_incidentPion_truePion_inc_initE, h_recoE_incidentPion_truePion_inc_interE] (double reco_firstEntryIncident, double reco_beam_interactingEnergy ) {

            h_recoE_incidentPion_truePion_inc_initE->Fill( reco_firstEntryIncident ); 
            h_recoE_incidentPion_truePion_inc_interE->Fill( reco_beam_interactingEnergy ); 
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});

   //True Abs available after beamCuts
   eventSel_post_beamCut
      //Events where trueInitE == trueInterE should not be available for this selection they are not a signal
      .Filter("true_absSignal && !true_initE_eq_interE")
      .Foreach( [h_recoE_beamCut_trueAbs_int] ( double reco_firstEntryIncident, double reco_beam_interactingEnergy ){

            h_recoE_beamCut_trueAbs_int->Fill(reco_beam_interactingEnergy);
            }            
            ,{"reco_firstEntryIncident","reco_interactingKE"});
   //True Abs available after Selection
   eventSel_abs
      //Events where trueInitE == trueInterE should not be available for this selection they are not a signal
      .Filter("true_absSignal && !true_initE_eq_interE")
      .Filter("true_absSignal")
      .Foreach( [h_recoE_selAbs_trueAbs_int] ( double reco_firstEntryIncident, double reco_beam_interactingEnergy ){

            h_recoE_selAbs_trueAbs_int->Fill(reco_beam_interactingEnergy);
            }            
            ,{"reco_firstEntryIncident","reco_interactingKE"});

   //True Pions trueE in Pion Selection
   //
   eventSel_incidentPion
      .Filter("true_primPionInel && !true_initE_eq_interE")
      .Foreach( [h_trueE_truePion_selPion_inc_initE, h_trueE_truePion_selPion_inc_interE] (double true_firstEntryIncident, double true_beam_interactingEnergy ) {

            h_trueE_truePion_selPion_inc_initE->Fill( true_firstEntryIncident ); 
            h_trueE_truePion_selPion_inc_interE->Fill( true_beam_interactingEnergy ); 
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP"});

   //Selected Pions trueE in Pion Selection
   //
   eventSel_incidentPion
      .Filter("!true_initE_eq_interE")
      .Foreach( [h_trueE_selPion_inc_initE, h_trueE_selPion_inc_interE] (double true_firstEntryIncident, double true_beam_interactingEnergy ) {

            h_trueE_selPion_inc_initE->Fill( true_firstEntryIncident ); 
            h_trueE_selPion_inc_interE->Fill( true_beam_interactingEnergy ); 
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP"});

   //True Abs trueE in Pion Selection
   //
   eventSel_abs
      .Filter("true_absSignal && !true_initE_eq_interE")      
      .Foreach( [h_trueE_trueAbs_selAbs_int] ( double true_beam_interactingEnergy ) {

            h_trueE_trueAbs_selAbs_int->Fill( true_beam_interactingEnergy ); 
            }
            ,{"true_KEint_fromEndP"});

   //Selected Abs trueE in Abs Selection
   //
   eventSel_abs
      .Filter("!true_initE_eq_interE")
      .Foreach( [h_trueE_selAbs_int] ( double true_beam_interactingEnergy ) {

            h_trueE_selAbs_int->Fill( true_beam_interactingEnergy ); 
            }
            ,{"true_KEint_fromEndP"});



   //------------------------------------------------------------
   //
   //  Compute eff and Purities for first unsmearing Step
   //------------------------------------------------------------
   //
   h_pur_removeBG_inc_initE->Divide( h_recoE_incidentPion_truePion_inc_initE, h_recoE_selPion_inc_initE );
   h_eff_eventSel_inc_initE->Divide( h_recoE_incidentPion_truePion_inc_initE, h_recoE_beamCut_truePion_inc_initE );

   h_pur_removeBG_inc_interE->Divide( h_recoE_incidentPion_truePion_inc_interE, h_recoE_selPion_inc_interE );
   h_eff_eventSel_inc_interE->Divide( h_recoE_incidentPion_truePion_inc_interE, h_recoE_beamCut_truePion_inc_interE );

   h_pur_removeBG_int->Divide( h_recoE_selAbs_trueAbs_int, h_recoE_selAbs_int );
   h_eff_eventSel_int->Divide( h_recoE_selAbs_trueAbs_int, h_recoE_beamCut_trueAbs_int );

   h_pur_removeBG_inc_initE->Write();
   h_eff_eventSel_inc_initE->Write();
   h_pur_removeBG_inc_interE->Write();
   h_eff_eventSel_inc_interE->Write();

   h_pur_removeBG_int->Write();
   h_eff_eventSel_int->Write();

   h_trueE_pur_removeBG_inc_initE->Divide( h_trueE_truePion_selPion_inc_initE, h_trueE_selPion_inc_initE);
   h_trueE_eff_eventSel_inc_initE->Divide( h_trueE_truePion_selPion_inc_initE, h_trueE_truePion_inc_initE);

   h_trueE_pur_removeBG_inc_interE->Divide( h_trueE_truePion_selPion_inc_interE, h_trueE_selPion_inc_interE);
   h_trueE_eff_eventSel_inc_interE->Divide( h_trueE_truePion_selPion_inc_interE, h_trueE_truePion_inc_interE);

   h_trueE_pur_removeBG_int->Divide( h_trueE_trueAbs_selAbs_int, h_trueE_selAbs_int);
   h_trueE_eff_eventSel_int->Divide( h_trueE_trueAbs_selAbs_int, h_trueE_trueAbs_int);

   h_trueE_pur_removeBG_inc_initE->Write();
   h_trueE_eff_eventSel_inc_initE->Write();
   h_trueE_pur_removeBG_inc_interE->Write();
   h_trueE_eff_eventSel_inc_interE->Write();

   h_trueE_pur_removeBG_int->Write();
   h_trueE_eff_eventSel_int->Write();

   //------------------------------------------------------------
   //
   //  Now build Inverse Smearing Matrix to translate back to recoE j --> trueE i
   //
   //  matrix Tij*Nj=Mi
   //  Tij has true on y and reco on x axis, normalised on reco
   //------------------------------------------------------------
   //

   TH2D* h2_smearing_interacting = new TH2D("h2_smearing_interacting", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_smearing_incident_initE = new TH2D("h2_smearing_incident_initE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_smearing_incident_interE = new TH2D("h2_smearing_incident_interE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   h2_smearing_interacting->SetTitle("Smearing Matrix for trueAbs in Sample after BeamCuts; reco Energy [MeV]; true Energy [MeV]");
   h2_smearing_incident_initE->SetTitle("Smearing Matrix for truePi initialE in Sample after BeamCuts; reco Energy [MeV]; true Energy [MeV]");
   h2_smearing_incident_interE->SetTitle("Smearing Matrix for truePi interE in Sample after BeamCuts; reco Energy [MeV]; true Energy [MeV]");

   TH2D* h2_inverse_smearing_interacting = new TH2D("h2_inverse_smearing_interacting", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_inverse_smearing_incident_initE = new TH2D("h2_inverse_smearing_incident_initE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_inverse_smearing_incident_interE = new TH2D("h2_inverse_smearing_incident_interE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   h2_inverse_smearing_interacting->SetTitle("Inverse smearing Matrix for trueAbs in Sample after BeamCuts; true Energy [MeV]; reco Energy [MeV]");
   h2_inverse_smearing_incident_initE->SetTitle("Inverse smearing Matrix for truePi initialE in Sample after BeamCuts; true Energy [MeV]; reco Energy [MeV]");
   h2_inverse_smearing_incident_interE->SetTitle("Inverse smearing Matrix for truePi interE in Sample after BeamCuts; true Energy [MeV]; reco Energy [MeV]");

   //Fill smearing matrix with entry 1 on the diagonal to avoid them being "uninvertible"
   for(int i = 1; i <= nBin_int; i++){
      h2_smearing_incident_initE->SetBinContent(i,i,1);
      h2_smearing_incident_interE->SetBinContent(i,i,1);
      h2_smearing_interacting->SetBinContent(i,i,1);
   };
   //========================================================
   //Build the smearing Incident Histogram initE
   //IGNORE events that have true initE-bin == interE-bin
   //---------

   eventSel_post_beamCut 
      //.Filter("true_primPionInel && !reco_initE_eq_interE")
      .Filter("true_primPionInel && !true_initE_eq_interE")
      //.Filter("true_primPionInel && true_reco_initE_eq_interE")
      .Foreach( [h2_smearing_incident_initE, h2_smearing_incident_interE] (double true_initE, double true_interE, double reco_initE, double reco_interE ){

            h2_smearing_incident_initE->Fill( reco_initE, true_initE );
            h2_smearing_incident_interE->Fill( reco_interE, true_interE ); 
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP", "reco_firstEntryIncident", "reco_interactingKE"}); 


   //NORMALISATION of initE and interE
   ////
   for(int i = 1; i <= nBin_int; i++){

      //sum for same x (recoBin
      int sum_true_initE = h2_smearing_incident_initE->Integral( 1, nBin_int, i, i ); //sum up all recoE signals for one trueE signal
      int sum_true_interE = h2_smearing_incident_interE->Integral( 1, nBin_int, i, i ); //sum up all recoE signals for one trueE signal

      if(sum_true_initE !=0){
         for(int j = 1; j <= nBin_int; j++){

            h2_smearing_incident_initE->SetBinContent( j , i, h2_smearing_incident_initE->GetBinContent(j,i) / sum_true_initE);

         };
      };

      if(sum_true_interE !=0){
         for(int j = 1; j <= nBin_int; j++){

            h2_smearing_incident_interE->SetBinContent( j , i, h2_smearing_incident_interE->GetBinContent(j,i) / sum_true_interE);

         };
      };
   };

   h2_smearing_incident_initE->Sumw2(0);
   h2_smearing_incident_initE->Write();

   h2_smearing_incident_interE->Sumw2(0);
   h2_smearing_incident_interE->Write();

   //TMatrix initE
   //have to exclude underflow and overflowbin bc otherwise matrix is singular
   TMatrixD matrix_smearing_incident_initE_pre(nBin_int + 2, nBin_int + 2, h2_smearing_incident_initE->GetArray(), "D");
   TMatrixD matrix_smearing_incident_initE = matrix_smearing_incident_initE_pre.GetSub(1, nBin_int, 1, nBin_int); //avoid underflow and overflowbin

   //matrix_smearing_incident_initE.Print();

   Double_t det1;
   TMatrixD matrix_inverse_smearing_incident_initE = matrix_smearing_incident_initE;
   matrix_inverse_smearing_incident_initE.Invert(&det1);

   TMatrixD U1(matrix_inverse_smearing_incident_initE, TMatrixD::kMult, matrix_smearing_incident_initE);
   TMatrixDDiag diag1(U1); diag1 = 0.0;
   const Double_t U1_max_offdiag = (U1.Abs()).Max();
   std::cout << "  Maximum off-diagonal = " << U1_max_offdiag << std::endl;
   std::cout << "  Determinant          = " << det1 << std::endl;   

   //put inverse into TH2
   for (int i = 1; i <= nBin_int; i++){
      for (int j= 1; j <= nBin_int; j++){
         h2_inverse_smearing_incident_initE->SetBinContent(j, i, matrix_inverse_smearing_incident_initE(i-1,j-1)); //vector indices style for matrix 
      };
   };

   h2_inverse_smearing_incident_initE->Write();

   //Check that matrix and inverse are unity
   TMatrixD unity_smearing_incident_initE = matrix_inverse_smearing_incident_initE * matrix_smearing_incident_initE ;
   TH2D* h2_prod_smearing_incident_initE = new TH2D("h2_prod_smearing_incident_initE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   for (int i = 1; i <= nBin_int; i++){
      for (int j= 1; j <= nBin_int; j++){
         h2_prod_smearing_incident_initE->SetBinContent(j, i, unity_smearing_incident_initE(i-1,j-1)); 
      }
   }

   h2_prod_smearing_incident_initE->Write();



   //TMatrix interE
   //have to exclude underflow and overflowbin bc otherwise matrix is singular
   TMatrixD matrix_smearing_incident_interE_pre(nBin_int + 2, nBin_int + 2, h2_smearing_incident_interE->GetArray(), "D");
   TMatrixD matrix_smearing_incident_interE = matrix_smearing_incident_interE_pre.GetSub(1, nBin_int, 1, nBin_int); //avoid underflow and overflowbin

   //matrix_smearing_incident_interE.Print();

   Double_t det3;
   TMatrixD matrix_inverse_smearing_incident_interE = matrix_smearing_incident_interE;
   matrix_inverse_smearing_incident_interE.Invert(&det3);

   TMatrixD U3(matrix_inverse_smearing_incident_interE, TMatrixD::kMult, matrix_smearing_incident_interE);
   TMatrixDDiag diag3(U3); diag3 = 0.0;
   const Double_t U3_max_offdiag = (U3.Abs()).Max();
   std::cout << "  Maximum off-diagonal = " << U3_max_offdiag << std::endl;
   std::cout << "  Determinant          = " << det3 << std::endl;   

   //put inverse into TH2
   for (int i = 1; i <= nBin_int; i++){
      for (int j= 1; j <= nBin_int; j++){
         h2_inverse_smearing_incident_interE->SetBinContent(j, i, matrix_inverse_smearing_incident_interE(i-1,j-1)); //vector indices style for matrix 
      };
   };

   h2_inverse_smearing_incident_interE->Write();

   //Check that matrix and inverse are unity
   TMatrixD unity_smearing_incident_interE = matrix_inverse_smearing_incident_interE * matrix_smearing_incident_interE ;
   TH2D* h2_prod_smearing_incident_interE = new TH2D("h2_prod_smearing_incident_interE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   for (int i = 1; i <= nBin_int; i++){
      for (int j= 1; j <= nBin_int; j++){
         h2_prod_smearing_incident_interE->SetBinContent(j, i, unity_smearing_incident_interE(i-1,j-1)); 
      }
   }

   h2_prod_smearing_incident_interE->Write();


   //========================================================
   //Build the smearing Interacting Histogram
   //---------
   //eventSel_abs
   //Ignore trueE init == trueE int as they are NOT a signal
   eventSel_post_beamCut
      .Filter("true_absSignal && !true_initE_eq_interE")
      //.Filter("true_absSignal && true_reco_initE_eq_interE")
      .Foreach( [h2_smearing_interacting] (double true_initE, double true_interE, double reco_initE, double reco_interE){

            h2_smearing_interacting->Fill( reco_interE, true_interE );             
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP", "reco_firstEntryIncident", "reco_interactingKE"}); 

   //NORMALISATION Normalise to recoE column
   //Go through true columns and normalise the entries 
   for(int i = 1; i <= nBin_int; i++){

      int sum_true_i = h2_smearing_interacting->Integral( 1 , nBin_int, i, i ) ; //sum up all recoE signals for one trueE signal

      if(sum_true_i !=0){
         for(int j = 1; j <= nBin_int; j++){

            h2_smearing_interacting->SetBinContent( j, i, h2_smearing_interacting->GetBinContent(j,i) / sum_true_i);

         };
      }
   };

   h2_smearing_interacting->Sumw2(0);
   h2_smearing_interacting->Write();


   //have to exclude underflow and overflowbin bc otherwise matrix is singular
   TMatrixD matrix_smearing_interacting_pre(nBin_int + 2, nBin_int + 2, h2_smearing_interacting->GetArray(), "D");
   TMatrixD matrix_smearing_interacting = matrix_smearing_interacting_pre.GetSub(1, nBin_int, 1, nBin_int); //avoid underflow and overflowbin

   //matrix_smearing_interacting.Print();

   Double_t det2;
   TMatrixD matrix_inverse_smearing_interacting = matrix_smearing_interacting;
   matrix_inverse_smearing_interacting.Invert(&det2);

   TMatrixD U2(matrix_inverse_smearing_interacting, TMatrixD::kMult, matrix_smearing_interacting);
   TMatrixDDiag diag2(U2); diag2 = 0.0;
   const Double_t U2_max_offdiag = (U2.Abs()).Max();
   std::cout << "  Maximum off-diagonal = " << U2_max_offdiag << std::endl;
   std::cout << "  Determinant          = " << det2 << std::endl;   

   //put inverse into TH2
   for (int i = 1; i <= nBin_int; i++){
      for (int j= 1; j <= nBin_int; j++){
         h2_inverse_smearing_interacting->SetBinContent(j, i, matrix_inverse_smearing_interacting(i-1,j-1)); //vector indices style for matrix 
      };
   };

   h2_inverse_smearing_interacting->Write();

   //Check that matrix and inverse are unity
   TMatrixD unity_smearing_interacting = matrix_inverse_smearing_interacting * matrix_smearing_interacting ;
   TH2D* h2_prod_smearing_interacting = new TH2D("h2_prod_smearing_interacting", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   for (int i = 1; i <= nBin_int; i++){
      for (int j= 1; j <= nBin_int; j++){
         h2_prod_smearing_interacting->SetBinContent(j, i, unity_smearing_interacting(i-1,j-1)); 
      }
   }

   h2_prod_smearing_interacting->Write();

   //=====================================================
   //------------------------------------------------------
   //Start Unsmearing Incident Init E
   //------------------------------------------------------
   //
   //    Remove the BGs and eventSelection inefficiency
   //    multiply by purity and divide by efficency
   h_help_unsmear_inc_initE->Multiply(h_recoE_selPion_inc_initE, h_pur_removeBG_inc_initE );
   h_help_unsmear_inc_initE->Divide( h_eff_eventSel_inc_initE );

   //   Undo the smearing by multiplying Tij*Nj, so every row of true percentage with every entry of help_unsmeared
   //   Row's are i(y), columns are j(x) 
   //
   //Loop through smearing matrix Incident
   for(int i = 1; i <= nBin_int; i++){
      double help_sum = 0;
      //rows
      for(int j = 1; j <= nBin_int; j++){
         help_sum += h_help_unsmear_inc_initE->GetBinContent(j)*h2_inverse_smearing_incident_initE->GetBinContent( i, j);
      }; 
      h_unsmeared_inc_initE->SetBinContent( i , help_sum);
   };

   h_unsmeared_inc_initE->Write();

   //------------------------------------------------------
   //Start Unsmearing Incident inter E
   //------------------------------------------------------
   //
   //    Remove the BGs and eventSelection inefficiency
   //    multiply by purity and divide by efficency
   h_help_unsmear_inc_interE->Multiply(h_recoE_selPion_inc_interE, h_pur_removeBG_inc_interE );
   h_help_unsmear_inc_interE->Divide( h_eff_eventSel_inc_interE );

   //   Undo the smearing by multiplying Tij*Nj, so every row of true percentage with every entry of help_unsmeared
   //   Row's are i(y), columns are j(x) 
   //
   //Loop through smearing matrix Interacting
   for(int i = 1; i <= nBin_int; i++){
      double help_sum = 0;
      //rows
      for(int j = 1; j <= nBin_int; j++){
         help_sum += h_help_unsmear_inc_interE->GetBinContent(j)*h2_inverse_smearing_incident_interE->GetBinContent( i, j);
      }; 
      h_unsmeared_inc_interE->SetBinContent( i , help_sum);
   };

   h_unsmeared_inc_interE->Write();

   //=====================================================
   //    Rebuild from initial and interacting Distribution the full incident unsmeared histogram
   //=====================================================

   TH1* h_cum_initE = h_unsmeared_inc_initE->GetCumulative(false);
   TH1* h_cum_interE = h_unsmeared_inc_interE->GetCumulative(false);
   //h_cum_initE->Write();
   //h_cum_interE->Write();

   //h_unsmeared_incident->SetBinContent(nBin_int, h_cum_initE->GetBinContent(nBin_int));
   for(int i=nBin_int-1; i >= 1; i--){
      double diff = h_cum_initE->GetBinContent(i+1) - h_cum_interE->GetBinContent(i+1);
      h_unsmeared_incident->SetBinContent( i, diff);
   };

   h_unsmeared_incident->Write();

   //------------------------------------------------------
   //Start Unsmearing Interacting
   //------------------------------------------------------
   //
   //    Remove the BGs and eventSelection inefficiency
   //    multiply by purity and divide by efficency
   h_help_unsmear_int->Multiply(h_recoE_selAbs_int, h_pur_removeBG_int );
   h_help_unsmear_int->Divide( h_eff_eventSel_int );

   //   Undo the smearing by multiplying Tij*Nj, so every row of true percentage with every entry of help_unsmeared
   //   Row's are i(y), columns are j(x) 
   //
   //Loop through smearing matrix Interacting
   for(int i = 1; i <= nBin_int; i++){
      double help_sum = 0;
      //rows
      for(int j = 1; j <= nBin_int; j++){
         help_sum += h_help_unsmear_int->GetBinContent(j)*h2_inverse_smearing_interacting->GetBinContent( i, j);
      }; 
      h_unsmeared_int->SetBinContent( i , help_sum);
   };

   h_unsmeared_int->Write();

   //------------------------------------------------------
   //          CANVAS
   //------------------------------------------------------
   //

   TCanvas *c_comp_int = new TCanvas("c_comp_int", "");

   h_trueE_trueAbs_int->SetLineColor(kBlue);
   h_trueE_trueAbs_int->SetLineWidth(2);
   h_trueE_trueAbs_int->SetMarkerColor(kBlue);
   h_trueE_trueAbs_int->SetBarOffset(0.5);
   h_trueE_trueAbs_int->GetXaxis()->SetTitle("KE (MeV)");

   h_unsmeared_int->SetLineColor(kRed);
   h_unsmeared_int->SetLineWidth(2);
   h_unsmeared_int->SetFillColorAlpha( kRed, 0.2);
   h_unsmeared_int->SetMarkerColor(kRed);
   h_unsmeared_int->SetBarOffset(-0.5);

   h_recoE_selAbs_int->SetLineStyle(2);
   h_recoE_selAbs_int->SetLineWidth(3);
   h_recoE_selAbs_int->SetFillColorAlpha(kBlack, 0.1);

   h_trueE_trueAbs_int->Draw("HIST TEXT00");
   h_unsmeared_int->Draw("HIST SAME");
   h_recoE_selAbs_int->Draw("HIST SAME");

   auto legend_int = new TLegend(0.1,0.7,0.48,0.9);
   legend_int->AddEntry(h_trueE_trueAbs_int,"True Absorption, trueE, after BeamCuts");
   legend_int->AddEntry(h_unsmeared_int,"Unsmeared Interacting");
   legend_int->AddEntry(h_recoE_selAbs_int,"Selected Interacting");
   legend_int->Draw();
   c_comp_int->Write();

   TCanvas *c_comp_inc = new TCanvas("c_comp_inc", "");

   h_trueE_truePion_incident->SetLineColor(kBlue);
   h_trueE_truePion_incident->SetLineWidth(2);
   h_trueE_truePion_incident->SetMarkerColor(kBlue);
   h_trueE_truePion_incident->SetMarkerSize(0.8);
   h_trueE_truePion_incident->SetBarOffset(0.5);
   h_trueE_truePion_incident->GetXaxis()->SetTitle("KE (MeV)");

   h_unsmeared_incident->SetLineColor(kRed);
   h_unsmeared_incident->SetLineWidth(2);
   h_unsmeared_incident->SetFillColorAlpha( kRed, 0.2);
   h_unsmeared_incident->SetMarkerColor(kRed);
   h_unsmeared_incident->SetBarOffset(-0.5);

   h_recoE_selPion_incident->SetLineStyle(2);
   h_recoE_selPion_incident->SetLineWidth(3);
   h_recoE_selPion_incident->SetFillColorAlpha(kBlack, 0.1);

   h_trueE_truePion_incident->Draw("HIST TEXT00");
   h_unsmeared_incident->Draw("HIST SAME");
   h_recoE_selPion_incident->Draw("HIST SAME");

   auto legend_inc = new TLegend(0.1,0.7,0.48,0.9);
   legend_inc->AddEntry(h_trueE_truePion_incident,"True Pions, trueE, after BeamCuts");
   legend_inc->AddEntry(h_unsmeared_incident,"Unsmeared Incident");
   legend_inc->AddEntry(h_recoE_selPion_incident,"Selected Incident");
   legend_inc->Draw();
   c_comp_inc->Write();

   TCanvas *c_comp_initE = new TCanvas("c_comp_initE", "");

   h_trueE_truePion_inc_initE->SetLineColor(kBlue);
   h_trueE_truePion_inc_initE->SetLineWidth(2);
   h_trueE_truePion_inc_initE->SetMarkerColor(kBlue);
   h_trueE_truePion_inc_initE->SetMarkerSize(0.8);
   h_trueE_truePion_inc_initE->SetBarOffset(0.5);
   h_trueE_truePion_inc_initE->GetXaxis()->SetTitle("KE (MeV)");

   h_unsmeared_inc_initE->SetLineColor(kRed);
   h_unsmeared_inc_initE->SetLineWidth(2);
   h_unsmeared_inc_initE->SetFillColorAlpha( kRed, 0.2);
   h_unsmeared_inc_initE->SetMarkerColor(kRed);
   h_unsmeared_inc_initE->SetBarOffset(-0.5);

   h_recoE_selPion_inc_initE->SetLineStyle(2);
   h_recoE_selPion_inc_initE->SetLineWidth(3);
   h_recoE_selPion_inc_initE->SetFillColorAlpha(kBlack, 0.1);

   h_trueE_truePion_inc_initE->Draw("HIST TEXT00");
   h_unsmeared_inc_initE->Draw("HIST SAME");
   h_recoE_selPion_inc_initE->Draw("HIST SAME");

   auto legend_initE = new TLegend(0.1,0.7,0.48,0.9);
   legend_initE->AddEntry(h_trueE_truePion_inc_initE,"True Pions, trueE, after BeamCuts");
   legend_initE->AddEntry(h_unsmeared_inc_initE,"Unsmeared inc_initE");
   legend_initE->AddEntry(h_recoE_selPion_inc_initE,"Selected inc_initE");
   legend_initE->Draw();

   c_comp_initE->Write();

   TCanvas *c_comp_interE = new TCanvas("c_comp_interE", "");

   h_trueE_truePion_inc_interE->SetLineColor(kBlue);
   h_trueE_truePion_inc_interE->SetLineWidth(2);
   h_trueE_truePion_inc_interE->SetMarkerColor(kBlue);
   h_trueE_truePion_inc_interE->SetMarkerSize(0.8);
   h_trueE_truePion_inc_interE->SetBarOffset(0.5);
   h_trueE_truePion_inc_interE->GetXaxis()->SetTitle("KE (MeV)");

   h_unsmeared_inc_interE->SetLineColor(kRed);
   h_unsmeared_inc_interE->SetLineWidth(2);
   h_unsmeared_inc_interE->SetFillColorAlpha( kRed, 0.2);
   h_unsmeared_inc_interE->SetMarkerColor(kRed);
   h_unsmeared_inc_interE->SetBarOffset(-0.5);

   h_recoE_selPion_inc_interE->SetLineStyle(2);
   h_recoE_selPion_inc_interE->SetLineWidth(3);
   h_recoE_selPion_inc_interE->SetFillColorAlpha(kBlack, 0.1);

   h_trueE_truePion_inc_interE->Draw("HIST TEXT00");
   h_unsmeared_inc_interE->Draw("HIST SAME");
   h_recoE_selPion_inc_interE->Draw("HIST SAME");

   auto legend_interE = new TLegend(0.1,0.7,0.48,0.9);
   legend_interE->AddEntry(h_trueE_truePion_inc_interE,"True Pions, trueE, after BeamCuts");
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

   h_trueE_effPurCorr_inc_interE->Multiply(h_trueE_selPion_inc_interE, h_trueE_pur_removeBG_inc_interE );
   h_trueE_effPurCorr_inc_interE->Divide( h_trueE_eff_eventSel_inc_interE );

   //    Rebuild from initial and interacting Distribution the full incident unsmeared histogram
   //=====================================================

   TH1* trueE_cum_initE = h_trueE_effPurCorr_inc_initE->GetCumulative(false);
   TH1* trueE_cum_interE = h_trueE_effPurCorr_inc_interE->GetCumulative(false);
   //trueE_cum_initE->Write();
   //trueE_cum_interE->Write();

   //h_trueE_effPurCorr_incident->SetBinContent(nBin_int, trueE_cum_initE->GetBinContent(nBin_int));
   for(int i=nBin_int-1; i >= 1; i--){
      double diff = trueE_cum_initE->GetBinContent(i+1) - trueE_cum_interE->GetBinContent(i+1);
      h_trueE_effPurCorr_incident->SetBinContent( i, diff);
   };

   h_trueE_effPurCorr_incident->Write();

   h_trueE_effPurCorr_int->Multiply(h_trueE_selAbs_int, h_trueE_pur_removeBG_int );
   h_trueE_effPurCorr_int->Divide( h_trueE_eff_eventSel_int );

   h_trueE_effPurCorr_int->Write();

   TCanvas *c_trueE_effPur_compare_trueE_int = new TCanvas("c_trueE_effPur_compare_trueE_int", "");

   h_trueE_trueAbs_int->SetLineColor(kBlue);
   h_trueE_trueAbs_int->SetLineWidth(2);
   h_trueE_trueAbs_int->SetMarkerColor(kBlue);
   h_trueE_trueAbs_int->SetBarOffset(0.5);
   h_trueE_trueAbs_int->GetXaxis()->SetTitle("KE (MeV)");

   h_trueE_effPurCorr_int->SetLineColor(kGreen);
   h_trueE_effPurCorr_int->SetLineWidth(2);
   h_trueE_effPurCorr_int->SetFillColorAlpha( kGreen, 0.2);
   h_trueE_effPurCorr_int->SetMarkerColor(kGreen);
   h_trueE_effPurCorr_int->SetBarOffset(-0.5);

   h_trueE_trueAbs_int->Draw("HIST TEXT00");
   h_trueE_effPurCorr_int->Draw("HIST SAME");

   auto leg_effPurCorr_int = new TLegend(0.1,0.7,0.48,0.9);
   leg_effPurCorr_int->AddEntry(h_trueE_trueAbs_int,"True Absorption, trueE, after BeamCuts");
   leg_effPurCorr_int->AddEntry(h_trueE_effPurCorr_int,"Corrected Selected Abs in trueE");
   leg_effPurCorr_int->Draw();
   c_trueE_effPur_compare_trueE_int->Write();

   TCanvas *c_trueE_effPur_compare_trueE_inc = new TCanvas("c_trueE_effPur_compare_trueE_inc", "");

   h_trueE_truePion_incident->SetLineColor(kBlue);
   h_trueE_truePion_incident->SetLineWidth(2);
   h_trueE_truePion_incident->SetMarkerColor(kBlue);
   h_trueE_truePion_incident->SetMarkerSize(0.8);
   h_trueE_truePion_incident->SetBarOffset(0.5);
   h_trueE_truePion_incident->GetXaxis()->SetTitle("KE (MeV)");

   h_trueE_effPurCorr_incident->SetLineColor(kGreen);
   h_trueE_effPurCorr_incident->SetLineWidth(2);
   h_trueE_effPurCorr_incident->SetFillColorAlpha( kGreen, 0.2);
   h_trueE_effPurCorr_incident->SetMarkerColor(kGreen);
   h_trueE_effPurCorr_incident->SetBarOffset(-0.5);

   h_trueE_truePion_incident->Draw("HIST TEXT00");
   h_trueE_effPurCorr_incident->Draw("HIST SAME");

   auto leg_effPurCorr_inc = new TLegend(0.1,0.7,0.48,0.9);
   leg_effPurCorr_inc->AddEntry(h_trueE_truePion_incident,"True Pions, trueE, after BeamCuts");
   leg_effPurCorr_inc->AddEntry(h_trueE_effPurCorr_incident,"Corrected Selected Pion Incident in trueE");
   leg_effPurCorr_inc->Draw();
   c_trueE_effPur_compare_trueE_inc->Write();

   TCanvas *c_trueE_effPur_compare_trueE_initE = new TCanvas("c_trueE_effPur_compare_trueE_initE", "");

   h_trueE_truePion_inc_initE->SetLineColor(kBlue);
   h_trueE_truePion_inc_initE->SetLineWidth(2);
   h_trueE_truePion_inc_initE->SetMarkerColor(kBlue);
   h_trueE_truePion_inc_initE->SetMarkerSize(0.8);
   h_trueE_truePion_inc_initE->SetBarOffset(0.5);
   h_trueE_truePion_inc_initE->GetXaxis()->SetTitle("KE (MeV)");

   h_trueE_effPurCorr_inc_initE->SetLineColor(kGreen);
   h_trueE_effPurCorr_inc_initE->SetLineWidth(2);
   h_trueE_effPurCorr_inc_initE->SetFillColorAlpha( kGreen, 0.2);
   h_trueE_effPurCorr_inc_initE->SetMarkerColor(kGreen);
   h_trueE_effPurCorr_inc_initE->SetBarOffset(-0.5);

   h_trueE_truePion_inc_initE->Draw("HIST TEXT00");
   h_trueE_effPurCorr_inc_initE->Draw("HIST SAME");

   auto leg_effPurCorr_initE = new TLegend(0.1,0.7,0.48,0.9);
   leg_effPurCorr_initE->AddEntry(h_trueE_truePion_inc_initE,"True Pions, trueE, after BeamCuts");
   leg_effPurCorr_initE->AddEntry(h_trueE_effPurCorr_inc_initE,"Correcte Selected Pion initE in trueE");
   leg_effPurCorr_initE->Draw();

   c_trueE_effPur_compare_trueE_initE->Write();

   TCanvas *c_trueE_effPur_compare_trueE_interE = new TCanvas("c_trueE_effPur_compare_trueE_interE", "");

   h_trueE_truePion_inc_interE->SetLineColor(kBlue);
   h_trueE_truePion_inc_interE->SetLineWidth(2);
   h_trueE_truePion_inc_interE->SetMarkerColor(kBlue);
   h_trueE_truePion_inc_interE->SetMarkerSize(0.8);
   h_trueE_truePion_inc_interE->SetBarOffset(0.5);
   h_trueE_truePion_inc_interE->GetXaxis()->SetTitle("KE (MeV)");

   h_trueE_effPurCorr_inc_interE->SetLineColor(kGreen);
   h_trueE_effPurCorr_inc_interE->SetLineWidth(2);
   h_trueE_effPurCorr_inc_interE->SetFillColorAlpha( kGreen, 0.2);
   h_trueE_effPurCorr_inc_interE->SetMarkerColor(kGreen);
   h_trueE_effPurCorr_inc_interE->SetBarOffset(-0.5);

   h_trueE_truePion_inc_interE->Draw("HIST TEXT00");
   h_trueE_effPurCorr_inc_interE->Draw("HIST SAME");

   auto leg_effPurCorr_interE = new TLegend(0.1,0.7,0.48,0.9);
   leg_effPurCorr_interE->AddEntry(h_trueE_truePion_inc_interE,"True Pions, trueE, after BeamCuts");
   leg_effPurCorr_interE->AddEntry(h_trueE_effPurCorr_inc_interE,"Correcte Selected Pion interE in trueE");
   leg_effPurCorr_interE->Draw();

   c_trueE_effPur_compare_trueE_interE->Write();

   return 0;
}

