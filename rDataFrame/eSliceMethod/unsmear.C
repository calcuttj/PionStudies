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
   ROOT::RDataFrame frame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");

   TFile *output = new TFile ( "100MeV/output_unsmear.root" , "RECREATE");

   //TrueProcess and TrueE Int and Inc Histos
   TH1D* h_trueE_trueAbs_int = new TH1D("h_trueE_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_inc_initE = new TH1D("h_trueE_truePion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_inc_interE = new TH1D("h_trueE_truePion_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_incident = new TH1D("h_trueE_truePion_incident", "", nBin_int, eEnd, eStart);
   //Selected Process and RecoE Int and Inc Histos
   TH1D* h_recoE_selAbs_int = new TH1D("h_recoE_selAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_inc_initE = new TH1D("h_recoE_selPion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_inc_interE = new TH1D("h_recoE_selPion_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_incident = new TH1D("h_recoE_selPion_incident", "", nBin_int, eEnd, eStart);
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

   //Build the True Process and TrueE Int and Inc Histograms that we need to compare unsmeared things to
   //
   //all available after beamCuts
   auto eventSel_post_beamCut = frame.Filter("primary_isBeamType && passBeamCut && passBeamCutBI");
   //selected incident Pions & selected absorption
   auto eventSel_incidentPion = frame.Filter("selected_incidentPion");
   auto eventSel_abs = frame.Filter("selected_abs");

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
      .Filter("true_primPionInel") //should also take into account pions that decay

      .Foreach( [h_trueE_truePion_inc_initE] (double true_firstEntryIncident) { 
            h_trueE_truePion_inc_initE->Fill( true_firstEntryIncident ); 
            }
            ,{"true_firstEntryIncident"});

   eventSel_post_beamCut
      .Filter("true_primPionInel") //should also take into account pions that decay

      .Foreach( [h_trueE_truePion_inc_interE] (double true_beam_interactingEnergy) { 
            h_trueE_truePion_inc_interE->Fill( true_beam_interactingEnergy ); 
            }
            ,{"true_KEint_fromEndP"});

   //Build trueIncident for comparison
    TH1* h_cum_trueE_initE = h_trueE_truePion_inc_initE->GetCumulative(false);
    TH1* h_cum_trueE_interE = h_trueE_truePion_inc_interE->GetCumulative(false);
    //h_cum_trueE_initE->Write();
    //h_cum_trueE_interE->Write();
    h_trueE_truePion_incident->SetBinContent(nBin_int, h_cum_trueE_initE->GetBinContent(nBin_int));
    for(int i=nBin_int; i > 1; i--){
       double diff = h_cum_trueE_initE->GetBinContent(i-1) - h_cum_trueE_interE->GetBinContent(i);
       h_trueE_truePion_incident->SetBinContent( i-1, diff);
    };

    h_trueE_truePion_incident->Write();

   //Interacting Histo
   eventSel_post_beamCut
      .Filter("true_absSignal") //should also take into account pions that decay

      .Foreach( [h_trueE_trueAbs_int] (double true_beam_interactingEnergy){
            h_trueE_trueAbs_int->Fill(true_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});

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
   eventSel_incidentPion
      .Foreach( [h_recoE_selPion_inc_initE] (double reco_firstEntryIncident) {
            h_recoE_selPion_inc_initE->Fill( reco_firstEntryIncident );            
            }
            ,{"reco_firstEntryIncident"});

   eventSel_incidentPion
      .Foreach( [h_recoE_selPion_inc_interE] (double reco_beam_interactingEnergy) {
            h_recoE_selPion_inc_interE->Fill( reco_beam_interactingEnergy );  
            }
            ,{"reco_interactingKE"});

   //Build recoIncident for comparison
    TH1* h_cum_recoE_initE = h_recoE_selPion_inc_initE->GetCumulative(false);
    TH1* h_cum_recoE_interE = h_recoE_selPion_inc_interE->GetCumulative(false);
    //h_cum_recoE_initE->Write();
    //h_cum_recoE_interE->Write();
    h_recoE_selPion_incident->SetBinContent(nBin_int, h_cum_recoE_initE->GetBinContent(nBin_int));
    for(int i=nBin_int; i > 1; i--){
       double diff = h_cum_recoE_initE->GetBinContent(i-1) - h_cum_recoE_interE->GetBinContent(i);
       h_recoE_selPion_incident->SetBinContent( i-1, diff);
    };

    h_recoE_selPion_incident->Write();

   eventSel_abs
      .Foreach( [h_recoE_selAbs_int] (double reco_beam_interactingEnergy){
            h_recoE_selAbs_int->Fill(reco_beam_interactingEnergy);
            }
            ,{"reco_interactingKE"});


   h_recoE_selPion_inc_interE->Write();
   h_recoE_selPion_inc_initE->Write();
   h_recoE_selAbs_int->Write();

   //------------------------------------------------------------
   //
   //    Test to get back original incident histo
   //------------------------------------------------------------


   // TH1* h_trueE_cum_initE = h_trueE_truePion_inc_initE->GetCumulative(false);
   // TH1* h_trueE_cum_interE = h_trueE_truePion_inc_interE->GetCumulative(false);
   // h_trueE_cum_initE->Write();
   // h_trueE_cum_interE->Write();

   // TH1D* h_trueE_inc_cumWay = new TH1D("h_trueE_inc_cumWay", "", nBin_int, eEnd, eStart);

   // h_trueE_inc_cumWay->SetBinContent(nBin_int, h_trueE_cum_initE->GetBinContent(nBin_int));
   // for(int i=nBin_int; i > 1; i--){
   //    double diff = h_trueE_cum_initE->GetBinContent(i-1) - h_trueE_cum_interE->GetBinContent(i);
   //    h_trueE_inc_cumWay->SetBinContent( i-1, diff);
   // };

   // h_trueE_inc_cumWay->Write();




   //------------------------------------------------------------
   //
   //Start creating the efficiencies and purities to subtract the BG and account 
   //for lost events from eventSelection
   //
   //------------------------------------------------------------
   TH1D* h_recoE_incidentPion_truePion_inc_initE = new TH1D("h_recoE_incidentPion_truePion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_beamCut_truePion_inc_initE = new TH1D("h_recoE_beamCut_truePion_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_incidentPion_truePion_inc_interE = new TH1D("h_recoE_incidentPion_truePion_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_beamCut_truePion_inc_interE = new TH1D("h_recoE_beamCut_truePion_inc_interE", "", nBin_int, eEnd, eStart);

   TH1D* h_recoE_selAbs_trueAbs_int = new TH1D("h_recoE_selAbs_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_beamCut_trueAbs_int = new TH1D("h_recoE_beamCut_trueAbs_int", "", nBin_int, eEnd, eStart);

   //True Pions available after beamCuts
   eventSel_post_beamCut
      .Filter("true_primPionInel")
      .Foreach( [h_recoE_beamCut_truePion_inc_initE] (double reco_firstEntryIncident) {
            h_recoE_beamCut_truePion_inc_initE->Fill( reco_firstEntryIncident );            
            }
            ,{"reco_firstEntryIncident"});

   eventSel_post_beamCut
      .Filter("true_primPionInel")
      .Foreach( [h_recoE_beamCut_truePion_inc_interE] (double reco_beam_interactingEnergy) {
            h_recoE_beamCut_truePion_inc_interE->Fill( reco_beam_interactingEnergy );            
            }
            ,{"reco_interactingKE"});

   //True Pions in selected incident
   eventSel_incidentPion
      .Filter("true_primPionInel")
      .Foreach( [h_recoE_incidentPion_truePion_inc_initE] (double reco_firstEntryIncident) {
            h_recoE_incidentPion_truePion_inc_initE->Fill( reco_firstEntryIncident );            
            }
            ,{"reco_firstEntryIncident"});

   eventSel_incidentPion
      .Filter("true_primPionInel")
      .Foreach( [h_recoE_incidentPion_truePion_inc_interE] (double reco_beam_interactingEnergy) {
            h_recoE_incidentPion_truePion_inc_interE->Fill( reco_beam_interactingEnergy );            
            }
            ,{"reco_interactingKE"});



   //True Abs available after beamCuts
   eventSel_post_beamCut
      .Filter("true_absSignal")
      .Foreach( [h_recoE_beamCut_trueAbs_int] (double reco_beam_interactingEnergy){
            h_recoE_beamCut_trueAbs_int->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   //True Abs available after beamCuts
   eventSel_abs
      .Filter("true_absSignal")
      .Foreach( [h_recoE_selAbs_trueAbs_int] (double reco_beam_interactingEnergy){
            h_recoE_selAbs_trueAbs_int->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

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
   }

   //========================================================
   //Build the smearing Incident Histogram initE
   //---------

   eventSel_post_beamCut 
      .Filter("true_primPionInel")
      .Foreach( [h2_smearing_incident_initE] (double true_initE, double reco_initE){

            h2_smearing_incident_initE->Fill( reco_initE, true_initE );}

            ,{"true_firstEntryIncident", "reco_firstEntryIncident"}); 

   //NORMALISATION Normalise to trueE row
   //
   for(int i = 1; i <= nBin_int; i++){

      //sum for same x (recoBin
      int sum_true_i = h2_smearing_incident_initE->Integral( 1, nBin_int, i, i ); //sum up all recoE signals for one trueE signal

      if(sum_true_i !=0){
         for(int j = 1; j <= nBin_int; j++){

            h2_smearing_incident_initE->SetBinContent( j , i, h2_smearing_incident_initE->GetBinContent(j,i) / sum_true_i);

         };
      }
   };

   h2_smearing_incident_initE->Sumw2(0);
   h2_smearing_incident_initE->Write();

   //TMatrix
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

   //========================================================
   //Build the smearing Incident Histogram interE
   //---------

   eventSel_post_beamCut 
      .Filter("true_primPionInel")
      .Foreach( [h2_smearing_incident_interE] (double true_interE, double reco_interE){

            h2_smearing_incident_interE->Fill( reco_interE, true_interE );}

            ,{"true_KEint_fromEndP", "reco_interactingKE"}); 


   for(int i = 1; i <= nBin_int; i++){

      int sum_true_i = h2_smearing_incident_interE->Integral( 1, nBin_int, i, i); //sum up all recoE signals for one trueE signal
      if(sum_true_i !=0){
         for(int j = 1; j <= nBin_int; j++){

            h2_smearing_incident_interE->SetBinContent( j , i, h2_smearing_incident_interE->GetBinContent(j,i) / sum_true_i);

         };
      }
   };

   h2_smearing_incident_interE->Sumw2(0);
   h2_smearing_incident_interE->Write();

   //TMatrix
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
   eventSel_post_beamCut 
      .Filter("true_absSignal")
      .Foreach( [h2_smearing_interacting] (double true_int, double reco_int){

            h2_smearing_interacting->Fill( reco_int, true_int );}

            ,{"true_KEint_fromEndP", "reco_interactingKE"}); 

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
   //Loop through smearing matrix Interacting
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

    h_unsmeared_incident->SetBinContent(nBin_int, h_cum_initE->GetBinContent(nBin_int));
    for(int i=nBin_int; i > 1; i--){
       double diff = h_cum_initE->GetBinContent(i-1) - h_cum_interE->GetBinContent(i);
       h_unsmeared_incident->SetBinContent( i-1, diff);
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
   TCanvas *c_smear_int = new TCanvas("c_smear_int", "");
   h2_smearing_interacting->SetMarkerSize(0.9);
   h2_smearing_interacting->Draw("COLZ TEXT");
   c_smear_int->Write();

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
   c_smear_int->Write();

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

   return 0;
}

