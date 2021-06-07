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

int test_unsmear(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");

   TFile *output = new TFile ( "output_test_unsmear_1.root" , "RECREATE");

   //TrueProcess and TrueE Int and Inc Histos
   TH1D* h_trueE_trueAbs_int = new TH1D("h_trueE_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_trueE_truePion_inc = new TH1D("h_trueE_truePion_inc", "", nBin_int, eEnd, eStart);
   //Selected Process and RecoE Int and Inc Histos
   TH1D* h_recoE_selAbs_int = new TH1D("h_recoE_selAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selPion_inc = new TH1D("h_recoE_selPion_inc", "", nBin_int, eEnd, eStart);
   //Unsmeared Histo, compare to h_trueE_trueProc
   TH1D* h_help_unsmear_int = new TH1D("h_help_unsmear_int", "", nBin_int, eEnd, eStart);
   TH1D* h_help_unsmear_inc = new TH1D("h_help_unsmear_inc", "", nBin_int, eEnd, eStart);
   
   TH1D* h_unsmeared_int = new TH1D("h_unsmeared_int", "", nBin_int, eEnd, eStart);
   TH1D* h_unsmeared_inc = new TH1D("h_unsmeared_inc", "", nBin_int, eEnd, eStart);

   //From evSel --> back to Reco MC Nj --> Nj'
   TH1D* h_pur_removeBG_int = new TH1D("h_pur_removeBG_int", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_int = new TH1D("h_eff_eventSel_int", "", nBin_int, eEnd, eStart);

   TH1D* h_pur_removeBG_inc = new TH1D("h_pur_removeBG_inc", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_inc = new TH1D("h_eff_eventSel_inc", "", nBin_int, eEnd, eStart);

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
   //Incident Histo
   eventSel_post_beamCut
      .Filter("true_primPionInel") //should also take into account pions that decay

      .Foreach( [h_trueE_truePion_inc] (double true_firstEntryIncident, double true_beam_interactingEnergy) {
            //bin that energy falls into is (int) energy/nbins + 1
            int binNumber_initEnergy = (int) true_firstEntryIncident / bin_size_int + 1;
            int binNumber_interEnergy = (int) true_beam_interactingEnergy / bin_size_int + 1;
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){                 

            if(i <= nBin_int){
            h_trueE_truePion_inc->SetBinContent( i, h_trueE_truePion_inc->GetBinContent(i) + 1 ); 
            };
            };
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP"});


   //Interacting Histo
   eventSel_post_beamCut
      .Filter("true_absSignal") //should also take into account pions that decay

      .Foreach( [h_trueE_trueAbs_int] (double true_beam_interactingEnergy){
            h_trueE_trueAbs_int->Fill(true_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});

   h_trueE_trueAbs_int->Write();
   h_trueE_truePion_inc->Write();
   //------------------------------------------------------------
   //
   //Create the selected and Reco histos that need to be unsmeared back
   //
   //
   //------------------------------------------------------------
   //
   eventSel_incidentPion
      .Foreach( [h_recoE_selPion_inc] (double reco_firstEntryIncident, double reco_beam_interactingEnergy) {
            int binNumber_initEnergy = (int) reco_firstEntryIncident / bin_size_int + 1;
            int binNumber_interEnergy = (int) reco_beam_interactingEnergy / bin_size_int + 1;
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
            
            if( i > 0 && i <= nBin_int){ //make sure we don't go outside of bin range
            h_recoE_selPion_inc->SetBinContent( i, h_recoE_selPion_inc->GetBinContent(i) + 1 ); 
            };
            };
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});

   eventSel_abs
      .Foreach( [h_recoE_selAbs_int] (double reco_beam_interactingEnergy){
            h_recoE_selAbs_int->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   h_recoE_selPion_inc->Write();
   h_recoE_selAbs_int->Write();
   //------------------------------------------------------------
   //
   //Start creating the efficiencies and purities to subtract the BG and account 
   //for lost events from eventSelection
   //
   //------------------------------------------------------------
   TH1D* h_recoE_incidentPion_truePion_inc = new TH1D("h_recoE_incidentPion_truePion_inc", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_beamCut_truePion_inc = new TH1D("h_recoE_beamCut_truePion_inc", "", nBin_int, eEnd, eStart);

   TH1D* h_recoE_selAbs_trueAbs_int = new TH1D("h_recoE_selAbs_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_beamCut_trueAbs_int = new TH1D("h_recoE_beamCut_trueAbs_int", "", nBin_int, eEnd, eStart);

   //True Pions available after beamCuts
   eventSel_post_beamCut
      .Filter("true_primPionInel")
      .Foreach( [h_recoE_beamCut_truePion_inc] (double reco_firstEntryIncident, double reco_beam_interactingEnergy) {
            int binNumber_initEnergy = (int) reco_firstEntryIncident / bin_size_int + 1;
            int binNumber_interEnergy = (int) reco_beam_interactingEnergy / bin_size_int + 1;
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
            if( i > 0 && i <= nBin_int){ //make sure we don't go outside of bin range
            h_recoE_beamCut_truePion_inc->SetBinContent( i, h_recoE_beamCut_truePion_inc->GetBinContent(i) + 1 ); 
            };
            };
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});

   //True Pions in selected incident
   eventSel_incidentPion
      .Filter("true_primPionInel")
      .Foreach( [h_recoE_incidentPion_truePion_inc] (double reco_firstEntryIncident, double reco_beam_interactingEnergy) {
            int binNumber_initEnergy = (int) reco_firstEntryIncident / bin_size_int + 1;
            int binNumber_interEnergy = (int) reco_beam_interactingEnergy / bin_size_int + 1;
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
            if( i > 0 && i <= nBin_int){ //make sure we don't go outside of bin range
            h_recoE_incidentPion_truePion_inc->SetBinContent( i, h_recoE_incidentPion_truePion_inc->GetBinContent(i) + 1 ); 
            };
            };
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});

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
   h_pur_removeBG_inc->Divide( h_recoE_incidentPion_truePion_inc, h_recoE_selPion_inc );
   h_eff_eventSel_inc->Divide( h_recoE_incidentPion_truePion_inc, h_recoE_beamCut_truePion_inc );

   h_pur_removeBG_int->Divide( h_recoE_selAbs_trueAbs_int, h_recoE_selAbs_int );
   h_eff_eventSel_int->Divide( h_recoE_selAbs_trueAbs_int, h_recoE_beamCut_trueAbs_int );

   h_pur_removeBG_inc->Write();
   h_eff_eventSel_inc->Write();

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


   TH2D* h2_smearing_incident = new TH2D("h2_smearing_incident", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   TH2D* h2_smearing_interacting = new TH2D("h2_smearing_interacting", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);


   //========================================================
   //Build the smearing Incident Histogram
   //---------
   //eventSel_incidentPion
   eventSel_post_beamCut
      .Filter("true_primPionInel")
      .Foreach( [h2_smearing_incident] (double true_inc, double true_int, double reco_inc, double reco_int) {

            int nBin_true = (int) (true_inc - true_int) / bin_size_int + 1;
            int nBin_reco = (int) (reco_inc - reco_int) / bin_size_int + 1;

            int array_size = 0;
            if(nBin_true > nBin_reco) array_size = nBin_true; //comparing length of reco vs true Vec
            else array_size = nBin_reco;

            double array_trueE[array_size], array_recoE[array_size];

            //fill arrays of same length with either energyTruevalue, reco or 0
            for(int i = 0; i < array_size; i++){

            if( true_inc - i * bin_size_int >= true_int ) array_trueE[i] = true_inc - i * bin_size_int; //while incident energy above interacting, it's fine, fill for inc - i*40MeV
            else array_trueE[i] = 0;

            if( reco_inc - i * bin_size_int >= reco_int ) array_recoE[i] = reco_inc - i * bin_size_int;
            else array_recoE[i] = 0;

            h2_smearing_incident->Fill( array_recoE[i], array_trueE[i]);

            };

            /*std::cout << "TrueE Vec Inc " << std::endl;
              for(auto i : array_trueE) std::cout << "trueE = " << i << std::endl;

              std::cout << "recoE Vec Inc " << std::endl;
              for(auto i : array_recoE) std::cout << "recoE = " << i << std::endl;*/
      }
   ,{"true_firstEntryIncident", "true_KEint_fromEndP", "reco_firstEntryIncident", "reco_interactingKE"});


   //NORMALISATION Normalise to RecoE column
   //
   //Go through true columns and normalise the entries 
   for(int i = 1; i <= nBin_int; i++){

      int sum_true_i = h2_smearing_incident->Integral( i, i ,1,nBin_int); //sum up all smeared trueSignals by integrating on 1RecoBin(x) over all True Bins(y)

      if(sum_true_i !=0){
         for(int j = 1; j <= nBin_int; j++){

            h2_smearing_incident->SetBinContent( i , j, h2_smearing_incident->GetBinContent(i,j) / sum_true_i);

         };
      }
   };

   h2_smearing_incident->Sumw2(0);
   h2_smearing_incident->SetTitle("Smearing Matrix for truePion in Sample after BeamCuts; reco Energy [MeV]; true Energy [MeV]");
   h2_smearing_incident->Write();


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
   //
   //Go through true columns and normalise the entries 
   for(int i = 1; i <= nBin_int; i++){

      int sum_true_i = h2_smearing_interacting->Integral( i, i ,1,nBin_int); //sum up all smeared trueSignals by integrating on 1TrueBin(x) over allReco Bins(y)

      if(sum_true_i !=0){
         for(int j = 1; j <= nBin_int; j++){

            h2_smearing_interacting->SetBinContent( i , j, h2_smearing_interacting->GetBinContent(i,j) / sum_true_i);

         };
      }
   };

   h2_smearing_interacting->Sumw2(0);
   h2_smearing_interacting->SetTitle("Smearing Matrix for trueAbs in Sample after BeamCuts; reco Energy [MeV]; true Energy [MeV]");
   h2_smearing_interacting->Write();

   //=====================================================
   //------------------------------------------------------
   //Start Unsmearing Incident
   //------------------------------------------------------
   //
   //    Remove the BGs and eventSelection inefficiency
   //    multiply by purity and divide by efficency
   h_help_unsmear_inc->Multiply(h_recoE_selPion_inc, h_pur_removeBG_inc );
   h_help_unsmear_inc->Divide( h_eff_eventSel_inc );
   
   //   Undo the smearing by multiplying Tij*Nj, so every row of true percentage with every entry of help_unsmeared
   //   Row's are i(y), columns are j(x) 
   //
   //Loop through smearing matrix Incident
   for(int i = 1; i <= nBin_int; i++){
      double help_sum = 0;
      //columns
      for(int j = 1; j <= nBin_int; j++){
         help_sum = help_sum + h2_smearing_incident->GetBinContent( j , i ) * h_help_unsmear_inc->GetBinContent(j) ;
      }; 
      h_unsmeared_inc->SetBinContent( i , help_sum);
   };

   h_unsmeared_inc->Write();

   
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
      //columns
      for(int j = 1; j <= nBin_int; j++){
         help_sum = help_sum + h2_smearing_interacting->GetBinContent( j , i ) * h_help_unsmear_int->GetBinContent(j) ;
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

   TCanvas *c_smear_inc = new TCanvas("c_smear_inc", "");
   h2_smearing_incident->SetMarkerSize(0.9);
   h2_smearing_incident->Draw("COLZ TEXT");
   c_smear_inc->Write();

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

   h_trueE_truePion_inc->SetLineColor(kBlue);
   h_trueE_truePion_inc->SetLineWidth(2);
   h_trueE_truePion_inc->SetMarkerColor(kBlue);
   h_trueE_truePion_inc->SetMarkerSize(0.8);
   h_trueE_truePion_inc->SetBarOffset(0.5);
   h_trueE_truePion_inc->GetXaxis()->SetTitle("KE (MeV)");

   h_unsmeared_inc->SetLineColor(kRed);
   h_unsmeared_inc->SetLineWidth(2);
   h_unsmeared_inc->SetFillColorAlpha( kRed, 0.2);
   h_unsmeared_inc->SetMarkerColor(kRed);
   h_unsmeared_inc->SetBarOffset(-0.5);

   h_recoE_selPion_inc->SetLineStyle(2);
   h_recoE_selPion_inc->SetLineWidth(3);
   h_recoE_selPion_inc->SetFillColorAlpha(kBlack, 0.1);

   h_trueE_truePion_inc->Draw("HIST TEXT00");
   h_unsmeared_inc->Draw("HIST SAME");
   h_recoE_selPion_inc->Draw("HIST SAME");

   auto legend_inc = new TLegend(0.1,0.7,0.48,0.9);
   legend_inc->AddEntry(h_trueE_truePion_inc,"True Pions, trueE, after BeamCuts");
   legend_inc->AddEntry(h_unsmeared_inc,"Unsmeared Incident");
   legend_inc->AddEntry(h_recoE_selPion_inc,"Selected Incident");
   legend_inc->Draw();
   c_smear_inc->Write();


   return 0;
}

