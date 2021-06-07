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

int unsmear_test(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");

   nBin_int = 12; //100MeV bins
   TFile *output = new TFile ( "100MeV/output_unsmear_test.root" , "RECREATE");

   //TrueProcess and TrueE Int and Inc Histos
   TH1D* h_trueE_trueAbs_int = new TH1D("h_trueE_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_selAbs_int = new TH1D("h_recoE_selAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_help_unsmear_int = new TH1D("h_help_unsmear_int", "", nBin_int, eEnd, eStart);
   TH1D* h_unsmeared_int = new TH1D("h_unsmeared_int", "", nBin_int, eEnd, eStart);
   //From evSel --> back to Reco MC Nj --> Nj'
   TH1D* h_pur_removeBG_int = new TH1D("h_pur_removeBG_int", "", nBin_int, eEnd, eStart);
   TH1D* h_eff_eventSel_int = new TH1D("h_eff_eventSel_int", "", nBin_int, eEnd, eStart);

   //Build the True Process and TrueE Int and Inc Histograms that we need to compare unsmeared things to
   //
   //all available after beamCuts
   auto eventSel_post_beamCut = frame.Filter("primary_isBeamType && passBeamCut && passBeamCutBI");
   //selected incident Pions & selected absorption
   auto eventSel_incidentPion = frame.Filter("selected_incidentPion");
   auto eventSel_abs = frame.Filter("selected_abs");



   //Interacting Histo
   eventSel_post_beamCut
      .Filter("true_absSignal") //should also take into account pions that decay

      .Foreach( [h_trueE_trueAbs_int] (double true_beam_interactingEnergy){
            h_trueE_trueAbs_int->Fill(true_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});

   h_trueE_trueAbs_int->Write();

   eventSel_abs
      .Foreach( [h_recoE_selAbs_int] (double reco_beam_interactingEnergy){
            h_recoE_selAbs_int->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   h_recoE_selAbs_int->Write();
   //------------------------------------------------------------
   //
   //Start creating the efficiencies and purities to subtract the BG and account 
   //for lost events from eventSelection
   //
   //------------------------------------------------------------
   TH1D* h_recoE_selAbs_trueAbs_int = new TH1D("h_recoE_selAbs_trueAbs_int", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_beamCut_trueAbs_int = new TH1D("h_recoE_beamCut_trueAbs_int", "", nBin_int, eEnd, eStart);


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
   h_pur_removeBG_int->Divide( h_recoE_selAbs_trueAbs_int, h_recoE_selAbs_int );
   h_eff_eventSel_int->Divide( h_recoE_selAbs_trueAbs_int, h_recoE_beamCut_trueAbs_int );

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
   //Fill smearing matrix with entry 1 on the diagonal to avoid them being "uninvertible"
   for(int i = 1; i <= nBin_int; i++){
      h2_smearing_interacting->SetBinContent(i,i,1);
   }
   
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
   //Go through true rows and normalise the entries 
   for(int i = 1; i <= nBin_int; i++){

      int sum_true_i = h2_smearing_interacting->Integral( 1 , nBin_int, i, i ) ; //sum up all recoE signals for one trueE signal

      if(sum_true_i !=0){
         for(int j = 1; j <= nBin_int; j++){

            h2_smearing_interacting->SetBinContent( j, i, h2_smearing_interacting->GetBinContent(j,i) / sum_true_i);

         };
      }
   };

   h2_smearing_interacting->Sumw2(0);
   h2_smearing_interacting->SetTitle("Smearing Matrix for trueAbs in Sample after BeamCuts; reco Energy [MeV]; true Energy [MeV]");
   h2_smearing_interacting->Write();

   TMatrixD mymatrix_pre(nBin_int + 2, nBin_int + 2, h2_smearing_interacting->GetArray(), "D");
   
   //have to exclude underflow and overflowbin bc otherwise matrix is singular
   TMatrixD mymatrix = mymatrix_pre.GetSub(1, nBin_int, 1, nBin_int); //avoid underflow and overflowbin
   
   mymatrix.Print();

   Double_t det2;
   TMatrixD inverse = mymatrix;
   inverse.Invert(&det2);
 
   TMatrixD U2(inverse,TMatrixD::kMult, mymatrix);
   TMatrixDDiag diag2(U2); diag2 = 0.0;
   const Double_t U2_max_offdiag = (U2.Abs()).Max();
   std::cout << "  Maximum off-diagonal = " << U2_max_offdiag << std::endl;
   std::cout << "  Determinant          = " << det2 << std::endl;   

   //inverse.Print();
   TH2D *h2_inverted = new TH2D("h2_inverted", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
        for (int i = 1; i <= nBin_int; i++){
            for (int j= 1; j <= nBin_int; j++){
                h2_inverted->SetBinContent(j, i, inverse(i-1,j-1)); //vector indices style for matrix 
            };
        };
   
   inverse.Print();
   h2_inverted->Write();



   TMatrixD product_inv_mymatrix = inverse*mymatrix;
        
   TH2D *h2_prod_inv_mymatrix = new TH2D("h2_prod_inv_mymatrix", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
        for (int i = 1; i <= nBin_int; i++){
            for (int j= 1; j<= nBin_int; j++){
                h2_prod_inv_mymatrix->SetBinContent(j, i, product_inv_mymatrix(i-1,j-1)); //avoid overflow and underflow bin
            };
        };

        h2_prod_inv_mymatrix->Write();

   //------------------------------------------------------
   //Start Unsmearing Interacting
   //------------------------------------------------------
   //
   //    Remove the BGs and eventSelection inefficiency
   //    multiply by purity and divide by efficency
   h_help_unsmear_int->Multiply(h_recoE_selAbs_int, h_pur_removeBG_int );
   h_help_unsmear_int->Divide( h_eff_eventSel_int );

   //Loop through smearing matrix Interacting
   for(int i = 1; i <= nBin_int; i++){
      double help_sum = 0;
      //rows
      for(int j = 1; j <= nBin_int; j++){

         //help_sum += h_help_unsmear_int->GetBinContent(j)*inverse(j-1,i-1); //sum over rows inverse ==> rows are reco, columns are true)
         help_sum += h_help_unsmear_int->GetBinContent(j)*h2_inverted->GetBinContent( i, j);
      }; 
      h_unsmeared_int->SetBinContent( i , help_sum);
   };
   h_unsmeared_int->Write();

   //------------------------------------------------------
   // Smear the True and try to get the reco
   //------------------------------------------------------
   //
   TH1D* try_smear = new TH1D("try_smear", "", nBin_int, eEnd, eStart);
   TH1D* try_smear_unsmear = new TH1D("try_smear_unsmear", "", nBin_int, eEnd, eStart);

   for(int i=1; i<=nBin_int; i++){

      double sum_try=0;

      for(int j=1; j<= nBin_int; j++){
         
         sum_try += h_trueE_trueAbs_int->GetBinContent(j) * h2_smearing_interacting->GetBinContent( i, j);

      };
      try_smear->SetBinContent(i, sum_try);
      try_smear_unsmear->SetBinContent(i, sum_try);
   };

   //unsmear the smeared
   for(int i = 1; i <= nBin_int; i++){
      double help_sum = 0;
      //rows
      for(int j = 1; j <= nBin_int; j++){

         //help_sum += h_help_unsmear_int->GetBinContent(j)*inverse(j-1,i-1); //sum over rows inverse ==> rows are reco, columns are true)
         help_sum += try_smear->GetBinContent(j)*h2_inverted->GetBinContent( i, j);
      }; 
      try_smear_unsmear->SetBinContent( i , help_sum);
   };
   try_smear_unsmear->Write();



   //try_smear->Multiply(h_eff_eventSel_int );
   //try_smear->Divide(h_pur_removeBG_int);

   try_smear->Write();
   

   //------------------------------------------------------
   //          CANVAS
   //------------------------------------------------------
   //
  /*
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
*/

   return 0;
}

