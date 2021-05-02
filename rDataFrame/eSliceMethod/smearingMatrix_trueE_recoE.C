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

int smearingMatrix_trueE_recoE(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");

   TFile *output = new TFile ( "output_smearMatrix_energy.root" , "RECREATE");

   //--------------------------------------------------------
   // Initialise smearing Matrix TH2D
   // NOT interested in BG smearing, so only consider truePrimaryPions  and trueAbsSignal
   // true_primPionInel is also including the decays
   //--------------------------------------------------------

   //Incident Histogram and Interacting Histogram, interacting for different Processes
   //Energies are kinetic Energy of beam particle
   //x-axis should be TrueE
   //y-axis should be RecoE
   
   TH2D* h2_smearing_incident = new TH2D("h2_smearing_incident", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   TH2D* h2_smearing_interacting = new TH2D("h2_smearing_interacting", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   
   //Initial Filters for all events
   auto mcIncident_selected_primaryPi_truePi = frame
      .Filter("true_primPionInel")
      .Filter("selected_incidentPion");

   auto mcInteracting_selected_abs_trueAbs = frame 
      .Filter("true_absSignal")
      .Filter("selected_abs");

   //========================================================
   //Build the smearing Incident Histogram
   //---------
   //take the interacting energy, that slice is filled too
        mcIncident_selected_primaryPi_truePi
         //.Range(30)
         .Foreach( [h2_smearing_incident] (double true_inc, double true_int, double reco_inc, double reco_int) {
               
               int nBin_true = (int) (true_inc - true_int) / bin_size_int + 1;
               int nBin_reco = (int) (reco_inc - reco_int) / bin_size_int + 1;

               int array_size = 0;
               if(nBin_true > nBin_reco) array_size = nBin_true;
               else array_size = nBin_reco;

               double array_trueE[array_size], array_recoE[array_size];

               //fill arrays of same length with either energyTruevalue, reco or 0
               for(int i = 0; i < array_size; i++){
                  
                  if( true_inc - i * bin_size_int >= true_int ) array_trueE[i] = true_inc - i * bin_size_int; //while incident energy above interacting, it's fine, fill for inc - i*40MeV
                  else array_trueE[i] = 0;

                  if( reco_inc - i * bin_size_int >= reco_int ) array_recoE[i] = reco_inc - i * bin_size_int;
                  else array_recoE[i] = 0;

                  h2_smearing_incident->Fill( array_trueE[i], array_recoE[i]);

               };
               
               /*std::cout << "TrueE Vec Inc " << std::endl;
               for(auto i : array_trueE) std::cout << "trueE = " << i << std::endl;

               std::cout << "recoE Vec Inc " << std::endl;
               for(auto i : array_recoE) std::cout << "recoE = " << i << std::endl;*/
         }
        ,{"true_firstEntryIncident", "true_KEint_fromEndP", "reco_firstEntryIncident", "reco_interactingKE"});


   //NORMALISATION Normalise to TrueE column
   //
   //Go through true columns and normalise the entries 
   for(int i = 1; i <= nBin_int; i++){

      int sum_true_i = h2_smearing_incident->Integral( i, i ,1,nBin_int); //sum up all smeared trueSignals by integrating on 1TrueBin(x) over allReco Bins(y)

      if(sum_true_i !=0){
         for(int j = 1; j <= nBin_int; j++){

            h2_smearing_incident->SetBinContent( i , j, h2_smearing_incident->GetBinContent(i,j) / sum_true_i);

         };
      }
   };

   h2_smearing_incident->Sumw2(0);
   h2_smearing_incident->SetTitle("Smearing Matrix for truePion in selected Incident Sample; true Energy [MeV]; reco Energy [MeV]");
   h2_smearing_incident->Write();


   //=====================================================
   //------------------------------------------------------
   //Interacting selected samples
   //------------------------------------------------------
   //
   //interacting samples are considered within the first APA, include APA3Cut bc of bad Energy reco after gap
   // --> For incident we still need to take into account that some interact 
   // further back in the TPC, so the cut on APA3 is not applied for the incident sample
   //
      mcInteracting_selected_abs_trueAbs
         .Foreach( [h2_smearing_interacting] (double true_int, double reco_int){

               h2_smearing_interacting->Fill( true_int, reco_int );}

               ,{"true_KEint_fromEndP", "reco_interactingKE"}); 

   //NORMALISATION Normalise to TrueE column
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
   h2_smearing_interacting->SetTitle("Smearing Matrix for trueAbs in selected interacting Sample; true Energy [MeV]; reco Energy [MeV]");
   h2_smearing_interacting->Write();

   //------------------------------------------------------

   //=====================================================*
   //            Plotting and Style
   //=====================================================
   //
   //h_xs_RecoE_selected_abs->Scale(2);
   //
   //

/*
   TCanvas *c_RecoE_abs = new TCanvas("c_RecoE_abs", "c_RecoE_abs");
   gPad->SetGrid(1,1);
   h_xs_RecoE_selected_abs->SetTitle( "Selected Absorption; Reco Kinetic Energy (MeV); #sigma (mb)");
   h_xs_RecoE_selected_abs->GetXaxis()->SetRangeUser(400,900);
   h_xs_RecoE_selected_abs->GetXaxis()->SetNdivisions(1020);
   h_xs_RecoE_selected_abs->GetYaxis()->SetNdivisions(1020);

   string c_abs_title;
   if(isMC) c_abs_title = "Absorption MC; Reco kinetic Energy (MeV); #sigma (mbarn)";
   else c_abs_title = "Absorption Data; Reco kinetic Energy (MeV); #sigma (mbarn)";

   abs_KE->SetTitle( c_abs_title.c_str());
   abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   abs_KE->SetLineColor(kBlue);
   abs_KE->SetLineWidth(3);
   abs_KE->Draw("AC");
   h_xs_RecoE_selected_abs->SetMarkerSize(0.7);
   h_xs_RecoE_selected_abs->Draw("PE0 SAME");

   c_RecoE_abs->Write();


      TCanvas *c_TrueE_abs = new TCanvas("c_TrueE_abs", "c_TrueE_abs");
      gPad->SetGrid(1,1);
      h_xs_TrueE_selected_abs->SetTitle( "Selected Absorption;True Kinetic Energy (MeV); #sigma (mb)");
      h_xs_TrueE_selected_abs->GetXaxis()->SetRangeUser(400,900);
      h_xs_TrueE_selected_abs->GetXaxis()->SetNdivisions(1020);
      h_xs_TrueE_selected_abs->GetYaxis()->SetNdivisions(1020);

      abs_KE->SetTitle( "Absorption;True Kinetic Energy (MeV); #sigma (mb)");
      abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
      abs_KE->SetLineColor(kRed);
      abs_KE->SetLineWidth(3);
      abs_KE->Draw("AC");
      h_xs_TrueE_selected_abs->SetMarkerSize(0.7);
      h_xs_TrueE_selected_abs->Draw("PE0 SAME");

      c_TrueE_abs->Write();

  
*/

   //output->Write();
   //f1.Close();
   return 0;
}

