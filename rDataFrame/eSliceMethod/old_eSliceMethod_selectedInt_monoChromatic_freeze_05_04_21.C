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

int eSliceMethod_selectedInt_monoChromatic_freeze_05_04_21(const string mcFilepath, bool isMC){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);

   //output file
   string outputNameData = "output_eSliceMethod_selectedEvents_DATA_effPur_monoChromatic_freeze_05_04_21.root";
   string outputNameMC = "output_eSliceMethod_selectedEvents_MC_effPurr_monoChromatic_freeze_05_04_21.root";
   string outputName;
   if(isMC) outputName = outputNameMC;
   else outputName = outputNameData;

   
   TFile *output = new TFile ( outputName.c_str() , "RECREATE");

   //access Jakes GeantFile in folder
   TFile f1("exclusive_xsec.root");
   TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
   TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
   TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
   f1.Close();

   //Get Efficiencies and Purities for the selections bin-by-bin
   TFile f2("eSliceMethod_effPur_binByBin_monoChromatic_freeze_05_04_21.root");
   TH1D *h_eff_interacting_abs = (TH1D*)f2.Get("h_eff_interacting_abs");
   TH1D *h_pur_interacting_abs = (TH1D*)f2.Get("h_pur_interacting_abs");
   TH1D *h_eff_selection_abs = (TH1D*)f2.Get("h_eff_selection_abs");

   TH1D *h_eff_selection_incidentPion = (TH1D*)f2.Get("h_eff_selection_incidentPion");
   TH1D *h_eff_incidentPion = (TH1D*)f2.Get("h_eff_incidentPion");
   TH1D *h_pur_incidentPion = (TH1D*)f2.Get("h_pur_incidentPion");


   string mg_title;
   if(isMC) mg_title = "Cross-Section MC; Reco kinetic Energy (MeV); #sigma (mbarn)";
   else mg_title = "Cross-Section Data; Reco kinetic Energy (MeV); #sigma (mbarn)";
   TMultiGraph *mg = new TMultiGraph();
   mg->Add(totInel_KE);
   mg->Add(abs_KE);
   mg->Add(cex_KE);
   mg->SetTitle(mg_title.c_str());

   //switch to output-file
   output->cd();

   //--------------------------------------------------------
   // Initialise incident and interacting Histograms
   //
   //--------------------------------------------------------

   //Incident Histogram and Interacting Histogram, interacting for different Processes
   //Energies are kinetic Energy of beam particle

   TH1D* h_selected_pion_incidentRecoE = new TH1D("h_selected_pion_incidentRecoE", "Incident Selected Pion", nBin_inc, eEnd, eStart);
   h_selected_pion_incidentRecoE->GetXaxis()->SetTitle("Reco KE (MeV)");

   TH1D* h_selected_abs_interactingRecoE = new TH1D("h_selected_abs_interactingRecoE", "Interacting Selected ABS Reco Energy", nBin_int, eEnd, eStart);
   h_selected_abs_interactingRecoE->GetXaxis()->SetTitle("Reco KE (MeV)");

      TH1D* h_selected_pion_incidentTrueE = new TH1D("h_selected_pion_incidentTrueE", "Incident Selected Pion", nBin_int, eEnd, eStart);
      h_selected_pion_incidentTrueE->GetXaxis()->SetTitle("True KE (MeV)");

      TH1D* h_selected_abs_interactingTrueE = new TH1D("h_selected_abs_interactingTrueE", "Interacting Selected ABS True Energy", nBin_int, eEnd, eStart);
      h_selected_abs_interactingTrueE->GetXaxis()->SetTitle("True KE (MeV)");

   //Initial Filters for all events
   auto mcIncident_selected_primaryPi = frame
      //.Filter("true_beam_PDG == abs(211)")
      .Filter("selected_incidentPion");

   //Filter for different interaction types
   auto h2_wire_KEinc = mcIncident_selected_primaryPi
      .Histo2D({"recoKEinc_vs_wire", "reco first incident KEnergy at reco wire", 400,0,800, 200,0,1000}, "reco_interacting_wire", "reco_firstEntryIncident");

   h2_wire_KEinc->Write();

   auto h2_wire_KEint = mcIncident_selected_primaryPi
      .Histo2D({"recoKEint_vs_wire", "reco interacting KEnergy at reco wire", 400,0,800, 200,0,1000}, "reco_interacting_wire", "reco_interactingKE");

   h2_wire_KEint->Write();

   if(isMC){
      auto h2_reco_true_KE_primPi = mcIncident_selected_primaryPi
         .Histo2D({"h2_reco_true_KE_primPi", "true KE vs reco KE incident Prim Pi", 400,0,1200, 400,0,1200}, "true_KEint_fromEndP", "reco_interactingKE");
      h2_reco_true_KE_primPi->Write();
   }

   //========================================================
   //Build the Incident Histogram
   //---------
   //take the interacting energy, that slice is filled too
   if(isMC){
      mcIncident_selected_primaryPi
         .Foreach( [h_selected_pion_incidentTrueE] (double true_firstEntryIncident, double true_beam_interactingEnergy) {
               //bin that energy falls into is (int) energy/nbins + 1
               int binNumber_initEnergy = (int) true_firstEntryIncident / bin_size_int + 1;
               int binNumber_interEnergy = (int) true_beam_interactingEnergy / bin_size_int + 1;
               for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){                 

               if(i <= nBin_int){
               h_selected_pion_incidentTrueE->SetBinContent( i, h_selected_pion_incidentTrueE->GetBinContent(i) + 1 ); 
               //increment previous bin content by one. should handle number of entries in hist. unlike AddBinContent
               };
               //h_true_pion_true_incidentE->Fill( true_firstEntryIncident - cnt*bin_size);
               };
               }
               ,{"true_firstEntryIncident", "true_KEint_fromEndP"});
   }
   mcIncident_selected_primaryPi
      .Foreach( [h_selected_pion_incidentRecoE] (double reco_firstEntryIncident, double reco_beam_interactingEnergy) {
            //bin that energy falls into is (int) energy/nbins + 1
            int binNumber_initEnergy = (int) reco_firstEntryIncident / bin_size_inc + 1;
            int binNumber_interEnergy = (int) reco_beam_interactingEnergy / bin_size_inc + 1;
            //if(binNumber_initEnergy < 0 || binNumber_interEnergy < 0) return;
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
            if( i > 0 && i <= nBin_inc){ //make sure we don't go outside of bin range
            h_selected_pion_incidentRecoE->SetBinContent( i, h_selected_pion_incidentRecoE->GetBinContent(i) + 1 ); 
            //increment previous bin content by one. should handle number of entries in hist. unlike AddBinContent
            };
            };
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});
   if(isMC){
      h_selected_pion_incidentTrueE->Sumw2(0);
      h_selected_pion_incidentTrueE->Rebin( bin_size_int/bin_size_int );
      h_selected_pion_incidentTrueE->Scale( 1 / (bin_size_int/bin_size_int) );
      h_selected_pion_incidentTrueE->Sumw2(0);
      h_selected_pion_incidentTrueE->Write();
   }   
   h_selected_pion_incidentRecoE->Sumw2(0);
   h_selected_pion_incidentRecoE->Rebin( bin_size_int/bin_size_inc );
   h_selected_pion_incidentRecoE->Scale( 1 / (bin_size_int/bin_size_inc) );
   h_selected_pion_incidentRecoE->Sumw2(0);



   //=====================================================
   //------------------------------------------------------
   //Interacting selected samples
   //------------------------------------------------------
   //
   //interacting samples are considered within the first APA, include APA3Cut bc of bad Energy reco after gap
   // --> For incident we still need to take into account that some interact 
   // further back in the TPC, so the cut on APA3 is not applied for the incident sample
   //
   auto mcInteracting_selected_allPrimaryPi = mcIncident_selected_primaryPi
      .Filter("primary_ends_inAPA3");


   auto mcInteracting_selected_abs = mcInteracting_selected_allPrimaryPi
      .Filter("selected_abs");
   //.Filter("has_noPion_daughter && !(has_shower_nHits_distance)");
   if(isMC){
      mcInteracting_selected_abs
         .Foreach( [h_selected_abs_interactingTrueE] (double true_beam_interactingEnergy){
               h_selected_abs_interactingTrueE->Fill(true_beam_interactingEnergy);}
               ,{"true_KEint_fromEndP"}); 

      h_selected_abs_interactingTrueE->Sumw2(0);
      h_selected_abs_interactingTrueE->Write();
   }
   mcInteracting_selected_abs 
      .Foreach( [h_selected_abs_interactingRecoE] (double reco_beam_interactingEnergy){
            h_selected_abs_interactingRecoE->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   h_selected_abs_interactingRecoE->Sumw2(0);
   //h_selected_abs_interactingRecoE->Write();


   //=====================================================
   //            Prepare BetheBloch Mean for each Bin 
   //            QUESTION: take betheBloch of Pion or Muon?? Comparison to data fits better muon Bethe... 
   //            at hihger momentum ~400-800 they anway are almost the same
   //=====================================================
   TH1D* h_betheMean_muon = new TH1D("h_betheMean_muon", "Mean Energy Loss", nBin_int, eEnd, eStart);

   //fill histo with Mean dEdX of bin center
   for(int i = 1; i <= nBin_int; i++){
      h_betheMean_muon->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size_int  , mass_muon) );
   };

   h_betheMean_muon->Write();

   //------------------------------------------------------
   //    Multiply Scale the Incident and Interacting Histograms 
   //    with their purities and efficiencies bin-by-bin
   // 
   //    hist * (purity/eff) on a bin-by-bin basis
   //------------------------------------------------------
   //Test Multiply Incident histo by it's purity and divide by it's efficiency not bin-by-bin, just testing overall
   //h_selected_pion_incidentRecoE->Scale(0.83/1);
   //Multiply Interacting by it's purity and divide by efficiency wrt to incident slection (I think..)
   //h_selected_abs_interactingRecoE->Scale(0.54 / 0.54 );

   h_selected_pion_incidentRecoE->Multiply( h_pur_incidentPion );
   h_selected_pion_incidentRecoE->Divide( h_eff_incidentPion );
   //Take into account the efficiency of losing true pions wrt beamCuts
   h_selected_pion_incidentRecoE->Divide( h_eff_selection_incidentPion );

   //Take account for purity of TrueE == RecoE and Efficiency from Smearing
   h_selected_abs_interactingRecoE->Multiply( h_pur_interacting_abs );
   h_selected_abs_interactingRecoE->Divide( h_eff_interacting_abs );
   //Take into account efficiency coming from selection, i.e. abs available at incident Pion Sample
   h_selected_abs_interactingRecoE->Divide( h_eff_selection_abs );
   h_selected_pion_incidentRecoE->Sumw2(0);
   h_selected_pion_incidentRecoE->Write();
   
   h_selected_abs_interactingRecoE->Sumw2(0);
   h_selected_abs_interactingRecoE->Write();

   //=====================================================
   //             Computing the XS
   //=====================================================
   //
   // xs(Ebin) = (A / (Na*density*bin_size)) * dEdX(Ebin) * hInteracting / hIncident
   //
   // More Accurate use log( Ninc / (Ninc - Nint ))
   //
   double scale_factor = factor_mbarn * atomic_mass / ( density * N_avogadro * bin_size_int );

   //------------------------------------------------------
   //    Absorption, Selected Interactions Reconstrucetd Energy
   //------------------------------------------------------
   //
   string xs_abs_name, xs_abs_title;
   if(isMC) {
      xs_abs_name = "h_xs_recoE_selected_abs_mc";
      xs_abs_title = "Absorption MC";
   }
   else {
      xs_abs_name = "h_xs_recoE_selected_abs_data";
      xs_abs_title = "Absorption Data";
   }


   TH1D* h_xs_RecoE_selected_abs = (TH1D*) h_selected_pion_incidentRecoE->Clone(xs_abs_name.c_str());
   h_xs_RecoE_selected_abs->SetTitle(xs_abs_title.c_str());
   TH1D* dummy_recoE_abs = (TH1D*) h_selected_pion_incidentRecoE->Clone(xs_abs_name.c_str());
   dummy_recoE_abs->Add( h_selected_abs_interactingRecoE, -1);
   h_xs_RecoE_selected_abs->Divide( dummy_recoE_abs );
   for(int i = 1; i <= nBin_int; i++) h_xs_RecoE_selected_abs->SetBinContent(i, log( h_xs_RecoE_selected_abs->GetBinContent(i) ));
   h_xs_RecoE_selected_abs->Multiply( h_betheMean_muon );
   h_xs_RecoE_selected_abs->Scale( scale_factor );

   TH1D* h_xs_TrueE_selected_abs = (TH1D*) h_selected_pion_incidentTrueE->Clone("h_xs_TrueE_selected_abs");
   h_xs_TrueE_selected_abs->GetXaxis()->SetTitle("True Kinetic Energy (MeV)");
   TH1D* dummy_TrueE_abs = (TH1D*) h_selected_pion_incidentTrueE->Clone("h_xs_TrueE_selected_abs");
   
   if(isMC){
   dummy_TrueE_abs->Add( h_selected_abs_interactingTrueE, -1);
   h_xs_TrueE_selected_abs->Divide( dummy_TrueE_abs );
   for(int i = 1; i <= nBin_int; i++) h_xs_TrueE_selected_abs->SetBinContent(i, log( h_xs_TrueE_selected_abs->GetBinContent(i) ));
   h_xs_TrueE_selected_abs->Multiply( h_betheMean_muon );
   h_xs_TrueE_selected_abs->Scale( scale_factor );
   }
   /* Less Accurate when Nint !<< Ninc
      TH1D* h_xs_RecoE_selected_abs = (TH1D*) h_selected_abs_interactingRecoE->Clone("h_xs_RecoE_selected_abs");
      h_xs_RecoE_selected_abs->Divide( h_selected_pion_incidentRecoE );
      h_xs_RecoE_selected_abs->Multiply( h_betheMean_muon );
      h_xs_RecoE_selected_abs->Scale( factor_mbarn*scale_factor );
      */

   //Try Bin Errors from error propagation derivation of sigma(E)
   //error is like error of binomial factor*sqrt( p(1-p)/Ninc) where p=Nint/Ninc

   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_abs_interactingRecoE->GetBinContent(i) / h_selected_pion_incidentRecoE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentRecoE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_muon->GetBinContent(i);

      h_xs_RecoE_selected_abs->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_RecoE_selected_abs->Write();

   if(isMC){
   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_abs_interactingTrueE->GetBinContent(i) / h_selected_pion_incidentTrueE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentTrueE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_muon->GetBinContent(i);

      h_xs_TrueE_selected_abs->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_TrueE_selected_abs->Write();
   }

   //------------------------------------------------------

   //=====================================================*
   //            Plotting and Style
   //=====================================================
   //
   //h_xs_RecoE_selected_abs->Scale(2);

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

   /*
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

