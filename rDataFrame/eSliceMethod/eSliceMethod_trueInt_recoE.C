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

int eSliceMethod_trueInt_recoE(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);

   //output file
   string outputName = "output_eSliceMethod_trueInt_recoE_selection.root";
   TFile *output = new TFile ( outputName.c_str() , "RECREATE");

   //access Jakes GeantFile in folder
   TFile f1("exclusive_xsec.root");
   TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
   TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
   TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
   f1.Close();

   //Get Efficiencies and Purities for the selections bin-by-bin
   //TFile f2("eSliceMethod_effPur_binByBin_noSmear.root");
   //TH1D *h_eff_interacting_abs = (TH1D*)f2.Get("h_eff_interacting_abs");
   //TH1D *h_pur_interacting_abs = (TH1D*)f2.Get("h_pur_interacting_abs");
   //TH1D *h_eff_interacting_cex = (TH1D*)f2.Get("h_eff_interacting_cex");
   //TH1D *h_pur_interacting_cex = (TH1D*)f2.Get("h_pur_interacting_cex");
   //TH1D *h_eff_interacting_totInel = (TH1D*)f2.Get("h_eff_interacting_totInel");
   //TH1D *h_pur_interacting_totInel = (TH1D*)f2.Get("h_pur_interacting_totInel");
   //TH1D *h_eff_incidentPion = (TH1D*)f2.Get("h_eff_incidentPion");
   //TH1D *h_pur_incidentPion = (TH1D*)f2.Get("h_pur_incidentPion");


   string mg_title = "Cross-Section MC - true Abs and true Pion Incident; Reco kinetic Energy (MeV); #sigma (mbarn)";
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

   //Initial Filters for all events
   //FILTER for true_primaryPionInelastic
   auto mcIncident_selected_primaryPi = frame
      //.Filter("true_beam_PDG == abs(211)")
      .Filter("primary_isBeamType && passBeamCut && passBeamCutBI")
      .Filter("true_primPionInel")
      .Filter("selected_incidentPion");

   //========================================================
   //Build the Incident Histogram
   //---------
   //take the interacting energy, that slice is filled too
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
   h_selected_pion_incidentRecoE->Sumw2(0);
   h_selected_pion_incidentRecoE->Rebin( bin_size_int/bin_size_inc );
   h_selected_pion_incidentRecoE->Scale( 1 / (bin_size_int/bin_size_inc) );


   //=====================================================
   //------------------------------------------------------
   //Interacting selected samples
   //------------------------------------------------------
   //
   //interacting samples are considered within the first APA, include APA3Cut bc of bad Energy reco after gap
   // --> For incident we still need to take into account that some interact 
   // further back in the TPC, so the cut on APA3 is not applied for the incident sample
   //
   auto mcInteracting_selected_allPrimaryPi = mcIncident_selected_primaryPi;
   //.Filter("primary_ends_inAPA3");


   auto mcInteracting_selected_abs = mcInteracting_selected_allPrimaryPi
      .Filter("true_absSignal")
   .Filter("selected_abs");
   
   mcInteracting_selected_abs 
      .Foreach( [h_selected_abs_interactingRecoE] (double reco_beam_interactingEnergy){
            h_selected_abs_interactingRecoE->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});


   //Efficiencies and Purities Interacting
   //-------------------------------------
   TH1D* h_true_abs = new TH1D("h_true_abs", "Selected Abs, trueAbs", nBin_int, eEnd, eStart);
   TH1D* h_sel_abs = new TH1D("h_sel_abs", "Selected Abs", nBin_int, eEnd, eStart);
   TH1D* h_inc_avail_abs = new TH1D("h_inc_avail_abs", "Available True Abs in Incident", nBin_int, eEnd, eStart);

   TH1D* h_eff_interacting_abs = new TH1D("h_eff_interacting_abs", "Efficiency Interacting Abs Histo", nBin_int, eEnd, eStart);
   TH1D* h_pur_interacting_abs = new TH1D("h_pur_interacting_abs", "Purity Interacting Abs Histo", nBin_int, eEnd, eStart);

   mcIncident_selected_primaryPi
      .Filter("true_absSignal")
      .Foreach( [h_inc_avail_abs] (double reco_beam_interactingEnergy){
            h_inc_avail_abs->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   mcInteracting_selected_abs
      .Foreach( [h_true_abs] (double true_beam_interactingEnergy, double reco_beam_interactingEnergy){

            int binReco = (int) reco_beam_interactingEnergy / bin_size_int + 1;
            int binTrue = (int) true_beam_interactingEnergy / bin_size_int + 1;

            if (binReco == binTrue) h_true_abs->Fill(reco_beam_interactingEnergy); }
            ,{"true_KEint_fromEndP", "reco_interactingKE"});
   
   mcInteracting_selected_abs
      .Foreach( [h_sel_abs] (double reco_beam_interactingEnergy){
            h_sel_abs->Fill(reco_beam_interactingEnergy);}
            , {"reco_interactingKE"} );

   h_eff_interacting_abs->Divide( h_true_abs, h_inc_avail_abs );
   h_pur_interacting_abs->Divide( h_true_abs, h_sel_abs );
   h_eff_interacting_abs->Write();
   h_pur_interacting_abs->Write();

   //Efficiencies and Purities Incident
   //-------------------------------------
   //
   TH1D* h_true_incPi = new TH1D("h_true_incPi", "Selected incPi, trueincPi", nBin_int, eEnd, eStart);
   TH1D* h_sel_incPi = new TH1D("h_sel_incPi", "Selected incPi", nBin_int, eEnd, eStart);
   TH1D* h_inc_avail_incPi = new TH1D("h_inc_avail_incPi", "Available True incPi in Incident", nBin_int, eEnd, eStart);

   TH1D* h_eff_incPi = new TH1D("h_eff_interacting_incPi", "Efficiency Interacting incPi Histo", nBin_int, eEnd, eStart);
   TH1D* h_pur_incPi = new TH1D("h_pur_interacting_incPi", "Purity Interacting incPi Histo", nBin_int, eEnd, eStart);

   frame
      .Filter("primary_isBeamType && passBeamCut && passBeamCutBI")
      .Filter("true_primPionInel")
      .Foreach( [h_inc_avail_incPi] (double reco_beam_interactingEnergy){
            h_inc_avail_incPi->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   mcIncident_selected_primaryPi
      .Foreach( [h_true_incPi] (double true_beam_interactingEnergy, double reco_beam_interactingEnergy){

            int binReco = (int) reco_beam_interactingEnergy / bin_size_int + 1;
            int binTrue = (int) true_beam_interactingEnergy / bin_size_int + 1;

            if (binReco == binTrue) h_true_incPi->Fill(reco_beam_interactingEnergy); }
            ,{"true_KEint_fromEndP", "reco_interactingKE"});
   
   mcIncident_selected_primaryPi
      .Foreach( [h_sel_incPi] (double reco_beam_interactingEnergy){
            h_sel_incPi->Fill(reco_beam_interactingEnergy);}
            , {"reco_interactingKE"} );

   h_eff_incPi->Divide( h_true_incPi, h_inc_avail_incPi );
   h_pur_incPi->Divide( h_true_incPi, h_sel_incPi );
   h_eff_incPi->Write();
   h_pur_incPi->Write();

 

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

   h_selected_pion_incidentRecoE->Multiply( h_pur_incPi );
   h_selected_pion_incidentRecoE->Divide( h_eff_incPi );

   h_selected_abs_interactingRecoE->Multiply( h_pur_interacting_abs );
   h_selected_abs_interactingRecoE->Divide( h_eff_interacting_abs );
 
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
   string xs_abs_name = "h_xs_recoE_true_abs_mc";
   string xs_abs_title = "True Absorption, recoE - MC";


   TH1D* h_xs_RecoE_selected_abs = (TH1D*) h_selected_pion_incidentRecoE->Clone(xs_abs_name.c_str());
   h_xs_RecoE_selected_abs->SetTitle(xs_abs_title.c_str());
   TH1D* dummy_recoE_abs = (TH1D*) h_selected_pion_incidentRecoE->Clone(xs_abs_name.c_str());
   dummy_recoE_abs->Add( h_selected_abs_interactingRecoE, -1);
   h_xs_RecoE_selected_abs->Divide( dummy_recoE_abs );
   for(int i = 1; i <= nBin_int; i++) h_xs_RecoE_selected_abs->SetBinContent(i, log( h_xs_RecoE_selected_abs->GetBinContent(i) ));
   h_xs_RecoE_selected_abs->Multiply( h_betheMean_muon );
   h_xs_RecoE_selected_abs->Scale( scale_factor );

   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_abs_interactingRecoE->GetBinContent(i) / h_selected_pion_incidentRecoE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentRecoE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_muon->GetBinContent(i);

      h_xs_RecoE_selected_abs->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_RecoE_selected_abs->Write();


   //=====================================================*
   //            Plotting and Style
   //=====================================================
   //
   //h_xs_RecoE_selected_abs->Scale(2);
   //h_xs_RecoE_selected_cex->Scale(2);
   //h_xs_RecoE_selected_totInel->Scale(2);

   TCanvas *c_RecoE_abs = new TCanvas("c_RecoE_abs", "c_RecoE_abs");
   gPad->SetGrid(1,1);
   h_xs_RecoE_selected_abs->SetTitle( "True Absorption; Reco Kinetic Energy (MeV); #sigma (mb)");
   h_xs_RecoE_selected_abs->GetXaxis()->SetRangeUser(400,900);
   h_xs_RecoE_selected_abs->GetXaxis()->SetNdivisions(1020);
   h_xs_RecoE_selected_abs->GetYaxis()->SetNdivisions(1020);

   string c_abs_title;
   c_abs_title = "True Absorption MC; Reco kinetic Energy (MeV); #sigma (mbarn)";

   abs_KE->SetTitle( c_abs_title.c_str());
   abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   abs_KE->SetLineColor(kBlue);
   abs_KE->SetLineWidth(3);
   abs_KE->Draw("AC");
   h_xs_RecoE_selected_abs->SetMarkerSize(0.7);
   h_xs_RecoE_selected_abs->Draw("PE0 SAME");

   c_RecoE_abs->Write();

   //output->Write();
   //f1.Close();
   return 0;
}

