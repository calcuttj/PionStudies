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
   ROOT::RDataFrame frame(inputTree, mcFilepath);

   //access Jakes GeantFile in folder
   TFile f1("exclusive_xsec.root");
   TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
   TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
   TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
   f1.Close();

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(totInel_KE);
   mg->Add(abs_KE);
   mg->Add(cex_KE);
   mg->SetTitle("Cross-Section; true kinetic Energy (MeV); #sigma (mbarn)");
   
   //output file
   TFile *output = new TFile ("output_eSliceMethod_trueProcess_recoEnergy.root", "RECREATE");

   double bin_size = 40.;
   double eStart = 1200.;
   double eEnd = 0.;
   int nBin = (eStart - eEnd) / bin_size;
   double mass_pion = 139.; // MeV
   
   //Incident Histogram and Interacting Histogram, interacting for different Processes
   //Energies are kinetic Energy of beam particle

   TH1D* h_true_pion_reco_incidentE = new TH1D("h_true_pion_reco_incidentE", "Incident True Pion", nBin, eEnd, eStart);
   h_true_pion_reco_incidentE->GetXaxis()->SetTitle("Reco KE (MeV)");

   TH1D* h_true_abs_reco_interactingE = new TH1D("h_true_abs_reco_interactingE", "Interacting true ABS Energy", nBin, eEnd, eStart);
   h_true_abs_reco_interactingE->GetXaxis()->SetTitle("Reco KE (MeV)");
   TH1D* h_true_cex_reco_interactingE = new TH1D("h_true_cex_true_interactingE", "Interacting true CEX Energy", nBin, eEnd, eStart);
   h_true_cex_reco_interactingE->GetXaxis()->SetTitle("Reco KE (MeV)");
   TH1D* h_true_totInel_reco_interactingE = new TH1D("h_true_totInel_reco_interactingE", "Interacting true total INEL Energy", nBin, eEnd, eStart);
   h_true_totInel_reco_interactingE->GetXaxis()->SetTitle("Reco KE (MeV)");

   //Initial Filters for all events
   auto frame_filter = frame
      .Filter("true_beam_endZ > 0")
      .Filter("primary_isBeamType && passBeamCut");    //particles that make it into TPC, couldn't be considered otherwise


   //Filter for different interaction types
   //
   //******************************************************
   //          TRUE Interactions, RECO Energies
   //******************************************************
   //------------------------------------------------------
   //Incident true sample, reconstructed Energy
   //------------------------------------------------------
   //
   //

   auto mcIncident_true_primaryPi = frame_filter      
      .Filter("true_beam_PDG == 211")
      //.Filter("true_primPionInel && !isDecay")
      .Define("reco_firstEntryIncident", firstIncident, {"reco_beam_incidentEnergies"})
      //.Define("reco_interactingKE", reco_interactingE, {"reco_beam_incidentEnergies"})
      .Define("reco_interactingKE", reco_interactingE, {"reco_beam_incidentEnergies", "reco_beam_calibrated_dEdX", "reco_beam_TrkPitch", "reco_beam_interactingEnergy"})
      .Filter("reco_interactingKE > 0 ");
   //.Range(20);


   //Build the Incident Histogram
   //---------
   //take the interacting energy, that slice is filled too
   //first entry of true_beam_incidentEnergies is the first bin to be filled up until the interacting energy bin

   mcIncident_true_primaryPi
      .Foreach( [h_true_pion_reco_incidentE, &bin_size, &nBin] (double reco_firstEntryIncident, double reco_beam_interactingEnergy) {
            //bin that energy falls into is (int) energy/nbins + 1
            int binNumber_initEnergy = (int) reco_firstEntryIncident / bin_size + 1;
            int binNumber_interEnergy = (int) reco_beam_interactingEnergy / bin_size + 1;
            //how many bins need to be filled
            //starting at interacting energy (lower than first incident energy)
            int cnt = 0;

            //if(binNumber_initEnergy < 0 || binNumber_interEnergy < 0) return;
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      


              if( i > 0 && i <= nBin){ //make sure we don't go outside of bin range
                  h_true_pion_reco_incidentE->AddBinContent(i);
               };
            //h_true_pion_reco_incidentE->Fill( true_firstEntryIncident - cnt*bin_size);
               cnt++;
            };
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});
            //,{"true_firstEntryIncident", "true_beam_interactingEnergy"});

   h_true_pion_reco_incidentE->Sumw2(0);
   h_true_pion_reco_incidentE->Write();


   //******************************************************
   //------------------------------------------------------
   //Interacting true samples
   //------------------------------------------------------
   //
   //interacting samples are considered within the first APA, include APA3Cut bc of bad Energy reco after gap
   // --> For incident we still need to take into account that some interact 
   // further back in the TPC, so the cut on APA3 is not applied for the incident sample
   //
   auto mcInteracting_true_allPrimaryPi = mcIncident_true_primaryPi.Filter("primary_ends_inAPA3");

   //for the true samples no beam Cuts are applied as one should reproduce the geant XS

   auto mcInteracting_true_abs = mcInteracting_true_allPrimaryPi
      .Filter("true_absSignal && true_pion_daughter == 0");
      //.Filter("true_absSignal");

   mcInteracting_true_abs
      .Foreach( [h_true_abs_reco_interactingE] (double true_beam_interactingEnergy){
            h_true_abs_reco_interactingE->Fill(true_beam_interactingEnergy);}
            ,{"reco_interactingKE"});
            //,{"true_beam_interactingEnergy"});

   h_true_abs_reco_interactingE->Sumw2(0);
   h_true_abs_reco_interactingE->Write();

   auto mcInteracting_true_cex = mcInteracting_true_allPrimaryPi
      .Filter("true_chexSignal && true_pion_daughter == 0");
      //.Filter("true_chexSignal");

   mcInteracting_true_cex
      .Foreach( [h_true_cex_reco_interactingE] (double true_beam_interactingEnergy){
            h_true_cex_reco_interactingE->Fill(true_beam_interactingEnergy);}
            ,{"reco_interactingKE"});
            //,{"true_beam_interactingEnergy"});

   h_true_cex_reco_interactingE->Sumw2(0);
   h_true_cex_reco_interactingE->Write();

   auto mcInteracting_true_totInel = mcInteracting_true_allPrimaryPi;

   mcInteracting_true_totInel
      .Foreach( [h_true_totInel_reco_interactingE] (double true_beam_interactingEnergy){
            h_true_totInel_reco_interactingE->Fill(true_beam_interactingEnergy);}
            ,{"reco_interactingKE"});
            //,{"true_beam_interactingEnergy"});

   h_true_totInel_reco_interactingE->Sumw2(0);
   h_true_totInel_reco_interactingE->Write();

   //******************************************************
   //            Prepare BetheBloch Mean for each Bin 
   //******************************************************
   TH1D* h_betheMean_pion = new TH1D("h_betheMean_pion", "Mean Energy Loss", nBin, eEnd, eStart);

   //fill histo with Mean dEdX of bin center
   for(int i = 1; i <= nBin; i++){
      h_betheMean_pion->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size  , mass_pion) );
   };

   h_betheMean_pion->Write();

   //******************************************************
   //             Computing the XS
   //******************************************************
   //
   // xs(Ebin) = (A / (Na*density*bin_size)) * dEdX(Ebin) * hInteracting / hIncident
   //
   // More Accurate use log( Ninc / (Ninc - Nint ))
   //
   double scale_factor = atomic_mass / ( density * N_avogadro * bin_size );
   
   //------------------------------------------------------
   //    Absorption, True Interactions, True Energy
   //------------------------------------------------------
   //
   

   TH1D* h_xs_true_abs_reco_E = (TH1D*) h_true_pion_reco_incidentE->Clone("h_xs_true_abs_reco_E");
   TH1D* dummy_abs = (TH1D*) h_true_pion_reco_incidentE->Clone("h_xs_true_abs_reco_E");
   dummy_abs->Add( h_true_abs_reco_interactingE, -1);
   h_xs_true_abs_reco_E->Divide( dummy_abs );
   for(int i = 1; i <= nBin; i++) h_xs_true_abs_reco_E->SetBinContent(i, log( h_xs_true_abs_reco_E->GetBinContent(i) ));
   h_xs_true_abs_reco_E->Multiply( h_betheMean_pion );
   h_xs_true_abs_reco_E->Scale( factor_mbarn*scale_factor );
   
   /* Less Accurate when Nint !<< Ninc
   TH1D* h_xs_true_abs_reco_E = (TH1D*) h_true_abs_reco_interactingE->Clone("h_xs_true_abs_reco_E");
   h_xs_true_abs_reco_E->Divide( h_true_pion_reco_incidentE );
   h_xs_true_abs_reco_E->Multiply( h_betheMean_pion );
   h_xs_true_abs_reco_E->Scale( factor_mbarn*scale_factor );
   */

   //Try Bin Errors from error propagation derivation of sigma(E)
   //error is like error of binomial factor*sqrt( p(1-p)/Ninc) where p=Nint/Ninc
   
   for(int i=1; i <= nBin; i++){

      double p = h_true_abs_reco_interactingE->GetBinContent(i) / h_true_pion_reco_incidentE->GetBinContent(i);
      double nInc_i = h_true_pion_reco_incidentE->GetBinContent(i);
      double help_factor = scale_factor*factor_mbarn*h_betheMean_pion->GetBinContent(i);
      
      h_xs_true_abs_reco_E->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_true_abs_reco_E->Write();


   //------------------------------------------------------
   //    ChargeExchange, True Interactions, True Energy
   //------------------------------------------------------
   //
   
   TH1D* h_xs_true_cex_reco_E = (TH1D*) h_true_pion_reco_incidentE->Clone("h_xs_true_cex_reco_E");
   TH1D* dummy_cex = (TH1D*) h_true_pion_reco_incidentE->Clone("h_xs_true_cex_reco_E");
   dummy_cex->Add( h_true_cex_reco_interactingE, -1);
   h_xs_true_cex_reco_E->Divide( dummy_cex );
   for(int i = 1; i <= nBin; i++) h_xs_true_cex_reco_E->SetBinContent(i, log( h_xs_true_cex_reco_E->GetBinContent(i) ));
   h_xs_true_cex_reco_E->Multiply( h_betheMean_pion );
   h_xs_true_cex_reco_E->Scale( factor_mbarn*scale_factor );
   
   /* Less Accurate when Nint !<< Ninc
   TH1D* h_xs_true_cex_reco_E = (TH1D*) h_true_cex_reco_interactingE->Clone("h_xs_true_cex_reco_E");
   h_xs_true_cex_reco_E->Divide( h_true_pion_reco_incidentE );
   h_xs_true_cex_reco_E->Multiply( h_betheMean_pion );
   h_xs_true_cex_reco_E->Scale( factor_mbarn*scale_factor );
   */

  //Bin Errors
   for(int i=1; i <= nBin; i++){

      double p = h_true_cex_reco_interactingE->GetBinContent(i) / h_true_pion_reco_incidentE->GetBinContent(i);
      double nInc_i = h_true_pion_reco_incidentE->GetBinContent(i);
      double help_factor = scale_factor*factor_mbarn*h_betheMean_pion->GetBinContent(i);
      
      h_xs_true_cex_reco_E->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_true_cex_reco_E->Write();
   //------------------------------------------------------
   //    Inelastic, True Interactions, True Energy
   //------------------------------------------------------
   //
   
  TH1D* h_xs_true_totInel_reco_E = (TH1D*) h_true_pion_reco_incidentE->Clone("h_xs_true_totInel_reco_E");
   TH1D* dummy_totInel = (TH1D*) h_true_pion_reco_incidentE->Clone("h_xs_true_totInel_reco_E");
   dummy_totInel->Add( h_true_totInel_reco_interactingE, -1);
   h_xs_true_totInel_reco_E->Divide( dummy_totInel );
   for(int i = 1; i <= nBin; i++) h_xs_true_totInel_reco_E->SetBinContent(i, log( h_xs_true_totInel_reco_E->GetBinContent(i) ));
   h_xs_true_totInel_reco_E->Multiply( h_betheMean_pion );
   h_xs_true_totInel_reco_E->Scale( factor_mbarn*scale_factor );
   
   /* Less Accurate when Nint !<< Ninc
   TH1D* h_xs_true_totInel_reco_E = (TH1D*) h_true_totInel_reco_interactingE->Clone("h_xs_true_totInel_reco_E");
   h_xs_true_totInel_reco_E->Divide( h_true_pion_reco_incidentE );
   h_xs_true_totInel_reco_E->Multiply( h_betheMean_pion );
   h_xs_true_totInel_reco_E->Scale( factor_mbarn*scale_factor );
   */

   //Bin Errors
   for(int i=1; i <= nBin; i++){

      double p = h_true_totInel_reco_interactingE->GetBinContent(i) / h_true_pion_reco_incidentE->GetBinContent(i);
      double nInc_i = h_true_pion_reco_incidentE->GetBinContent(i);
      double help_factor = scale_factor*factor_mbarn*h_betheMean_pion->GetBinContent(i);
      
      h_xs_true_totInel_reco_E->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_true_totInel_reco_E->Write();
 

   //******************************************************
   //            Plotting and Style
   //******************************************************
   //
   
   TCanvas *c_abs = new TCanvas("c_abs", "c_abs");
   gPad->SetGrid(1,1);
   h_xs_true_abs_reco_E->SetTitle( "Absorption;Kinetic Energy (MeV); #sigma (mb)");
   h_xs_true_abs_reco_E->GetXaxis()->SetRangeUser(400,900);
   h_xs_true_abs_reco_E->GetXaxis()->SetNdivisions(1020);
   h_xs_true_abs_reco_E->GetYaxis()->SetNdivisions(1020);
   
   abs_KE->SetTitle( "Absorption;Kinetic Energy (MeV); #sigma (mb)");
   abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   abs_KE->SetLineColor(kRed);
   abs_KE->SetLineWidth(3);
   abs_KE->Draw("AC");
   h_xs_true_abs_reco_E->SetMarkerSize(0.7);
   h_xs_true_abs_reco_E->Draw("PE0 SAME");

   c_abs->Write();
   
   TCanvas *c_cex = new TCanvas("c_cex", "c_cex");
   gPad->SetGrid(1,1);
   h_xs_true_cex_reco_E->SetTitle( "Charge Exchange;Kinetic Energy (MeV); #sigma (mb)");
   h_xs_true_cex_reco_E->GetXaxis()->SetRangeUser(400,900);
   h_xs_true_cex_reco_E->GetXaxis()->SetNdivisions(1020);
   h_xs_true_cex_reco_E->GetYaxis()->SetNdivisions(1020);
   
   cex_KE->SetTitle( "Charge Exchange;Kinetic Energy (MeV); #sigma (mb)");
   cex_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   cex_KE->GetYaxis()->SetRangeUser(0, 200);
   cex_KE->SetLineColor(kRed);
   cex_KE->SetLineWidth(3);
   cex_KE->Draw("AC"); 
   h_xs_true_cex_reco_E->SetMarkerSize(0.7);
   h_xs_true_cex_reco_E->Draw("PE0 SAME");

   c_cex->Write();

   TCanvas *c_totInel = new TCanvas("c_totInel", "c_totInel");
   gPad->SetGrid(1,1);
   h_xs_true_totInel_reco_E->SetTitle( "Total Inelastic;Kinetic Energy (MeV); #sigma (mb)");
   h_xs_true_totInel_reco_E->GetXaxis()->SetRangeUser(400,900);
   h_xs_true_totInel_reco_E->GetXaxis()->SetNdivisions(1020);
   h_xs_true_totInel_reco_E->GetYaxis()->SetNdivisions(1020);
   
   totInel_KE->SetTitle( "Total Inelastic;Kinetic Energy (MeV); #sigma (mb)");
   totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   totInel_KE->SetLineColor(kRed);
   totInel_KE->SetLineWidth(3);
   totInel_KE->Draw("AC");
   h_xs_true_totInel_reco_E->SetMarkerSize(0.7);
   h_xs_true_totInel_reco_E->Draw("PE0 SAME");

   c_totInel->Write();

   TCanvas *c_all = new TCanvas("c_all", "c_all");
   gPad->SetGrid(1,1);
   
   cex_KE->SetLineColor(kGreen);
   abs_KE->SetLineColor(kBlue);
   mg->GetXaxis()->SetRangeUser(0,1100);
   mg->GetXaxis()->SetNdivisions(1020);
   mg->Draw("AC");
   h_xs_true_totInel_reco_E->Draw("PE0 SAME");
   h_xs_true_cex_reco_E->Draw("PE0 SAME");
   h_xs_true_abs_reco_E->Draw("PE0 SAME");

   c_all->Write();




   output->Write();
   //f1.Close();
   return 0;
}

