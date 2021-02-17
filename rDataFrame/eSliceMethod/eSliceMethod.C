#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
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

int eSliceMethod(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(inputTree, mcFilepath);

   //access Jakes GeantFile in folder
   TFile f1("exclusive_xsec.root");
   TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
   TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
   TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
   f1.Close();
   
   //output file
   TFile *output = new TFile ("output_eSliceMethod.root", "RECREATE");

   double bin_size = 20.;
   double eStart = 1000.;
   double eEnd = 0.;
   int nBin = (eStart - eEnd) / bin_size;

   //Incident Histogram and Interacting Histogram, interacting for different Processes
   //Energies are kinetic Energy of beam particle

   TH1D* h_true_pion_true_incidentE = new TH1D("h_true_pion_true_incidentE", "incident true Energy", nBin, eEnd, eStart);

   TH1D* h_true_abs_true_interactingE = new TH1D("h_true_abs_true_interactingE", "interacting true ABS Energy", nBin, eEnd, eStart);
   TH1D* h_true_cex_true_interactingE = new TH1D("h_true_cex_true_interactingE", "interacting true CEX Energy", nBin, eEnd, eStart);
   TH1D* h_true_totInel_true_interactingE = new TH1D("h_true_totInel_true_interactingE", "interacting true total INEL Energy", nBin, eEnd, eStart);

   //Initial Filters for all events
   auto frame_filter = frame
      .Filter("true_beam_endZ > 0")
      .Filter("primary_isBeamType && passBeamCut");    //particles that make it into TPC, couldn't be considered otherwise


   //Filter for different interaction types
   //
   //******************************************************
   //          TRUE Interactions, TRUE Energies
   //******************************************************
   //------------------------------------------------------
   //Incident true sample
   //------------------------------------------------------
   //
   //for the true samples no beam Cuts are applied as one should reproduce the geant XS

   auto mcIncident_true_primaryPi = frame_filter      
      .Filter("true_primPionInel && !isDecay")
      .Define("true_firstEntryIncident", firstIncident, {"true_beam_incidentEnergies"});
   //.Range(50);

   //Build the Incident Histogram
   //---------
   //take the interacting energy, that slice is filled too
   //first entry of true_beam_incidentEnergies is the first bin to be filled up until the interacting energy bin

   mcIncident_true_primaryPi
      .Foreach( [h_true_pion_true_incidentE, &bin_size, &nBin] (double true_firstEntryIncident, double true_beam_interactingEnergy) {
            //bin that energy falls into is (int) energy/nbins + 1
            int binNumber_initEnergy = (int) true_firstEntryIncident / bin_size + 1;
            int binNumber_interEnergy = (int) true_beam_interactingEnergy / bin_size + 1;
            //how many bins need to be filled
            //start from initial energy value and subtract nBin at every time
            int cnt = 0;
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){ 
               
              if(i <= nBin){
                  h_true_pion_true_incidentE->AddBinContent(i);
               };
            //h_true_pion_true_incidentE->Fill( true_firstEntryIncident - cnt*bin_size);
               cnt++;
            };
            }
            ,{"true_firstEntryIncident", "true_beam_interactingEnergy"});

   h_true_pion_true_incidentE->Sumw2(0);
   h_true_pion_true_incidentE->Write();



   //******************************************************
   //------------------------------------------------------
   //Interacting true samples
   //------------------------------------------------------
   //
   //interacting samples are considered within the first APA, include Cut to show where our measurement can measure XS
   // --> For incident we still need to take into account that some interact 
   // further back in the TPC, so the cut on APA3 is not applied for the incident sample
   //
   auto mcInteracting_true_allPrimaryPi = mcIncident_true_primaryPi;//.Filter("true_beam_endZ < 226");

   //for the true samples no beam Cuts are applied as one should reproduce the geant XS

   auto mcInteracting_true_abs = mcInteracting_true_allPrimaryPi
      .Filter("true_absSignal && true_pion_daughter == 0");
      //.Filter("true_absSignal");

   mcInteracting_true_abs
      .Foreach( [h_true_abs_true_interactingE] (double true_beam_interactingEnergy){
            h_true_abs_true_interactingE->Fill(true_beam_interactingEnergy);}
            ,{"true_beam_interactingEnergy"});

   h_true_abs_true_interactingE->Sumw2(0);
   h_true_abs_true_interactingE->Write();

   auto mcInteracting_true_cex = mcInteracting_true_allPrimaryPi
      .Filter("true_chexSignal && true_pion_daughter == 0");
      //.Filter("true_chexSignal");

   mcInteracting_true_cex
      .Foreach( [h_true_cex_true_interactingE] (double true_beam_interactingEnergy){
            h_true_cex_true_interactingE->Fill(true_beam_interactingEnergy);}
            ,{"true_beam_interactingEnergy"});

   h_true_cex_true_interactingE->Sumw2(0);
   h_true_cex_true_interactingE->Write();

   auto mcInteracting_true_totInel = mcInteracting_true_allPrimaryPi;

   mcInteracting_true_totInel
      .Foreach( [h_true_totInel_true_interactingE] (double true_beam_interactingEnergy){
            h_true_totInel_true_interactingE->Fill(true_beam_interactingEnergy);}
            ,{"true_beam_interactingEnergy"});

   h_true_totInel_true_interactingE->Sumw2(0);
   h_true_totInel_true_interactingE->Write();

   //******************************************************
   //            Prepare BetheBloch Mean for each Bin 
   //******************************************************
   TH1D* h_betheMean_pion = new TH1D("h_betheMean_pion", "Mean Energy Loss", nBin, eEnd, eStart);

   //fill histo with Mean dEdX of bin center
   for(int i = 1; i <= nBin; i++){
      h_betheMean_pion->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size  , mass_pion) );
   };

   h_betheMean_pion->Write();

   /*TH1D* h_test3 = (TH1D*) h_betheMean_pion->Clone("bethe/20");
   h_test3->Divide(h_test3);
   h_test3->Sumw2(0);
   h_test3->Divide(h_betheMean_pion);
   h_test3->Sumw2(0);
   h_test3->Scale(20);
   h_test3->Sumw2(0);
   h_test3->Write();*/
   //******************************************************
   //             Computing the XS
   //******************************************************
   //
   // xs(Ebin) = (A / (Na*density*bin_size)) * dEdX(Ebin) * hInteracting / hIncident
   //
   double scale_factor = atomic_mass / ( density * N_avogadro * bin_size );
   
   //------------------------------------------------------
   //    Absorption, True Interactions, True Energy
   //------------------------------------------------------
   //
   

   TH1D* h_xs_true_abs_true_E = (TH1D*) h_true_abs_true_interactingE->Clone("h_xs_true_abs_true_E");
   h_xs_true_abs_true_E->Divide( h_true_pion_true_incidentE );
   h_xs_true_abs_true_E->Multiply( h_betheMean_pion );
   h_xs_true_abs_true_E->Sumw2(0);
   h_xs_true_abs_true_E->Scale( factor_mbarn*scale_factor );
   h_xs_true_abs_true_E->Sumw2(0);
   h_xs_true_abs_true_E->Write();


   //------------------------------------------------------
   //    ChargeExchange, True Interactions, True Energy
   //------------------------------------------------------
   //
   
   TH1D* h_xs_true_cex_true_E = (TH1D*) h_true_cex_true_interactingE->Clone("h_xs_true_cex_true_E");
   h_xs_true_cex_true_E->Divide( h_true_pion_true_incidentE );
   h_xs_true_cex_true_E->Multiply( h_betheMean_pion );
   h_xs_true_cex_true_E->Sumw2(0);
   h_xs_true_cex_true_E->Scale( factor_mbarn*scale_factor );
   h_xs_true_cex_true_E->Sumw2(0);
   h_xs_true_cex_true_E->Write();
 
   //------------------------------------------------------
   //    Inelastic, True Interactions, True Energy
   //------------------------------------------------------
   //
   
   TH1D* h_xs_true_totInel_true_E = (TH1D*) h_true_totInel_true_interactingE->Clone("h_xs_true_totInel_true_E");
   h_xs_true_totInel_true_E->Divide( h_true_pion_true_incidentE );
   h_xs_true_totInel_true_E->Multiply( h_betheMean_pion );
   h_xs_true_totInel_true_E->Sumw2(0);
   h_xs_true_totInel_true_E->Scale( factor_mbarn*scale_factor );
   h_xs_true_totInel_true_E->Sumw2(0);
   h_xs_true_totInel_true_E->Write();
 

   //******************************************************
   //            Plotting and Style
   //******************************************************
   //
   
   TCanvas *c_abs = new TCanvas("c_abs", "c_abs");
   h_xs_true_abs_true_E->SetTitle( "Absorption;KE (MeV); #sigma (mb)");
   h_xs_true_abs_true_E->GetXaxis()->SetRangeUser(400,900);
   
   abs_KE->SetTitle( "Absorption;KE (MeV); #sigma (mb)");
   abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   abs_KE->SetLineColor(kRed);
   abs_KE->SetLineWidth(2);
   abs_KE->Draw("AC");
   h_xs_true_abs_true_E->SetMarkerColor(kGreen);
   h_xs_true_abs_true_E->Draw("PE0 SAME");

   c_abs->Write();
   
   TCanvas *c_cex = new TCanvas("c_cex", "c_cex");
   h_xs_true_cex_true_E->SetTitle( "Charge Exchange;KE (MeV); #sigma (mb)");
   h_xs_true_cex_true_E->GetXaxis()->SetRangeUser(400,900);
   
   cex_KE->SetTitle( "Charge Exchange;KE (MeV); #sigma (mb)");
   cex_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   cex_KE->SetLineColor(kRed);
   cex_KE->SetLineWidth(2);
   cex_KE->Draw("AC"); 
   h_xs_true_cex_true_E->Draw("PE0 SAME");

   c_cex->Write();

   TCanvas *c_totInel = new TCanvas("c_totInel", "c_totInel");
   h_xs_true_totInel_true_E->SetTitle( "Total Inelastic;KE (MeV); #sigma (mb)");
   h_xs_true_totInel_true_E->GetXaxis()->SetRangeUser(400,900);
   
   totInel_KE->SetTitle( "Total Inelastic;KE (MeV); #sigma (mb)");
   totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   totInel_KE->SetLineColor(kRed);
   totInel_KE->SetLineWidth(2);
   totInel_KE->Draw("AC");
   h_xs_true_totInel_true_E->Draw("PE0 SAME");

   c_totInel->Write();



   output->Write();
   //f1.Close();
   return 0;
}

