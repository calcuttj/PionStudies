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

int eSliceMethod_selectedInt_trueInt(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(inputTree, mcFilepath);
   gStyle->SetNdivisions(1020);

   //output file
   TFile *output = new TFile ("output_eSliceMethod_selectedEvents_onlyTrueInt.root", "RECREATE");

   //access Jakes GeantFile in folder
   TFile f1("exclusive_xsec.root");
   TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
   TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
   TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
   f1.Close();

   TFile f2("fit_mc_Prod4_dEdX_pitch_1_14_21.root");
   TH1D *fit_dEdX_lifetime_mpv = (TH1D*)f2.Get("dEdX_mpv_lifetime"); //mean value corrected for lifetime
   TH1D *fit_pitch_mean = (TH1D*)f2.Get("fit_mc_pitch_mean");

   //TH1D *fit_dEdX_lifetime_mpv = (TH1D*)f2.Get("fit_mc_dEdX_SCEcorr_mpv"); //mean value corrected for lifetime
   //TH1D *fit_pitch_mean = (TH1D*)f2.Get("fit_mc_pitch_SCEcorr_mean");

   //Datadriven corrected pitch
   //TH1D *fit_pitch_mean = (TH1D*)f2.Get("true_pitch_mean_lifetime");

   output->cd();
   fit_dEdX_lifetime_mpv->Write();
   fit_pitch_mean->Write();

   TMultiGraph *mg = new TMultiGraph();
   //mg->Add(totInel_KE);
   mg->Add(abs_KE);
   mg->Add(cex_KE);
   mg->SetTitle("Cross-Section; reco kinetic Energy (MeV); #sigma (mbarn)");

   //--------------------------------------------------------
   //Bethe Bloch MPV and Mean
   //--------------------------------------------------------
   TH1D* bethe_pi_mpv = new TH1D("betheMPV_pion", "", fit_pitch_mean->GetNbinsX(), 1, fit_pitch_mean->GetNbinsX() );
   bethe_pi_mpv->GetXaxis()->SetTitle("wire");  bethe_pi_mpv->GetYaxis()->SetTitle("dEdX (MeV/cm)");

   TH1D* bethe_mu_mpv = new TH1D("betheMPV_muon", "", fit_pitch_mean->GetNbinsX(), 1, fit_pitch_mean->GetNbinsX() );
   bethe_mu_mpv->GetXaxis()->SetTitle("wire");  bethe_mu_mpv->GetYaxis()->SetTitle("dEdX (MeV/cm)");

   TH1D* bethe_pi_mean = new TH1D("betheMean_pion", "", fit_pitch_mean->GetNbinsX(), 1, fit_pitch_mean->GetNbinsX() );
   bethe_pi_mean->GetXaxis()->SetTitle("wire"); bethe_pi_mean->GetYaxis()->SetTitle("dEdX (MeV/cm)");

   TH1D* bethe_mu_mean = new TH1D("betheMean_muon", "", fit_pitch_mean->GetNbinsX(), 1, fit_pitch_mean->GetNbinsX() );
   bethe_mu_mean->GetXaxis()->SetTitle("wire"); bethe_mu_mean->GetYaxis()->SetTitle("dEdX (MeV/cm)");

   hist_bethe_mpv( KE_in_pion, mass_pion, fit_pitch_mean, bethe_pi_mpv);
   hist_bethe_mpv( KE_in_muon, mass_muon, fit_pitch_mean, bethe_mu_mpv);
   hist_bethe_mean( KE_in_pion, mass_pion, fit_pitch_mean, bethe_pi_mean);
   hist_bethe_mean( KE_in_muon, mass_muon, fit_pitch_mean, bethe_mu_mean);

   bethe_pi_mpv->Write();
   bethe_mu_mpv->Write();

   bethe_pi_mean->Write();
   bethe_mu_mean->Write();

   //--------------------------------------------------------
   //PREP for INTERACTING ENERGY
   //Histogram with energy deposit along wire and running sum
   //
   //Strategy for MPV --> Mean of Data is to multiply fitted MPV (stable fit) by factor of betheMean/betheMPV
   TH1D* bethe_frac_mu = new TH1D("betheFrac_mu", "", fit_pitch_mean->GetNbinsX(), 1, fit_pitch_mean->GetNbinsX() );
   bethe_frac_mu->Divide(bethe_mu_mean, bethe_mu_mpv);
   bethe_frac_mu->Write();

   TH1D* dEdX_mean_calc_fit_bethe = (TH1D*)bethe_frac_mu->Clone("dEdX_mean_calc_fit_bethe");
   dEdX_mean_calc_fit_bethe->Multiply(fit_dEdX_lifetime_mpv);
   dEdX_mean_calc_fit_bethe->Write();

   TH1D dE_product_fit_dEdX_pitch = (*dEdX_mean_calc_fit_bethe) * (*fit_pitch_mean);
   dE_product_fit_dEdX_pitch.SetName("dE_product_fit_dEdX_pitch");

   TH1D *runningSum_dE = new TH1D("runningSum_dE", "", fit_pitch_mean->GetNbinsX(), 1, fit_pitch_mean->GetNbinsX() );

   double temp = 0;
   for(int i=0; i <= dE_product_fit_dEdX_pitch.GetNbinsX(); i++){
      temp += dE_product_fit_dEdX_pitch.GetBinContent(i);
      runningSum_dE->SetBinContent( i, temp);
   };

   dE_product_fit_dEdX_pitch.Write();
   runningSum_dE->Write();
   //--------------------------------------------------------

   //Incident Histogram and Interacting Histogram, interacting for different Processes
   //Energies are kinetic Energy of beam particle

   TH1D* h_selected_pion_incidentE = new TH1D("h_selected_pion_incidentE", "Incident Selected Pion", nBin_inc, eEnd, eStart);
   h_selected_pion_incidentE->GetXaxis()->SetTitle("Reco KE (MeV)");

   TH1D* h_selected_abs_interactingE = new TH1D("h_selected_abs_interactingE", "Interacting Selected ABS Energy", nBin_int, eEnd, eStart);
   h_selected_abs_interactingE->GetXaxis()->SetTitle("Reco KE (MeV)");
   TH1D* h_selected_cex_interactingE = new TH1D("h_selected_cex_interactingE", "Interacting Selected CEX Energy", nBin_int, eEnd, eStart);
   h_selected_cex_interactingE->GetXaxis()->SetTitle("Reco KE (MeV)");
   TH1D* h_selected_totInel_interactingE = new TH1D("h_selected_totInel_interactingE", "Interacting Selected total INEL Energy", nBin_int, eEnd, eStart);
   h_selected_totInel_interactingE->GetXaxis()->SetTitle("Reco KE (MeV)");

   //Initial Filters for all events
   auto frame_filter = frame
      .Filter("primary_isBeamType && passBeamCut && passBeamCutBI");    //particles that make it into TPC, couldn't be considered otherwise


   //Filter for different interaction types
   //
   //=====================================================
   //          SELECTED Interactions, RECO Energies
   //=====================================================
   //------------------------------------------------------
   //Incident selected sample, reconstructed Energy
   //------------------------------------------------------
   //
   //for Interacting energy, the wire the primary particle interacted at should be the last one in reco_beam_calo_wire
   //

   auto mcIncident_selected_primaryPi = frame_filter      
      //.Range(20)
      .Define("reco_firstEntryIncident", firstIncident, {"reco_beam_incidentEnergies"})

      .Define("reco_interactingKE", [runningSum_dE](const std::vector<double> &reco_beam_calo_wire, double incidentE){
            double interactingWire = reco_beam_calo_wire[ reco_beam_calo_wire.size() ];
            double interactingKE;
            if(interactingWire >= 1 && interactingWire < runningSum_dE->GetNbinsX()){
            interactingKE = incidentE - runningSum_dE->GetBinContent(interactingWire);
            }
            else interactingKE = -999;
            return interactingKE;
            }
            ,{"reco_beam_calo_wire", "reco_firstEntryIncident"})

      .Define("reco_incident_wire", [](std::vector<double> &reco_beam_calo_wire){
            return reco_beam_calo_wire[ reco_beam_calo_wire[0] ];
            },{"reco_beam_calo_wire"})

   .Define("reco_interacting_wire", [](std::vector<double> &reco_beam_calo_wire){
         return reco_beam_calo_wire[ reco_beam_calo_wire.size() ];
         },{"reco_beam_calo_wire"})

   .Filter("reco_interactingKE > 0 ")
      .Filter("true_beam_PDG == 211"); //select only true pions form event Selection

   //Build the Incident Histogram
   //---------
   //take the interacting energy, that slice is filled too

   mcIncident_selected_primaryPi
      .Foreach( [h_selected_pion_incidentE] (double reco_firstEntryIncident, double reco_beam_interactingEnergy) {
            //bin that energy falls into is (int) energy/nbins + 1
            int binNumber_initEnergy = (int) reco_firstEntryIncident / bin_size_inc + 1;
            int binNumber_interEnergy = (int) reco_beam_interactingEnergy / bin_size_inc + 1;
            //how many bins need to be filled
            //starting at interacting energy (lower than first incident energy)
            //std::cout << "Incident Histo" << std::endl;
            //std::cout << "Interacting Energy = " << reco_beam_interactingEnergy << std::endl;
            //std::cout << "starting to fill at bin = " << binNumber_interEnergy << std::endl;
            //std::cout << "First incident Energy = " << reco_firstEntryIncident << std::endl;
            //std::cout << "ending to fill at bin = " << binNumber_initEnergy << std::endl;

            //if(binNumber_initEnergy < 0 || binNumber_interEnergy < 0) return;
            
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
            if( i > 0 && i <= nBin_inc){ //make sure we don't go outside of bin range
               h_selected_pion_incidentE->SetBinContent( i, h_selected_pion_incidentE->GetBinContent(i) + 1 ); 
                  //increment previous bin content by one. should handle number of entries in hist. unlike AddBinContent
            };
            };
      }
   ,{"reco_firstEntryIncident", "reco_interactingKE"});

   h_selected_pion_incidentE->Sumw2(0);
   h_selected_pion_incidentE->Rebin( bin_size_int/bin_size_inc );
   h_selected_pion_incidentE->Scale( 1 / (bin_size_int/bin_size_inc) );
   h_selected_pion_incidentE->Write();



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
      .Filter("primary_ends_inAPA3")
      .Filter("true_primPionInel");


   auto mcInteracting_selected_abs = mcInteracting_selected_allPrimaryPi
      .Filter("has_noPion_daughter && !(has_shower_nHits_distance)")
      .Filter("true_absSignal && true_pion_daughter == 0");

   mcInteracting_selected_abs
      .Foreach( [h_selected_abs_interactingE] (double reco_beam_interactingEnergy){
            h_selected_abs_interactingE->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   h_selected_abs_interactingE->Sumw2(0);
   h_selected_abs_interactingE->Write();

   auto mcInteracting_selected_cex = mcInteracting_selected_allPrimaryPi
      .Filter("has_noPion_daughter && has_shower_nHits_distance")
      .Filter("true_chexSignal && true_pion_daughter == 0");

   mcInteracting_selected_cex
      .Foreach( [h_selected_cex_interactingE] (double reco_beam_interactingEnergy){
            h_selected_cex_interactingE->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   h_selected_cex_interactingE->Sumw2(0);
   h_selected_cex_interactingE->Write();

   auto mcInteracting_selected_totInel = mcInteracting_selected_allPrimaryPi;

   mcInteracting_selected_totInel
      .Foreach( [h_selected_totInel_interactingE] (double reco_beam_interactingEnergy){
            h_selected_totInel_interactingE->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   h_selected_totInel_interactingE->Sumw2(0);
   h_selected_totInel_interactingE->Write();

   //=====================================================
   //            Prepare BetheBloch Mean for each Bin 
   //            QUESTION: take betheBloch of Pion or Muon?? Comparison to data fits better muon Bethe... 
   //            at hihger momentum ~400-800 they anway are almost the same
   //=====================================================
   TH1D* h_betheMean_pion = new TH1D("h_betheMean_pion", "Mean Energy Loss", nBin_int, eEnd, eStart);

   //fill histo with Mean dEdX of bin center
   for(int i = 1; i <= nBin_int; i++){
      h_betheMean_pion->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size_int  , mass_pion) );
   };

   h_betheMean_pion->Write();

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
   //    Absorption, True Interactions, True Energy
   //------------------------------------------------------
   //


   TH1D* h_xs_selected_abs = (TH1D*) h_selected_pion_incidentE->Clone("h_xs_selected_abs");
   TH1D* dummy_abs = (TH1D*) h_selected_pion_incidentE->Clone("h_xs_selected_abs");
   dummy_abs->Add( h_selected_abs_interactingE, -1);
   h_xs_selected_abs->Divide( dummy_abs );
   for(int i = 1; i <= nBin_int; i++) h_xs_selected_abs->SetBinContent(i, log( h_xs_selected_abs->GetBinContent(i) ));
   h_xs_selected_abs->Multiply( h_betheMean_pion );
   h_xs_selected_abs->Scale( scale_factor );

   /* Less Accurate when Nint !<< Ninc
      TH1D* h_xs_selected_abs = (TH1D*) h_selected_abs_interactingE->Clone("h_xs_selected_abs");
      h_xs_selected_abs->Divide( h_selected_pion_incidentE );
      h_xs_selected_abs->Multiply( h_betheMean_pion );
      h_xs_selected_abs->Scale( factor_mbarn*scale_factor );
      */

   //Try Bin Errors from error propagation derivation of sigma(E)
   //error is like error of binomial factor*sqrt( p(1-p)/Ninc) where p=Nint/Ninc

   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_abs_interactingE->GetBinContent(i) / h_selected_pion_incidentE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_pion->GetBinContent(i);

      h_xs_selected_abs->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_selected_abs->Write();


   //------------------------------------------------------
   //    ChargeExchange,  Selected Interactions, Reco Energy
   //------------------------------------------------------
   //

   TH1D* h_xs_selected_cex = (TH1D*) h_selected_pion_incidentE->Clone("h_xs_selected_cex");
   TH1D* dummy_cex = (TH1D*) h_selected_pion_incidentE->Clone("h_xs_selected_cex");
   dummy_cex->Add( h_selected_cex_interactingE, -1);
   h_xs_selected_cex->Divide( dummy_cex );
   for(int i = 1; i <= nBin_int; i++) h_xs_selected_cex->SetBinContent(i, log( h_xs_selected_cex->GetBinContent(i) ));
   h_xs_selected_cex->Multiply( h_betheMean_pion );
   h_xs_selected_cex->Scale( scale_factor );

   /* Less Accurate when Nint !<< Ninc
      TH1D* h_xs_selected_cex = (TH1D*) h_selected_cex_interactingE->Clone("h_xs_selected_cex");
      h_xs_selected_cex->Divide( h_selected_pion_incidentE );
      h_xs_selected_cex->Multiply( h_betheMean_pion );
      h_xs_selected_cex->Scale( factor_mbarn*scale_factor );
      */

   //Bin Errors
   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_cex_interactingE->GetBinContent(i) / h_selected_pion_incidentE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_pion->GetBinContent(i);

      h_xs_selected_cex->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_selected_cex->Write();
   //------------------------------------------------------
   //    Inelastic, Selected Interactions, Reco Energy
   //------------------------------------------------------
   //

   TH1D* h_xs_selected_totInel = (TH1D*) h_selected_pion_incidentE->Clone("h_xs_selected_totInel");
   TH1D* dummy_totInel = (TH1D*) h_selected_pion_incidentE->Clone("h_xs_selected_totInel");
   dummy_totInel->Add( h_selected_totInel_interactingE, -1);
   h_xs_selected_totInel->Divide( dummy_totInel );
   for(int i = 1; i <= nBin_int; i++) h_xs_selected_totInel->SetBinContent(i, log( h_xs_selected_totInel->GetBinContent(i) ));
   h_xs_selected_totInel->Multiply( h_betheMean_pion );
   h_xs_selected_totInel->Scale( scale_factor );

   /* Less Accurate when Nint !<< Ninc
      TH1D* h_xs_selected_totInel = (TH1D*) h_selected_totInel_interactingE->Clone("h_xs_selected_totInel");
      h_xs_selected_totInel->Divide( h_selected_pion_incidentE );
      h_xs_selected_totInel->Multiply( h_betheMean_pion );
      h_xs_selected_totInel->Scale( factor_mbarn*scale_factor );
      */

   //Bin Errors
   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_totInel_interactingE->GetBinContent(i) / h_selected_pion_incidentE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_pion->GetBinContent(i);

      h_xs_selected_totInel->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_selected_totInel->Write();


   //=====================================================*
   //            Plotting and Style
   //=====================================================
   //

   TCanvas *c_abs = new TCanvas("c_abs", "c_abs");
   gPad->SetGrid(1,1);
   h_xs_selected_abs->SetTitle( "Selected Absorption;Kinetic Energy (MeV); #sigma (mb)");
   h_xs_selected_abs->GetXaxis()->SetRangeUser(400,900);
   h_xs_selected_abs->GetXaxis()->SetNdivisions(1020);
   h_xs_selected_abs->GetYaxis()->SetNdivisions(1020);

   abs_KE->SetTitle( "Absorption;Kinetic Energy (MeV); #sigma (mb)");
   abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   abs_KE->SetLineColor(kRed);
   abs_KE->SetLineWidth(3);
   abs_KE->Draw("AC");
   h_xs_selected_abs->SetMarkerSize(0.7);
   h_xs_selected_abs->Draw("PE0 SAME");

   c_abs->Write();

   TCanvas *c_cex = new TCanvas("c_cex", "c_cex");
   gPad->SetGrid(1,1);
   h_xs_selected_cex->SetTitle( "Selected Charge Exchange;Kinetic Energy (MeV); #sigma (mb)");
   h_xs_selected_cex->GetXaxis()->SetRangeUser(400,900);
   h_xs_selected_cex->GetXaxis()->SetNdivisions(1020);
   h_xs_selected_cex->GetYaxis()->SetNdivisions(1020);

   cex_KE->SetTitle( "Charge Exchange;Kinetic Energy (MeV); #sigma (mb)");
   cex_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   cex_KE->GetYaxis()->SetRangeUser(0, 200);
   cex_KE->SetLineColor(kRed);
   cex_KE->SetLineWidth(3);
   cex_KE->Draw("AC"); 
   h_xs_selected_cex->SetMarkerSize(0.7);
   h_xs_selected_cex->Draw("PE0 SAME");

   c_cex->Write();

   TCanvas *c_totInel = new TCanvas("c_totInel", "c_totInel");
   gPad->SetGrid(1,1);
   h_xs_selected_totInel->SetTitle( "Selected Total Inelastic;Kinetic Energy (MeV); #sigma (mb)");
   h_xs_selected_totInel->GetXaxis()->SetRangeUser(400,900);
   h_xs_selected_totInel->GetXaxis()->SetNdivisions(1020);
   h_xs_selected_totInel->GetYaxis()->SetNdivisions(1020);

   totInel_KE->SetTitle( "Total Inelastic;Kinetic Energy (MeV); #sigma (mb)");
   totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   totInel_KE->SetLineColor(kRed);
   totInel_KE->SetLineWidth(3);
   totInel_KE->Draw("AC");
   h_xs_selected_totInel->SetMarkerSize(0.7);
   h_xs_selected_totInel->Draw("PE0 SAME");

   c_totInel->Write();

   TCanvas *c_all = new TCanvas("c_all", "c_all");
   gPad->SetGrid(1,1);

   cex_KE->SetLineColor(kGreen);
   abs_KE->SetLineColor(kBlue);
   mg->GetXaxis()->SetRangeUser(0,1100);
   mg->GetXaxis()->SetNdivisions(1020);
   mg->Draw("AC");
   //h_xs_selected_totInel->Draw("PE0 SAME");
   h_xs_selected_cex->Draw("PE0 SAME");
   h_xs_selected_abs->Draw("PE0 SAME");

   c_all->Write();




   output->Write();
   //f1.Close();
   return 0;
}

