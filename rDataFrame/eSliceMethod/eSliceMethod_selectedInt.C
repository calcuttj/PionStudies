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

int eSliceMethod_selectedInt(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(inputTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   
   //output file
   TFile *output = new TFile ("output_eSliceMethod_selectedEvents.root", "RECREATE");

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
   mg->Add(totInel_KE);
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
   hist_bethe_mean( KE_in_pion, mass_pion, fit_pitch_mean, bethe_pi_mean);
   //hist_bethe_mpv( 885.7, mass_muon, fit_pitch_mean, bethe_mu_mpv); //885.7 is mean of dsitribution of reco_beam_incidentEnergies[0]
   //hist_bethe_mean( 885.7, mass_muon, fit_pitch_mean, bethe_mu_mean);
   hist_bethe_mpv( KE_in_muon, mass_muon, fit_pitch_mean, bethe_mu_mpv);
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
   for(int i=1; i <= dE_product_fit_dEdX_pitch.GetNbinsX(); i++){
      temp += dE_product_fit_dEdX_pitch.GetBinContent(i);
      runningSum_dE->SetBinContent( i, temp);
   };

   dE_product_fit_dEdX_pitch.Write();
   runningSum_dE->Write();
//--------------------------------------------------------
      
   //Incident Histogram and Interacting Histogram, interacting for different Processes
   //Energies are kinetic Energy of beam particle

   TH1D* h_selected_pion_incidentRecoE = new TH1D("h_selected_pion_incidentRecoE", "Incident Selected Pion", nBin_inc, eEnd, eStart);
   h_selected_pion_incidentRecoE->GetXaxis()->SetTitle("Reco KE (MeV)");

   TH1D* h_selected_abs_interactingRecoE = new TH1D("h_selected_abs_interactingRecoE", "Interacting Selected ABS Reco Energy", nBin_int, eEnd, eStart);
   h_selected_abs_interactingRecoE->GetXaxis()->SetTitle("Reco KE (MeV)");
   TH1D* h_selected_cex_interactingRecoE = new TH1D("h_selected_cex_interactingRecoE", "Interacting Selected CEX Reco Energy", nBin_int, eEnd, eStart);
   h_selected_cex_interactingRecoE->GetXaxis()->SetTitle("Reco KE (MeV)");
   TH1D* h_selected_totInel_interactingRecoE = new TH1D("h_selected_totInel_interactingRecoE", "Interacting Selected total INEL Reco Energy", nBin_int, eEnd, eStart);
   h_selected_totInel_interactingRecoE->GetXaxis()->SetTitle("Reco KE (MeV)");
   
   TH1D* h_selected_pion_incidentTrueE = new TH1D("h_selected_pion_incidentTrueE", "Incident Selected Pion", nBin_int, eEnd, eStart);
   h_selected_pion_incidentTrueE->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_selected_abs_interactingTrueE = new TH1D("h_selected_abs_interactingTrueE", "Interacting Selected ABS True Energy", nBin_int, eEnd, eStart);
   h_selected_abs_interactingTrueE->GetXaxis()->SetTitle("True KE (MeV)");
   TH1D* h_selected_cex_interactingTrueE = new TH1D("h_selected_cex_interactingTrueE", "Interacting Selected CEX True Energy", nBin_int, eEnd, eStart);
   h_selected_cex_interactingTrueE->GetXaxis()->SetTitle("True KE (MeV)");
   TH1D* h_selected_totInel_interactingTrueE = new TH1D("h_selected_totInel_interactingTrueE", "Interacting Selected total True INEL Energy", nBin_int, eEnd, eStart);
   h_selected_totInel_interactingTrueE->GetXaxis()->SetTitle("True KE (MeV)");


   //Initial Filters for all events
   auto frame_filter = frame
      .Filter("true_beam_endZ > 0")
      //.Filter("true_beam_PDG == abs(211)")
      //.Filter("primary_isBeamType && passBeamCut && passBeamQuality");
      .Filter("primary_isBeamType && passBeamCutBI");    //particles that make it into TPC, couldn't be considered otherwise


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
      //.Range(3)
      .Define("true_firstEntryIncident", firstIncident, {"true_beam_incidentEnergies"})
      .Define("true_KEint_fromEndP", [](double true_beam_endP){
            true_beam_endP = 1000*true_beam_endP; //convert GeV -> MeV
            double endKE = sqrt( pow(true_beam_endP,2) + pow(mass_pion,2)  ) - mass_pion;
            return endKE;}
            ,{"true_beam_endP"})
      
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
            },{"reco_beam_calo_wire"});

      //.Filter("reco_interactingKE > 0 ");
 
    auto h2_wire_KEinc = mcIncident_selected_primaryPi
       .Histo2D({"recoKEinc_vs_wire", "reco first incident KEnergy at reco wire", 400,0,800, 200,0,1000}, "reco_interacting_wire", "reco_firstEntryIncident");

   h2_wire_KEinc->Write();

   auto h2_wire_KEint = mcIncident_selected_primaryPi
       .Histo2D({"recoKEint_vs_wire", "reco interacting KEnergy at reco wire", 400,0,800, 200,0,1000}, "reco_interacting_wire", "reco_interactingKE");

   h2_wire_KEint->Write();

   auto h2_reco_true_KE_primPi = mcIncident_selected_primaryPi
       .Histo2D({"h2_reco_true_KE_primPi", "true KE vs reco KE incident Prim Pi", 400,0,1200, 400,0,1200}, "true_KEint_fromEndP", "reco_interactingKE");
   h2_reco_true_KE_primPi->Write();


   //========================================================
   //Build the Incident Histogram
   //---------
   //take the interacting energy, that slice is filled too

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

   h_selected_pion_incidentTrueE->Sumw2(0);
   h_selected_pion_incidentTrueE->Rebin( bin_size_int/bin_size_int );
   h_selected_pion_incidentTrueE->Scale( 1 / (bin_size_int/bin_size_int) );
   h_selected_pion_incidentTrueE->Write();
   
   h_selected_pion_incidentRecoE->Sumw2(0);
   h_selected_pion_incidentRecoE->Rebin( bin_size_int/bin_size_inc );
   h_selected_pion_incidentRecoE->Scale( 1 / (bin_size_int/bin_size_inc) );
   h_selected_pion_incidentRecoE->Write();



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
      .Filter("has_noPion_daughter && !(has_shower_nHits_distance)");

   mcInteracting_selected_abs
     .Foreach( [h_selected_abs_interactingTrueE] (double true_beam_interactingEnergy){
            h_selected_abs_interactingTrueE->Fill(true_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"}); 

    mcInteracting_selected_abs 
       .Foreach( [h_selected_abs_interactingRecoE] (double reco_beam_interactingEnergy){
            h_selected_abs_interactingRecoE->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   h_selected_abs_interactingTrueE->Sumw2(0);
   h_selected_abs_interactingTrueE->Write();
   
   h_selected_abs_interactingRecoE->Sumw2(0);
   h_selected_abs_interactingRecoE->Write();

   //DEBUG

   auto h2_reco_true_KE_abs= mcInteracting_selected_abs
       .Histo2D({"h2_reco_true_KE_abs", "true KE vs reco KE selected ABS", 400,0,1.2, 400,0,1200}, "true_beam_endP", "reco_interactingKE");
   h2_reco_true_KE_abs->Write();

   auto mcInteracting_selected_cex = mcInteracting_selected_allPrimaryPi
      .Filter("has_noPion_daughter && has_shower_nHits_distance");

   mcInteracting_selected_cex
     .Foreach( [h_selected_cex_interactingTrueE] (double true_beam_interactingEnergy){
            h_selected_cex_interactingTrueE->Fill(true_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});
   
   mcInteracting_selected_cex
     .Foreach( [h_selected_cex_interactingRecoE] (double reco_beam_interactingEnergy){
            h_selected_cex_interactingRecoE->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   h_selected_cex_interactingTrueE->Sumw2(0);
   h_selected_cex_interactingTrueE->Write();
   
   h_selected_cex_interactingRecoE->Sumw2(0);
   h_selected_cex_interactingRecoE->Write();

   auto mcInteracting_selected_totInel = mcInteracting_selected_allPrimaryPi;

   mcInteracting_selected_totInel
     .Foreach( [h_selected_totInel_interactingTrueE] (double true_beam_interactingEnergy){
            h_selected_totInel_interactingTrueE->Fill(true_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});
   
   mcInteracting_selected_totInel
     .Foreach( [h_selected_totInel_interactingRecoE] (double reco_beam_interactingEnergy){
            h_selected_totInel_interactingRecoE->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   h_selected_totInel_interactingTrueE->Sumw2(0);
   h_selected_totInel_interactingTrueE->Write();
   
   h_selected_totInel_interactingRecoE->Sumw2(0);
   h_selected_totInel_interactingRecoE->Write();

   //=====================================================
   //            Prepare BetheBloch Mean for each Bin 
   //            QUESTION: take betheBloch of Pion or Muon?? Comparison to data fits better muon Bethe... 
   //            at hihger momentum ~400-800 they anway are almost the same
   //=====================================================
   TH1D* h_betheMean_pion = new TH1D("h_betheMean_pion", "Mean Energy Loss", nBin_int, eEnd, eStart);

   //fill histo with Mean dEdX of bin center
   for(int i = 1; i <= nBin_int; i++){
      h_betheMean_pion->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size_int  , mass_muon) );
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
   //    Absorption, Selected Interactions Reconstrucetd Energy
   //------------------------------------------------------
   //
   

   TH1D* h_xs_RecoE_selected_abs = (TH1D*) h_selected_pion_incidentRecoE->Clone("h_xs_RecoE_selected_abs");
   TH1D* dummy_recoE_abs = (TH1D*) h_selected_pion_incidentRecoE->Clone("h_xs_RecoE_selected_abs");
   dummy_recoE_abs->Add( h_selected_abs_interactingRecoE, -1);
   h_xs_RecoE_selected_abs->Divide( dummy_recoE_abs );
   for(int i = 1; i <= nBin_int; i++) h_xs_RecoE_selected_abs->SetBinContent(i, log( h_xs_RecoE_selected_abs->GetBinContent(i) ));
   h_xs_RecoE_selected_abs->Multiply( h_betheMean_pion );
   h_xs_RecoE_selected_abs->Scale( scale_factor );

   TH1D* h_xs_TrueE_selected_abs = (TH1D*) h_selected_pion_incidentTrueE->Clone("h_xs_TrueE_selected_abs");
   TH1D* dummy_TrueE_abs = (TH1D*) h_selected_pion_incidentTrueE->Clone("h_xs_TrueE_selected_abs");
   dummy_TrueE_abs->Add( h_selected_abs_interactingTrueE, -1);
   h_xs_TrueE_selected_abs->Divide( dummy_TrueE_abs );
   for(int i = 1; i <= nBin_int; i++) h_xs_TrueE_selected_abs->SetBinContent(i, log( h_xs_TrueE_selected_abs->GetBinContent(i) ));
   h_xs_TrueE_selected_abs->Multiply( h_betheMean_pion );
   h_xs_TrueE_selected_abs->Scale( scale_factor );
   
   /* Less Accurate when Nint !<< Ninc
   TH1D* h_xs_RecoE_selected_abs = (TH1D*) h_selected_abs_interactingRecoE->Clone("h_xs_RecoE_selected_abs");
   h_xs_RecoE_selected_abs->Divide( h_selected_pion_incidentRecoE );
   h_xs_RecoE_selected_abs->Multiply( h_betheMean_pion );
   h_xs_RecoE_selected_abs->Scale( factor_mbarn*scale_factor );
   */

   //Try Bin Errors from error propagation derivation of sigma(E)
   //error is like error of binomial factor*sqrt( p(1-p)/Ninc) where p=Nint/Ninc
   
   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_abs_interactingRecoE->GetBinContent(i) / h_selected_pion_incidentRecoE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentRecoE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_pion->GetBinContent(i);
      
      h_xs_RecoE_selected_abs->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_RecoE_selected_abs->Write();
   
   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_abs_interactingTrueE->GetBinContent(i) / h_selected_pion_incidentTrueE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentTrueE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_pion->GetBinContent(i);
      
      h_xs_TrueE_selected_abs->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_TrueE_selected_abs->Write();


   //------------------------------------------------------
   //    ChargeExchange,  Selected Interactions, Reco Energy
   //------------------------------------------------------
   //
   
   TH1D* h_xs_RecoE_selected_cex = (TH1D*) h_selected_pion_incidentRecoE->Clone("h_xs_RecoE_selected_cex");
   TH1D* dummy_recoE_cex = (TH1D*) h_selected_pion_incidentRecoE->Clone("h_xs_RecoE_selected_cex");
   dummy_recoE_cex->Add( h_selected_cex_interactingRecoE, -1);
   h_xs_RecoE_selected_cex->Divide( dummy_recoE_cex );
   for(int i = 1; i <= nBin_int; i++) h_xs_RecoE_selected_cex->SetBinContent(i, log( h_xs_RecoE_selected_cex->GetBinContent(i) ));
   h_xs_RecoE_selected_cex->Multiply( h_betheMean_pion );
   h_xs_RecoE_selected_cex->Scale( scale_factor );
   
   TH1D* h_xs_TrueE_selected_cex = (TH1D*) h_selected_pion_incidentTrueE->Clone("h_xs_TrueE_selected_cex");
   TH1D* dummy_TrueE_cex = (TH1D*) h_selected_pion_incidentTrueE->Clone("h_xs_TrueE_selected_cex");
   dummy_TrueE_cex->Add( h_selected_cex_interactingTrueE, -1);
   h_xs_TrueE_selected_cex->Divide( dummy_TrueE_cex );
   for(int i = 1; i <= nBin_int; i++) h_xs_TrueE_selected_cex->SetBinContent(i, log( h_xs_TrueE_selected_cex->GetBinContent(i) ));
   h_xs_TrueE_selected_cex->Multiply( h_betheMean_pion );
   h_xs_TrueE_selected_cex->Scale( scale_factor );
   
   /* Less Accurate when Nint !<< Ninc
   TH1D* h_xs_RecoE_selected_cex = (TH1D*) h_selected_cex_interactingRecoE->Clone("h_xs_RecoE_selected_cex");
   h_xs_RecoE_selected_cex->Divide( h_selected_pion_incidentRecoE );
   h_xs_RecoE_selected_cex->Multiply( h_betheMean_pion );
   h_xs_RecoE_selected_cex->Scale( factor_mbarn*scale_factor );
   */

  //Bin Errors
   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_cex_interactingRecoE->GetBinContent(i) / h_selected_pion_incidentRecoE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentRecoE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_pion->GetBinContent(i);
      
      h_xs_RecoE_selected_cex->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_RecoE_selected_cex->Write();

   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_cex_interactingTrueE->GetBinContent(i) / h_selected_pion_incidentTrueE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentTrueE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_pion->GetBinContent(i);
      
      h_xs_TrueE_selected_cex->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_TrueE_selected_cex->Write();
   //------------------------------------------------------
   //    Inelastic, Selected Interactions, Reco Energy
   //------------------------------------------------------
   //
   
   TH1D* h_xs_RecoE_selected_totInel = (TH1D*) h_selected_pion_incidentRecoE->Clone("h_xs_RecoE_selected_totInel");
   TH1D* dummy_recoE_totInel = (TH1D*) h_selected_pion_incidentRecoE->Clone("h_xs_RecoE_selected_totInel");
   dummy_recoE_totInel->Add( h_selected_totInel_interactingRecoE, -1);
   h_xs_RecoE_selected_totInel->Divide( dummy_recoE_totInel );
   for(int i = 1; i <= nBin_int; i++) h_xs_RecoE_selected_totInel->SetBinContent(i, log( h_xs_RecoE_selected_totInel->GetBinContent(i) ));
   h_xs_RecoE_selected_totInel->Multiply( h_betheMean_pion );
   h_xs_RecoE_selected_totInel->Scale( scale_factor );

   TH1D* h_xs_TrueE_selected_totInel = (TH1D*) h_selected_pion_incidentTrueE->Clone("h_xs_TrueE_selected_totInel");
   TH1D* dummy_TrueE_totInel = (TH1D*) h_selected_pion_incidentTrueE->Clone("h_xs_TrueE_selected_totInel");
   dummy_TrueE_totInel->Add( h_selected_totInel_interactingTrueE, -1);
   h_xs_TrueE_selected_totInel->Divide( dummy_TrueE_totInel );
   for(int i = 1; i <= nBin_int; i++) h_xs_TrueE_selected_totInel->SetBinContent(i, log( h_xs_TrueE_selected_totInel->GetBinContent(i) ));
   h_xs_TrueE_selected_totInel->Multiply( h_betheMean_pion );
   h_xs_TrueE_selected_totInel->Scale( scale_factor );
   
   /* Less Accurate when Nint !<< Ninc
   TH1D* h_xs_RecoE_selected_totInel = (TH1D*) h_selected_totInel_interactingRecoE->Clone("h_xs_RecoE_selected_totInel");
   h_xs_RecoE_selected_totInel->Divide( h_selected_pion_incidentRecoE );
   h_xs_RecoE_selected_totInel->Multiply( h_betheMean_pion );
   h_xs_RecoE_selected_totInel->Scale( factor_mbarn*scale_factor );
   */

   //Bin Errors
   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_totInel_interactingRecoE->GetBinContent(i) / h_selected_pion_incidentRecoE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentRecoE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_pion->GetBinContent(i);
      
      h_xs_RecoE_selected_totInel->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_RecoE_selected_totInel->Write();

   for(int i=1; i <= nBin_int; i++){

      double p = h_selected_totInel_interactingTrueE->GetBinContent(i) / h_selected_pion_incidentTrueE->GetBinContent(i);
      double nInc_i = h_selected_pion_incidentTrueE->GetBinContent(i);
      double help_factor = scale_factor*h_betheMean_pion->GetBinContent(i);
      
      h_xs_TrueE_selected_totInel->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
   };
   h_xs_TrueE_selected_totInel->Write();
 

   //=====================================================*
   //            Plotting and Style
   //=====================================================
   //
   //h_xs_RecoE_selected_abs->Scale(2);
   //h_xs_RecoE_selected_cex->Scale(2);
   //h_xs_RecoE_selected_totInel->Scale(2);
   
   TCanvas *c_RecoE_abs = new TCanvas("c_RecoE_abs", "c_RecoE_abs");
   gPad->SetGrid(1,1);
   h_xs_RecoE_selected_abs->SetTitle( "Selected Absorption;Reco Kinetic Energy (MeV); #sigma (mb)");
   h_xs_RecoE_selected_abs->GetXaxis()->SetRangeUser(400,900);
   h_xs_RecoE_selected_abs->GetXaxis()->SetNdivisions(1020);
   h_xs_RecoE_selected_abs->GetYaxis()->SetNdivisions(1020);
   
   abs_KE->SetTitle( "Absorption;Reco Kinetic Energy (MeV); #sigma (mb)");
   abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   abs_KE->SetLineColor(kRed);
   abs_KE->SetLineWidth(3);
   abs_KE->Draw("AC");
   h_xs_RecoE_selected_abs->SetMarkerSize(0.7);
   h_xs_RecoE_selected_abs->Draw("PE0 SAME");

   c_RecoE_abs->Write();
   
   TCanvas *c_RecoE_cex = new TCanvas("c_RecoE_cex", "c_RecoE_cex");
   gPad->SetGrid(1,1);
   h_xs_RecoE_selected_cex->SetTitle( "Selected Charge Exchange;Reco Kinetic Energy (MeV); #sigma (mb)");
   h_xs_RecoE_selected_cex->GetXaxis()->SetRangeUser(400,900);
   h_xs_RecoE_selected_cex->GetXaxis()->SetNdivisions(1020);
   h_xs_RecoE_selected_cex->GetYaxis()->SetNdivisions(1020);
   
   cex_KE->SetTitle( "Charge Exchange;Reco Kinetic Energy (MeV); #sigma (mb)");
   cex_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   cex_KE->GetYaxis()->SetRangeUser(0, 200);
   cex_KE->SetLineColor(kRed);
   cex_KE->SetLineWidth(3);
   cex_KE->Draw("AC"); 
   h_xs_RecoE_selected_cex->SetMarkerSize(0.7);
   h_xs_RecoE_selected_cex->Draw("PE0 SAME");

   c_RecoE_cex->Write();

   TCanvas *c_RecoE_totInel = new TCanvas("c_RecoE_totInel", "c_RecoE_totInel");
   gPad->SetGrid(1,1);
   h_xs_RecoE_selected_totInel->SetTitle( "Selected Total Inelastic;Reco Kinetic Energy (MeV); #sigma (mb)");
   h_xs_RecoE_selected_totInel->GetXaxis()->SetRangeUser(400,900);
   h_xs_RecoE_selected_totInel->GetXaxis()->SetNdivisions(1020);
   h_xs_RecoE_selected_totInel->GetYaxis()->SetNdivisions(1020);
   
   totInel_KE->SetTitle( "Total Inelastic;Reco Kinetic Energy (MeV); #sigma (mb)");
   totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   totInel_KE->SetLineColor(kRed);
   totInel_KE->SetLineWidth(3);
   totInel_KE->Draw("AC");
   h_xs_RecoE_selected_totInel->SetMarkerSize(0.7);
   h_xs_RecoE_selected_totInel->Draw("PE0 SAME");

   c_RecoE_totInel->Write();

   TCanvas *c_RecoE_all = new TCanvas("c_RecoE_all", "c_RecoE_all");
   gPad->SetGrid(1,1);
   
   cex_KE->SetLineColor(kGreen);
   abs_KE->SetLineColor(kBlue);
   mg->GetXaxis()->SetRangeUser(0,1100);
   mg->GetXaxis()->SetNdivisions(1020);
   mg->Draw("AC");
   h_xs_RecoE_selected_totInel->Draw("PE0 SAME");
   h_xs_RecoE_selected_cex->Draw("PE0 SAME");
   h_xs_RecoE_selected_abs->Draw("PE0 SAME");

   c_RecoE_all->Write();

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
   
   TCanvas *c_TrueE_cex = new TCanvas("c_TrueE_cex", "c_TrueE_cex");
   gPad->SetGrid(1,1);
   h_xs_TrueE_selected_cex->SetTitle( "Selected Charge Exchange;True Kinetic Energy (MeV); #sigma (mb)");
   h_xs_TrueE_selected_cex->GetXaxis()->SetRangeUser(400,900);
   h_xs_TrueE_selected_cex->GetXaxis()->SetNdivisions(1020);
   h_xs_TrueE_selected_cex->GetYaxis()->SetNdivisions(1020);
   
   cex_KE->SetTitle( "Charge Exchange;True Kinetic Energy (MeV); #sigma (mb)");
   cex_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   cex_KE->GetYaxis()->SetRangeUser(0, 200);
   cex_KE->SetLineColor(kRed);
   cex_KE->SetLineWidth(3);
   cex_KE->Draw("AC"); 
   h_xs_TrueE_selected_cex->SetMarkerSize(0.7);
   h_xs_TrueE_selected_cex->Draw("PE0 SAME");

   c_TrueE_cex->Write();

   TCanvas *c_TrueE_totInel = new TCanvas("c_TrueE_totInel", "c_TrueE_totInel");
   gPad->SetGrid(1,1);
   h_xs_TrueE_selected_totInel->SetTitle( "Selected Total Inelastic;True Kinetic Energy (MeV); #sigma (mb)");
   h_xs_TrueE_selected_totInel->GetXaxis()->SetRangeUser(400,900);
   h_xs_TrueE_selected_totInel->GetXaxis()->SetNdivisions(1020);
   h_xs_TrueE_selected_totInel->GetYaxis()->SetNdivisions(1020);
   
   totInel_KE->SetTitle( "Total Inelastic;True Kinetic Energy (MeV); #sigma (mb)");
   totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   totInel_KE->SetLineColor(kRed);
   totInel_KE->SetLineWidth(3);
   totInel_KE->Draw("AC");
   h_xs_TrueE_selected_totInel->SetMarkerSize(0.7);
   h_xs_TrueE_selected_totInel->Draw("PE0 SAME");

   c_TrueE_totInel->Write();

   TCanvas *c_TrueE_all = new TCanvas("c_TrueE_all", "c_TrueE_all");
   gPad->SetGrid(1,1);
   
   cex_KE->SetLineColor(kGreen);
   abs_KE->SetLineColor(kBlue);
   mg->GetXaxis()->SetRangeUser(0,1100);
   mg->GetXaxis()->SetNdivisions(1020);
   mg->Draw("AC");
   h_xs_TrueE_selected_totInel->Draw("PE0 SAME");
   h_xs_TrueE_selected_cex->Draw("PE0 SAME");
   h_xs_TrueE_selected_abs->Draw("PE0 SAME");

   c_TrueE_all->Write();
  


   //output->Write();
   //f1.Close();
   return 0;
}

