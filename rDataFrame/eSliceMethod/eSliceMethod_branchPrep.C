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

//This Macro prepares the interacting and incident energy-branches for the eSliceMethod files.
//Provide output name

//***********************
//Main Function

int eSliceMethod_branchPrep(const string mcFilepath, const string outputName){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);

   TFile f2("fit_mc_Prod4_dEdX_pitch_05_19_21.root", "UPDATE");
   TH1D *fit_dEdX_lifetime_mpv = (TH1D*)f2.Get("dEdX_mpv_lifetime"); //mean value corrected for lifetime
   TH1D *fit_pitch_mean = (TH1D*)f2.Get("fit_mc_pitch_mean");

   TH1D *fit_dEdX_lifetime_mpv_SCEcorr = (TH1D*)f2.Get("fit_mc_dEdX_SCEcorr_mpv"); //mean value corrected for lifetime
   TH1D *fit_pitch_mean_SCEcorr = (TH1D*)f2.Get("fit_mc_pitch_SCEcorr_mean");
   TFile *output = new TFile ("eSliceMethod_energyDeposit_05_19_21.root", "RECREATE");

   output->cd();
   fit_dEdX_lifetime_mpv->Write();
   fit_pitch_mean->Write();


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

   //unclibrated

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

   //SCEcorr
   TH1D* dEdX_mean_calc_fit_bethe_SCEcorr = (TH1D*)bethe_frac_mu->Clone("dEdX_mean_calc_fit_bethe_SCEcorr");
   dEdX_mean_calc_fit_bethe_SCEcorr->Multiply(fit_dEdX_lifetime_mpv_SCEcorr);
   dEdX_mean_calc_fit_bethe_SCEcorr->Write();

   TH1D dE_product_fit_dEdX_pitch_SCEcorr = (*dEdX_mean_calc_fit_bethe_SCEcorr) * (*fit_pitch_mean_SCEcorr);
   dE_product_fit_dEdX_pitch_SCEcorr.SetName("dE_product_fit_dEdX_pitch_SCEcorr");

   TH1D *runningSum_dE_SCEcorr = new TH1D("runningSum_dE_SCEcorr", "", fit_pitch_mean_SCEcorr->GetNbinsX(), 1, fit_pitch_mean_SCEcorr->GetNbinsX() );

   temp = 0;
   for(int i=1; i <= dE_product_fit_dEdX_pitch_SCEcorr.GetNbinsX(); i++){
      temp += dE_product_fit_dEdX_pitch_SCEcorr.GetBinContent(i);
      runningSum_dE_SCEcorr->SetBinContent( i, temp);
   };

   dE_product_fit_dEdX_pitch_SCEcorr.Write();
   runningSum_dE_SCEcorr->Write();
   //--------------------------------------------------------

   //Initial Filters for all events
   auto frame_filter = frame
      .Filter("true_beam_endZ > 0");
 

   //Filter for different interaction types
   //
   //=====================================================
   //         BUILD Branches for true/reco Interacting / Incident energies
   //=====================================================
   //------------------------------------------------------
   //
   //for Interacting energy, the wire the primary particle interacted at should be the last one in reco_beam_calo_wire
   //

   auto mcIncident_selected_primaryPi = frame_filter      
      //.Range(100)
      .Define("true_firstEntryIncident", [](std::vector<double> &trajKE){
            return trajKE[0] - eLoss_mc_trueE;
            }
      ,{"true_beam_traj_KE"})
      
      .Define("true_interactingKE", [](double true_beam_endP){
            true_beam_endP = 1000*true_beam_endP; //convert GeV -> MeV
            double endKE = sqrt( pow(true_beam_endP,2) + pow(mass_pion,2)  ) - mass_pion;
            return endKE;}
            ,{"true_beam_endP"})
      //uncalibrated
      .Define("reco_firstEntryIncident", firstIncident, {"reco_beam_incidentEnergies"})
      //For now in Reco not subtratcting any energy Loss... can do later on if needed

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

      //SCEcorr only need to change interactingKE
      .Define("reco_interactingKE_SCEcorr", [runningSum_dE_SCEcorr](const std::vector<double> &reco_beam_calo_wire, double incidentE){
            double interactingWire = reco_beam_calo_wire[ reco_beam_calo_wire.size() ];
            double interactingKE;
            if(interactingWire >= 1 && interactingWire < runningSum_dE_SCEcorr->GetNbinsX()){
            interactingKE = incidentE - runningSum_dE_SCEcorr->GetBinContent(interactingWire);
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

   delete fit_dEdX_lifetime_mpv;
   delete fit_pitch_mean;
   mcIncident_selected_primaryPi.Snapshot("pionana/beamana", outputName);
   
   mcIncident_selected_primaryPi.Range(0,31700).Snapshot("pionana/beamana", "eSliceMethod_Prod4_mc_1GeV_part1_05_19_21.root");
   mcIncident_selected_primaryPi.Range(31701,0).Snapshot("pionana/beamana", "eSliceMethod_Prod4_mc_1GeV_part2_05_19_21.root");


   

   f2.Close();

   return 0;
}

