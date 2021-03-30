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

//Computing efficiencies after event Selection for a bin-by-bin for the incident and interacting histograms
//These efficiencies are produced from MC
//
//Eff = N true selected / N true available 
//Pur = N true selected / N total selected

//***********************
//Main Function

int effPur_eventSelection(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);

   //--------------------------------------------------------
   //
   //Read In relevant Histos from Fit for reco interacting and incident Energy estimation
   //
   //--------------------------------------------------------

   TFile f2("fit_mc_Prod4_dEdX_pitch_03_14_21.root");
   TH1D *fit_dEdX_lifetime_mpv = (TH1D*)f2.Get("dEdX_mpv_lifetime"); //mean value corrected for lifetime
   TH1D *fit_pitch_mean = (TH1D*)f2.Get("fit_mc_pitch_mean");

   //output file
   TFile *output = new TFile ("eventSelection_eff_pur.root", "RECREATE");
   
   //--------------------------------------------------------
   //Bethe Bloch MPV and Mean
   //--------------------------------------------------------

   TH1D* bethe_mu_mpv = new TH1D("betheMPV_muon", "", fit_pitch_mean->GetNbinsX(), 1, fit_pitch_mean->GetNbinsX() );
   bethe_mu_mpv->GetXaxis()->SetTitle("wire");  bethe_mu_mpv->GetYaxis()->SetTitle("dEdX (MeV/cm)");
   
   TH1D* bethe_mu_mean = new TH1D("betheMean_muon", "", fit_pitch_mean->GetNbinsX(), 1, fit_pitch_mean->GetNbinsX() );
   bethe_mu_mean->GetXaxis()->SetTitle("wire"); bethe_mu_mean->GetYaxis()->SetTitle("dEdX (MeV/cm)");
   
   hist_bethe_mpv( KE_in_muon, mass_muon, fit_pitch_mean, bethe_mu_mpv);
   hist_bethe_mean( KE_in_muon, mass_muon, fit_pitch_mean, bethe_mu_mean);

   bethe_mu_mpv->Write();
   bethe_mu_mean->Write();
   
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
   //
   //Create Relevant Branches that are used in eSlice Method
   //
   //Note: should somehow create prep-macro that does these branches and then they are stored in rDataFrame
   //
   //--------------------------------------------------------
   //
   auto frame_branches = frame
      //Filter only what is available for the selection, thus particles that really make it into the TPC (is that the only and sufficient Criterion?)
      .Filter("true_beam_endZ > 0")
      //MC true values
      .Define("true_firstEntryIncident", firstIncident, {"true_beam_incidentEnergies"})
      .Define("true_KEint_fromEndP", [](double true_beam_endP){
            true_beam_endP = 1000*true_beam_endP; //convert GeV -> MeV
            double endKE = sqrt( pow(true_beam_endP,2) + pow(mass_pion,2)  ) - mass_pion;
            return endKE;}
            ,{"true_beam_endP"})
      .Define("true_interacting_wire", "true_beam_endZ / 0.48")

      .Define("reco_firstEntryIncident", firstIncident, {"reco_beam_incidentEnergies"})
      //RECO values  
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

   //--------------------------------------------------------
   //
   //For the efficiency Calculation filter frame_branches for true absorption events 
   //
   //       pandora eff = pandora_trueSelected / N true available
   //       beamcut eff = beamcut_trueSelected / N true available from pandora
   //       selection eff = selection_trueSelected / N true available in incident primary pi sample
   //
   //       the inefficiencies that are for primary incident pions and abs will cancel anyway for the XS
   //--------------------------------------------------------
   //
   auto all_avail_abs = frame_branches
      .Filter("true_absSignal");

   auto pandora_all = frame_branches
      .Filter("primary_isBeamType");
   
   auto pandora_avail_abs = all_avail_abs 
      .Filter("primary_isBeamType");

   auto incidentPion_all = frame_branches
      .Filter("selected_incidentPion");

   auto incidentPion_avail_abs = all_avail_abs
      .Filter("selected_incidentPion");

   auto frame_branches
      .Filter("selected_abs");

   auto abs_select_avail_abs = incidentPion_avail_abs
      .Filter("selected_abs");

   //--------------------------------------------------------
   //
   //Start by building the interacting histos
   //Reco energy
   //Available = true_beam_endZ > 0 they interact IN TPC
   //--------------------------------------------------------
   //
   TH1D* h_trueAvailable = new TH1D("abs_trueAvailable", "Available (trueEndZ > 0) true abs events Interacting Reco Energy", nBin_int, eEnd, eStart);
   
   TH1D* h_sel_pandora = new TH1D("abs_sel_pandora", "Pandora reconstructed, trueAbs", nBin_int, eEnd, eStart);
   TH1D* h_sel_incidentPion_avail_abs = new TH1D("abs_sel_incidentPion_avail_abs", "Selected incident Pion, trueAbs", nBin_int, eEnd, eStart);
   TH1D* h_sel_abs_evSel = new TH1D("abs_sel_evSel", "Pandora reconstructed, trueAbs", nBin_int, eEnd, eStart);
   
   TH1D* h_eff_pandora = new TH1D("eff_pandora", "Efficiency Pandora trueAbs", nBin_int, eEnd, eStart);
   TH1D* h_eff_incidentPion_avail_abs = new TH1D("eff_incidentPion_avail_abs", "Efficiency selected incident Pion trueAbs", nBin_int, eEnd, eStart);
   TH1D* h_eff_abs_evSel = new TH1D("eff_abs_evSel", "Efficiency evSel trueAbs", nBin_int, eEnd, eStart);

   //Fill Interacting Histograms
   all_avail_abs
      .Foreach( [h_trueAvailable] (double reco_beam_interactingEnergy){
            h_trueAvailable->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   pandora_avail_abs
      .Foreach( [h_sel_pandora] (double reco_beam_interactingEnergy){
            h_sel_pandora->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});
   
   incidentPion_avail_abs
      .Foreach( [h_sel_incidentPion_avail_abs] (double reco_beam_interactingEnergy){
            h_sel_incidentPion_avail_abs->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});
   
   abs_select_avail_abs
      .Foreach( [h_sel_abs_evSel] (double reco_beam_interactingEnergy){
            h_sel_abs_evSel->Fill(reco_beam_interactingEnergy);}
            ,{"reco_interactingKE"});

   //Make Efficiency Histos

   h_eff_pandora->Divide(h_sel_pandora, h_trueAvailable);
   h_eff_incidentPion_avail_abs->Divide(h_sel_incidentPion_avail_abs, h_sel_pandora);
   h_eff_abs_evSel->Divide(h_sel_abs_evSel, h_sel_incidentPion_avail_abs);


   //--------------------------------------------------------
   // Save Histos
   //--------------------------------------------------------
   output->cd();
   h_trueAvailable->Write(); h_sel_pandora->Write(); h_sel_incidentPion_avail_abs->Write(); h_sel_abs_evSel->Write();
   h_eff_pandora->Write(); h_eff_incidentPion_avail_abs->Write(); h_eff_abs_evSel->Write();

   //frame_branches.Snapshot("pionana/beamana", "test_snap.root");

   return 0;
}

