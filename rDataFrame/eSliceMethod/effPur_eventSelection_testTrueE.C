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

int effPur_eventSelection_testTrueE(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);


   //output file
   TFile *output = new TFile ("eSliceMethod_effPur_binByBin_trueE.root", "RECREATE");

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
   //all available True Signals
   auto all_avail_abs = frame.Filter("true_absSignal");

   auto all_avail_cex = frame.Filter("true_chexSignal");

   auto all_avail_totInel = frame.Filter("true_primPionInel");

   //Pandora Recognised
   auto pandora_all = frame.Filter("primary_isBeamType");

   auto pandora_avail_abs = frame.Filter("primary_isBeamType && true_absSignal");

   auto pandora_avail_cex = frame.Filter("primary_isBeamType && true_chexSignal");

   auto pandora_avail_totinel = frame.Filter("primary_isBeamType && passBeamCut && passBeamCutBI && true_primPionInel");

   //After incident Pion Selection (Pandora && BeamCuts && muon removal)
   auto incidentPion_all = frame.Filter("selected_incidentPion");

   auto incidentPion_avail_abs = frame.Filter("selected_incidentPion && true_absSignal");

   auto incidentPion_avail_cex = frame.Filter("selected_incidentPion && true_chexSignal");

   auto incidentPion_avail_totInel = frame.Filter("selected_incidentPion && true_primPionInel");

   //Final selected Events
   auto selected_abs = frame.Filter("selected_abs");

   auto abs_select_avail_abs = frame.Filter("selected_abs && true_absSignal");

   auto selected_cex = frame.Filter("selected_cex");

   auto cex_select_avail_cex = frame.Filter("selected_cex && true_chexSignal");

   auto selected_totInel = frame.Filter("primary_ends_inAPA3 && selected_incidentPion");

   auto totInel_select_avail_totInel = frame
      .Filter("primary_ends_inAPA3 && selected_incidentPion && true_primPionInel");

   //--------------------------------------------------------
   //
   //For the interacting Histogram Build the efficiency and purity bin-by-bin histo
   //SELECTION EFF AND PUR
   //
   //Eff = N True Selected Abs / N True Available Abs in incident Selection
   //Pur = N True Selected Abs / N Selected Abs
   //
   //--------------------------------------------------------
   //             ABSORPTION
   //
   //--------------------------------------------------------


   TH1D* h_evSel_true_abs = new TH1D("h_evSel_true_abs", "Selected Abs, trueAbs", nBin_int, eEnd, eStart);
   TH1D* h_incidentPion_avail_abs = new TH1D("h_incidentPion_avail_abs", "Available Abs after incident Pion, trueAbs", nBin_int, eEnd, eStart);
   TH1D* h_evSel_selected_abs = new TH1D("h_evSel_selected_abs", "Selected Abs", nBin_int, eEnd, eStart);

   TH1D* h_eff_interacting_abs = new TH1D("h_eff_interacting_abs", "Efficiency Interacting Abs Histo", nBin_int, eEnd, eStart);
   TH1D* h_pur_interacting_abs = new TH1D("h_pur_interacting_abs", "Purity Interacting Abs Histo", nBin_int, eEnd, eStart);


   incidentPion_avail_abs
      .Foreach( [h_incidentPion_avail_abs] (double reco_beam_interactingEnergy){
            h_incidentPion_avail_abs->Fill(reco_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});

   abs_select_avail_abs
      .Foreach( [h_evSel_true_abs] (double reco_beam_interactingEnergy){
            h_evSel_true_abs->Fill(reco_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});

   selected_abs
      .Foreach( [h_evSel_selected_abs] (double reco_beam_interactingEnergy){
            h_evSel_selected_abs->Fill(reco_beam_interactingEnergy);}
            , {"true_KEint_fromEndP"} );

   h_eff_interacting_abs->Divide( h_evSel_true_abs, h_incidentPion_avail_abs );
   h_pur_interacting_abs->Divide( h_evSel_true_abs, h_evSel_selected_abs );

   output->cd();
   h_eff_interacting_abs->Write();
   h_pur_interacting_abs->Write();

   //--------------------------------------------------------
   //             CEX
   //
   //--------------------------------------------------------

   TH1D* h_evSel_true_cex = new TH1D("h_evSel_true_cex", "Selected cex, truecex", nBin_int, eEnd, eStart);
   TH1D* h_incidentPion_avail_cex = new TH1D("h_incidentPion_avail_cex", "Available cex after incident Pion, truecex", nBin_int, eEnd, eStart);
   TH1D* h_evSel_selected_cex = new TH1D("h_evSel_selected_cex", "Selected cex", nBin_int, eEnd, eStart);

   TH1D* h_eff_interacting_cex = new TH1D("h_eff_interacting_cex", "Efficiency Interacting cex Histo", nBin_int, eEnd, eStart);
   TH1D* h_pur_interacting_cex = new TH1D("h_pur_interacting_cex", "Purity Interacting cex Histo", nBin_int, eEnd, eStart);


   incidentPion_avail_cex
      .Foreach( [h_incidentPion_avail_cex] (double reco_beam_interactingEnergy){
            h_incidentPion_avail_cex->Fill(reco_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});

   cex_select_avail_cex
      .Foreach( [h_evSel_true_cex] (double reco_beam_interactingEnergy){
            h_evSel_true_cex->Fill(reco_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});

   selected_cex
      .Foreach( [h_evSel_selected_cex] (double reco_beam_interactingEnergy){
            h_evSel_selected_cex->Fill(reco_beam_interactingEnergy);}
            , {"true_KEint_fromEndP"} );

   h_eff_interacting_cex->Divide( h_evSel_true_cex, h_incidentPion_avail_cex );
   h_pur_interacting_cex->Divide( h_evSel_true_cex, h_evSel_selected_cex );

   output->cd();
   h_eff_interacting_cex->Write();
   h_pur_interacting_cex->Write();

   //--------------------------------------------------------
   //             totInel
   //
   //--------------------------------------------------------

   TH1D* h_evSel_true_totInel = new TH1D("h_evSel_true_totInel", "Selected totInel, truetotInel", nBin_int, eEnd, eStart);
   TH1D* h_incidentPion_avail_totInel = new TH1D("h_incidentPion_avail_totInel", "Available totInel after incident Pion, truetotInel", nBin_int, eEnd, eStart);
   TH1D* h_evSel_selected_totInel = new TH1D("h_evSel_selected_totInel", "Selected totInel", nBin_int, eEnd, eStart);

   TH1D* h_eff_interacting_totInel = new TH1D("h_eff_interacting_totInel", "Efficiency Interacting totInel Histo", nBin_int, eEnd, eStart);
   TH1D* h_pur_interacting_totInel = new TH1D("h_pur_interacting_totInel", "Purity Interacting totInel Histo", nBin_int, eEnd, eStart);


   incidentPion_avail_totInel
      .Foreach( [h_incidentPion_avail_totInel] (double reco_beam_interactingEnergy){
            h_incidentPion_avail_totInel->Fill(reco_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});

   totInel_select_avail_totInel
      .Foreach( [h_evSel_true_totInel] (double reco_beam_interactingEnergy){
            h_evSel_true_totInel->Fill(reco_beam_interactingEnergy);}
            ,{"true_KEint_fromEndP"});

   selected_totInel
      .Foreach( [h_evSel_selected_totInel] (double reco_beam_interactingEnergy){
            h_evSel_selected_totInel->Fill(reco_beam_interactingEnergy);}
            , {"true_KEint_fromEndP"} );

   h_eff_interacting_totInel->Divide( h_evSel_true_totInel, h_incidentPion_avail_totInel );
   h_pur_interacting_totInel->Divide( h_evSel_true_totInel, h_evSel_selected_totInel );

   output->cd();
   h_eff_interacting_totInel->Write();
   h_pur_interacting_totInel->Write();


   //--------------------------------------------------------
   //
   //For the incident Histogram Build the efficiency and purity bin-by-bin histo
   //SELECTION EFF AND PUR
   //
   //Eff = N True Selected IncPi / N True Available IncPi after Pandora recognition
   //Pur = N True Selected IncPi / N Selected incPi
   //
   //--------------------------------------------------------
   //
   TH1D* h_incidentPion_true_incPi = new TH1D("h_incidentPion_true_incPi", "Selected Incident Pions, true incident Pions", nBin_inc, eEnd, eStart);
   TH1D* h_selected_incidentPion = new TH1D("h_selected_incidentPion", "Selected Incident Pions", nBin_inc, eEnd, eStart);
   //true inc pi available that enter detector and are reconstructed by Pandore
   TH1D* h_available_true_incPi = new TH1D("h_available_true_incPi", "Available true incident Pions after Pandora Reconstruction", nBin_inc, eEnd, eStart);

   TH1D* h_eff_incidentPion = new TH1D("h_eff_incidentPion", "Efficiency incident Abs Histo", nBin_int, eEnd, eStart);
   TH1D* h_pur_incidentPion = new TH1D("h_pur_incidentPion", "Purity incident Abs Histo", nBin_int, eEnd, eStart);


   incidentPion_avail_totInel //N true incident pions after incident pion selection
      .Foreach( [h_incidentPion_true_incPi] (double true_firstEntryIncident, double reco_beam_interactingEnergy) {
            //bin that energy falls into is (int) energy/nbins + 1
            int binNumber_initEnergy = (int) true_firstEntryIncident / bin_size_inc + 1;
            int binNumber_interEnergy = (int) reco_beam_interactingEnergy / bin_size_inc + 1;
            //if(binNumber_initEnergy < 0 || binNumber_interEnergy < 0) return;
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
            if( i > 0 && i <= nBin_inc){ //make sure we don't go outside of bin range
            h_incidentPion_true_incPi->SetBinContent( i, h_incidentPion_true_incPi->GetBinContent(i) + 1 ); 
            //increment previous bin content by one. should handle number of entries in hist. unlike AddBinContent
            };
            };
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP"});

   incidentPion_all //N selected incident Pions
      .Foreach( [h_selected_incidentPion] (double true_firstEntryIncident, double reco_beam_interactingEnergy) {
            //bin that energy falls into is (int) energy/nbins + 1
            int binNumber_initEnergy = (int) true_firstEntryIncident / bin_size_inc + 1;
            int binNumber_interEnergy = (int) reco_beam_interactingEnergy / bin_size_inc + 1;
            //if(binNumber_initEnergy < 0 || binNumber_interEnergy < 0) return;
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
            if( i > 0 && i <= nBin_inc){ //make sure we don't go outside of bin range
            h_selected_incidentPion->SetBinContent( i, h_selected_incidentPion->GetBinContent(i) + 1 ); 
            //increment previous bin content by one. should handle number of entries in hist. unlike AddBinContent
            };
            };
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP"});

   pandora_avail_totinel //N true incident pions available after BeamCuts & Pandora
      .Foreach( [h_available_true_incPi] (double true_firstEntryIncident, double reco_beam_interactingEnergy) {
            //bin that energy falls into is (int) energy/nbins + 1
            int binNumber_initEnergy = (int) true_firstEntryIncident / bin_size_inc + 1;
            int binNumber_interEnergy = (int) reco_beam_interactingEnergy / bin_size_inc + 1;
            //if(binNumber_initEnergy < 0 || binNumber_interEnergy < 0) return;
            for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
            if( i > 0 && i <= nBin_inc){ //make sure we don't go outside of bin range
            h_available_true_incPi->SetBinContent( i, h_available_true_incPi->GetBinContent(i) + 1 ); 
            //increment previous bin content by one. should handle number of entries in hist. unlike AddBinContent
            };
            };
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP"});

   h_incidentPion_true_incPi->Sumw2(0);
   h_incidentPion_true_incPi->Rebin( bin_size_int/bin_size_inc );
   h_incidentPion_true_incPi->Scale( 1 / (bin_size_int/bin_size_inc) );
   h_incidentPion_true_incPi->Sumw2(0);
   
   h_selected_incidentPion->Sumw2(0);
   h_selected_incidentPion->Rebin( bin_size_int/bin_size_inc );
   h_selected_incidentPion->Scale( 1 / (bin_size_int/bin_size_inc) );
   h_selected_incidentPion->Sumw2(0);
  
   h_available_true_incPi->Sumw2(0);
   h_available_true_incPi->Rebin( bin_size_int/bin_size_inc );
   h_available_true_incPi->Scale( 1 / (bin_size_int/bin_size_inc) );
   h_available_true_incPi->Sumw2(0);


   h_pur_incidentPion->Divide( h_incidentPion_true_incPi, h_selected_incidentPion );
   h_eff_incidentPion->Divide( h_incidentPion_true_incPi, h_available_true_incPi );

   h_pur_incidentPion->Write();
   h_eff_incidentPion->Write();

   //--------------------------------------------------------
   // Save Histos
   //--------------------------------------------------------
   
   //frame_branches.Snapshot("pionana/beamana", "test_snap.root");

   return 0;
}

