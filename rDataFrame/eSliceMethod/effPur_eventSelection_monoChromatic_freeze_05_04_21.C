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

int effPur_eventSelection_monoChromatic_freeze_05_04_21(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);


   //output file
   TFile *output = new TFile ("eSliceMethod_effPur_binByBin_monoChromatic_freeze_05_04_21.root", "RECREATE");

   //--------------------------------------------------------
   //
   //--------------------------------------------------------
   //

   //After incident Pion Selection (Pandora && BeamCuts && muon removal)
   auto selected_incidentPion = frame.Filter("selected_incidentPion");

   auto selected_incidentPion_truePion = frame.Filter("selected_incidentPion && true_primPionInel");

   auto selected_beamCut_truePion = frame.Filter("primary_isBeamType && passBeamCut && passBeamCutBI && true_primPionInel");

   //Final selected Events
   auto selected_abs = frame.Filter("selected_abs");

   auto selected_abs_trueAbs = frame.Filter("selected_abs && true_absSignal");

   //Available Absorptions at incidentPion Sample (--> need to also provide a selection efficiency!)
   auto selected_incidentPion_trueAbs = frame.Filter("selected_incidentPion && true_absSignal");

   //--------------------------------------------------------
   //
   //For the interacting Histogram Build the efficiency and purity bin-by-bin histo
   //SELECTION EFF AND PUR
   //
   //ASSUME The MONOCHROMATIC CASE where the smearing into other energy bins is not included 
   //--> stat killer, can be developed further to take into account the smeared signals and not consider them as lost
   //
   //Eff = Number of True Signal in RecoE == TrueE bin / Number of True Signal in TrueE bin
   //Pur = Number of True Signal in RecoE == TrueE bin / Number of Selected Abs Candidates in that RecoBin ( Signal + BG and BG includes also the smeared Abs from other E-bins)
   //
   //
   //--------------------------------------------------------
   //             ABSORPTION
   //
   //--------------------------------------------------------


   //only filling now for the selected where trueBin == recoBin
   
   //Number of True Signal in RecoE == TrueE bin
   TH1D* h_sel_trueAbs_trueE_equal_recoE = new TH1D("h_sel_trueAbs_trueE_equal_recoE", "", nBin_int, eEnd, eStart);
   
   //Number of True Signal in TrueE bin 
   TH1D* h_sel_trueAbs_trueE = new TH1D("h_sel_trueAbs_trueE", "", nBin_int, eEnd, eStart);
   
   //Number of Selected Abs Candidates in that RecoBin
   TH1D* h_sel_absCandidates_recoE = new TH1D("h_sel_absCandidates_recoE", "", nBin_int, eEnd, eStart);
   
   //Number of Available Absorption Events in Incident Sample
   TH1D* h_trueAbs_atIncident = new TH1D("h_trueAbs_atIncident", "", nBin_int, eEnd, eStart);

   TH1D* h_eff_selection_abs = new TH1D("h_eff_selection_abs", "Efficiency From Selection Abs Histo", nBin_int, eEnd, eStart);
   TH1D* h_eff_interacting_abs = new TH1D("h_eff_interacting_abs", "Efficiency Interacting Abs Histo", nBin_int, eEnd, eStart);
   TH1D* h_pur_interacting_abs = new TH1D("h_pur_interacting_abs", "Purity Interacting Abs Histo", nBin_int, eEnd, eStart);

   //check here that binReco and binTrue are the same for the interacting KE
   selected_abs_trueAbs
      .Foreach( [h_sel_trueAbs_trueE_equal_recoE] (double true_beam_interactingEnergy, double reco_beam_interactingEnergy){

            int binReco = (int) reco_beam_interactingEnergy / bin_size_int + 1;
            int binTrue = (int) true_beam_interactingEnergy / bin_size_int + 1;

            if (binReco == binTrue) h_sel_trueAbs_trueE_equal_recoE->Fill(reco_beam_interactingEnergy); }
            ,{"true_KEint_fromEndP", "reco_interactingKE"});


   selected_abs_trueAbs
      .Foreach( [h_sel_trueAbs_trueE] (double true_beam_interactingEnergy){

            h_sel_trueAbs_trueE->Fill(true_beam_interactingEnergy); }
            
            ,{"true_KEint_fromEndP"});
//test
   selected_incidentPion_trueAbs      
      .Foreach( [h_trueAbs_atIncident] (double true_beam_interactingEnergy){

            h_trueAbs_atIncident->Fill(true_beam_interactingEnergy); }
            
            ,{"true_KEint_fromEndP"});

    selected_abs
      .Foreach( [h_sel_absCandidates_recoE] (double reco_beam_interactingEnergy){
            h_sel_absCandidates_recoE->Fill(reco_beam_interactingEnergy);}
            
            , {"reco_interactingKE"} );

   //Efficiency due to Smearing within selected Absorption
   h_eff_interacting_abs->Divide( h_sel_trueAbs_trueE_equal_recoE, h_sel_trueAbs_trueE );
   
   //Selection Efficiency wrt trueAbs in Incident Sample
   h_eff_selection_abs->Divide( h_sel_trueAbs_trueE, h_trueAbs_atIncident );

   h_pur_interacting_abs->Divide( h_sel_trueAbs_trueE_equal_recoE, h_sel_absCandidates_recoE );


   h_sel_trueAbs_trueE_equal_recoE->Write();
   h_sel_trueAbs_trueE->Write();
   h_sel_absCandidates_recoE->Write();
   h_sel_trueAbs_trueE->Write();
   
   output->cd();
   h_eff_interacting_abs->Write();
   h_eff_selection_abs->Write();
   h_pur_interacting_abs->Write();

   //--------------------------------------------------------
   //--------------------------------------------------------
   //
   //For the incident Histogram Build the efficiency and purity bin-by-bin histo
   //SELECTION EFF AND PUR
   //
   //
   //ASSUME The MONOCHROMATIC CASE where the smearing into other energy bins is not included 
   //--> stat killer, can be developed further to take into account the smeared signals and not consider them as lost
   //
   //Eff = Number of True Signal in RecoE == TrueE bin / Number of True Signal in TrueE bin
   //Pur = Number of True Signal in RecoE == TrueE bin / Number of Selected IncPi Candidates in that RecoBin ( Signal + BG and BG includes also the smeared incPi from other E-bins)
   //
   //--------------------------------------------------------
   //
   //Number of True Signal in RecoE == TrueE bin
   TH1D* h_sel_truePi_trueE_equal_recoE = new TH1D("h_sel_truePi_trueE_equal_recoE", "", nBin_int, eEnd, eStart);
   //Number of True Signal in TrueE bin
   TH1D* h_sel_truePi_trueE = new TH1D("h_sel_truePi_trueE", "", nBin_int, eEnd, eStart);
   //Number of Selected IncPi Candidates in that RecoBin
   TH1D* h_sel_piCandidates_recoE = new TH1D("h_sel_piCandidates_recoE", "", nBin_int, eEnd, eStart);
   
   TH1D* h_truePi_atBeamCut = new TH1D("h_truePi_atBeamCut", "", nBin_int, eEnd, eStart);

   TH1D* h_eff_selection_incidentPion = new TH1D("h_eff_selection_incidentPion", "", nBin_int, eEnd, eStart);

   TH1D* h_eff_incidentPion = new TH1D("h_eff_incidentPion", "", nBin_int, eEnd, eStart);
   TH1D* h_pur_incidentPion = new TH1D("h_pur_incidentPion", "", nBin_int, eEnd, eStart);


   selected_incidentPion_truePion //N true incident pions after incident pion selection
      .Foreach( [h_sel_truePi_trueE_equal_recoE] (double true_first, double true_int, double reco_first, double reco_int) {
            //bin that energy falls into is (int) energy/nbins + 1
            //
            int true_binHigh_init = (int) true_first / bin_size_int + 1;
            int reco_binHigh_init = (int) reco_first / bin_size_int + 1;

            int true_binLow_inter = (int) true_int / bin_size_int + 1;
            int reco_binLow_inter = (int) reco_int / bin_size_int + 1;

            //compare reco bins and true bins, find overlapping area
            //for initial bin / high bin take the smaller of the two
            //for interacting bin / lower bin take the higher of the two
            int binNumber_initEnergy, binNumber_interEnergy;
            
            if(true_binHigh_init > reco_binHigh_init) binNumber_initEnergy = reco_binHigh_init;
            else binNumber_initEnergy = true_binHigh_init;

            if(true_binLow_inter < reco_binLow_inter) binNumber_interEnergy = reco_binLow_inter;
            else binNumber_interEnergy = true_binLow_inter;

                      for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
            if( i > 0 && i <= nBin_int){ //make sure we don't go outside of bin range
            h_sel_truePi_trueE_equal_recoE->SetBinContent( i, h_sel_truePi_trueE_equal_recoE->GetBinContent(i) + 1 ); 
                        };
            };
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP","reco_firstEntryIncident", "reco_interactingKE"});
   
   selected_incidentPion_truePion //N selected incident Pions
      .Foreach( [h_sel_truePi_trueE] (double true_firstEntryIncident, double true_beam_interactingEnergy) {
            //bin that energy falls into is (int) energy/nbins + 1
            int binNumber_initEnergy = (int) true_firstEntryIncident / bin_size_int + 1;
            int binNumber_interEnergy = (int) true_beam_interactingEnergy / bin_size_int + 1;
               for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
                  if( i > 0 && i <= nBin_int){
                     h_sel_truePi_trueE->SetBinContent( i, h_sel_truePi_trueE->GetBinContent(i) + 1 );
                  };
               };
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP"});

   selected_beamCut_truePion //N selected incident Pions
      .Foreach( [h_truePi_atBeamCut] (double true_firstEntryIncident, double true_beam_interactingEnergy) {
            //bin that energy falls into is (int) energy/nbins + 1
            int binNumber_initEnergy = (int) true_firstEntryIncident / bin_size_int + 1;
            int binNumber_interEnergy = (int) true_beam_interactingEnergy / bin_size_int + 1;
               for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
                  if( i > 0 && i <= nBin_int){
                     h_truePi_atBeamCut->SetBinContent( i, h_truePi_atBeamCut->GetBinContent(i) + 1 );
                  };
               };
            }
            ,{"true_firstEntryIncident", "true_KEint_fromEndP"});




   selected_incidentPion 
      .Foreach( [h_sel_piCandidates_recoE] (double reco_int, double reco_inc) {
            //bin that energy falls into is (int) energy/nbins + 1
            //
            int reco_bin_interacting = (int) reco_int / bin_size_int + 1;
            int reco_bin_incident = (int) reco_inc / bin_size_int + 1;

           
            for(int i = reco_bin_interacting; i <= reco_bin_incident; i++){      
            if( i > 0 && i <= nBin_int){ //make sure we don't go outside of bin range
               h_sel_piCandidates_recoE->SetBinContent( i, h_sel_piCandidates_recoE->GetBinContent(i) + 1 );
            };

            };
            
            }
            ,{"reco_interactingKE", "reco_firstEntryIncident"});


 
   h_eff_incidentPion->Divide( h_sel_truePi_trueE_equal_recoE, h_sel_truePi_trueE );

   //Selection Efficiency wrt trueIncPi after Beam Cuts sample
   h_eff_selection_incidentPion->Divide( h_sel_truePi_trueE, h_truePi_atBeamCut );

   h_pur_incidentPion->Divide( h_sel_truePi_trueE_equal_recoE, h_sel_piCandidates_recoE );

   h_sel_truePi_trueE_equal_recoE->Write();
   h_sel_truePi_trueE->Write();
   h_sel_piCandidates_recoE->Write();

   h_pur_incidentPion->Write();
   h_eff_incidentPion->Write();
   h_eff_selection_incidentPion->Write();

   //--------------------------------------------------------
   // Save Histos
   //--------------------------------------------------------

   //frame_branches.Snapshot("pionana/beamana", "test_snap.root");

   return 0;
}

