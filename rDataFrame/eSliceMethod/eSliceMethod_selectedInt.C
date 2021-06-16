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

int eSliceMethod_selectedInt(const string mcFilepath, bool isMC, bool doEffPur, bool doUnsmear, bool doRecoEff){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);

   //output file
   string outputNameData = "output_eSliceMethod_selectedEvents_DATA_effPur_francesco.root";
   string outputNameMC = "output_eSliceMethod_selectedEvents_MC_effPur_unsmear_halfMCmatrix.root";
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
   string unsmearFile;
   if(isMC) unsmearFile = "unsmear_" + std::to_string((int) bin_size_int) + "MeV.root";
   else unsmearFile = "unsmear_" + std::to_string((int) bin_size_int) + "MeV.root";
   TFile f2(unsmearFile.c_str());

   TH1D *h_eff_eventSel_int = (TH1D*)f2.Get("h_eff_eventSel_int");
   TH1D *h_pur_removeBG_int = (TH1D*)f2.Get("h_pur_removeBG_int");
   TH1D *h_trueE_eff_reconstruction_int = (TH1D*)f2.Get("h_trueE_eff_reconstruction_int");
   TH2D *h2_inverse_smearing_interacting = (TH2D*)f2.Get("h2_inverse_smearing_interacting");

   TH1D *h_eff_eventSel_inc_initE = (TH1D*)f2.Get("h_eff_eventSel_inc_initE");
   TH1D *h_pur_removeBG_inc_initE = (TH1D*)f2.Get("h_pur_removeBG_inc_initE");
   TH1D *h_trueE_eff_reconstruction_inc_initE = (TH1D*)f2.Get("h_trueE_eff_reconstruction_inc_initE");
   TH2D *h2_inverse_smearing_incident_initE = (TH2D*)f2.Get("h2_inverse_smearing_incident_initE");

   TH1D *h_eff_eventSel_inc_interE = (TH1D*)f2.Get("h_eff_eventSel_inc_interE");
   TH1D *h_pur_removeBG_inc_interE = (TH1D*)f2.Get("h_pur_removeBG_inc_interE");
   TH1D *h_trueE_eff_reconstruction_inc_interE = (TH1D*)f2.Get("h_trueE_eff_reconstruction_inc_interE");
   TH2D *h2_inverse_smearing_incident_interE = (TH2D*)f2.Get("h2_inverse_smearing_incident_interE");

   //Get Reco Efficiencies to apply on trueE

   //Eff & Pur on TrueE bins
   TH1D *h_trueE_eff_eventSel_int = (TH1D*)f2.Get("h_trueE_eff_eventSel_int");
   TH1D *h_trueE_pur_removeBG_int = (TH1D*)f2.Get("h_trueE_pur_removeBG_int");

   TH1D *h_trueE_eff_eventSel_inc_initE = (TH1D*)f2.Get("h_trueE_eff_eventSel_inc_initE");
   TH1D *h_trueE_pur_removeBG_inc_initE = (TH1D*)f2.Get("h_trueE_pur_removeBG_inc_initE");

   TH1D *h_trueE_eff_eventSel_inc_interE = (TH1D*)f2.Get("h_trueE_eff_eventSel_inc_interE");
   TH1D *h_trueE_pur_removeBG_inc_interE = (TH1D*)f2.Get("h_trueE_pur_removeBG_inc_interE");

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
   // may 21
   // build incident histo now in a different way after having unsmeared start and end distributions
   //
   //--------------------------------------------------------

   //Incident Histogram and Interacting Histogram, interacting for different Processes
   //Energies are kinetic Energy of beam particle

   TH1D* h_recoE_selPion_inc_initE = new TH1D("h_recoE_selPion_inc_initE", "Incident Selected Pion reco_initE", nBin_int, eEnd, eStart);
   h_recoE_selPion_inc_initE->GetXaxis()->SetTitle("Reco KE (MeV)");

   TH1D* h_trueE_selPion_inc_initE = new TH1D("h_trueE_selPion_inc_initE", "Incident Selected Pion true_initE", nBin_int, eEnd, eStart);
   h_trueE_selPion_inc_initE->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_recoE_selPion_inc_interE = new TH1D("h_recoE_selPion_inc_interE", "Incident Selected Pion reco_interE", nBin_int, eEnd, eStart);
   h_recoE_selPion_inc_interE->GetXaxis()->SetTitle("Reco KE (MeV)");

   TH1D* h_trueE_selPion_inc_interE = new TH1D("h_trueE_selPion_inc_interE", "Incident Selected Pion true_interE", nBin_int, eEnd, eStart);
   h_trueE_selPion_inc_interE->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_trueE_selPion_incident = new TH1D("h_trueE_selPion_incident", "Incident Selected Pion trueE", nBin_int, eEnd, eStart);
   h_trueE_selPion_incident->GetXaxis()->SetTitle("True KE (MeV)");



   TH1D* h_recoE_selAbs_interacting = new TH1D("h_recoE_selAbs_interacting", "Interacting Selected ABS Reco Energy", nBin_int, eEnd, eStart);
   h_recoE_selAbs_interacting->GetXaxis()->SetTitle("Reco KE (MeV)");

   TH1D* h_trueE_selAbs_interacting = new TH1D("h_trueE_selAbs_interacting", "Interacting Selected ABS True Energy", nBin_int, eEnd, eStart);
   h_trueE_selAbs_interacting->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_recoE_unsmeared_inc_initE = new TH1D("h_recoE_unsmeared_inc_initE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_unsmeared_inc_interE = new TH1D("h_recoE_unsmeared_inc_interE", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_unsmeared_incident = new TH1D("h_recoE_unsmeared_incident", "", nBin_int, eEnd, eStart);
   TH1D* h_recoE_unsmeared_interacting = new TH1D("h_recoE_unsmeared_interacting", "", nBin_int, eEnd, eStart);

   //Initial Filters for all events
   auto mcIncident_selected_primaryPi = frame
      .Define("true_initKE", "true_firstEntryIncident")
      .Define("true_interKE", "true_interactingKE_fromLength")
      .Filter("selected_incidentPion");

   //========================================================
   //Build the Incident Histogram
   //---------
   //take the interacting energy, that slice is filled too
   if(isMC){
      mcIncident_selected_primaryPi
         .Foreach( [h_trueE_selPion_inc_initE, h_trueE_selPion_inc_interE] (double init_KE, double inter_KE) { 
               //make sure incident Pion does not interact in bin it was born
               int bin_initE = (int) init_KE / bin_size_inc + 1;
               int bin_interE = (int) inter_KE / bin_size_inc + 1;
               if( checkBins(init_KE, inter_KE, bin_initE, bin_interE) ){

               h_trueE_selPion_inc_initE->SetBinContent( bin_initE, h_trueE_selPion_inc_initE->GetBinContent( bin_initE ) + 1); 
               h_trueE_selPion_inc_interE->SetBinContent( bin_interE, h_trueE_selPion_inc_interE->GetBinContent( bin_interE ) + 1); 

               };
               }
               ,{"true_initKE", "true_interKE"});
   }

   mcIncident_selected_primaryPi
      .Foreach( [h_recoE_selPion_inc_initE, h_recoE_selPion_inc_interE] (double init_KE, double inter_KE) { 
            //make sure incident Pion does not interact in bin it was born
            int bin_initE = (int) init_KE / bin_size_inc + 1;
            int bin_interE = (int) inter_KE / bin_size_inc + 1;
            if( checkBins(init_KE, inter_KE, bin_initE, bin_interE) ){

            h_recoE_selPion_inc_initE->SetBinContent( bin_initE,  h_recoE_selPion_inc_initE->GetBinContent( bin_initE ) + 1); 
            h_recoE_selPion_inc_interE->SetBinContent( bin_interE, h_recoE_selPion_inc_interE->GetBinContent( bin_interE ) + 1); 

            };
            }
            ,{"reco_firstEntryIncident", "reco_interactingKE"});


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

         .Foreach( [h_trueE_selAbs_interacting] (double init_KE, double inter_KE) { 
               //make sure incident Pion does not interact in bin it was born
               int bin_initE = (int) init_KE / bin_size_inc + 1;
               int bin_interE = (int) inter_KE / bin_size_inc + 1;
               if( checkBins(init_KE, inter_KE, bin_initE, bin_interE) ){

               h_trueE_selAbs_interacting->SetBinContent( bin_interE, h_trueE_selAbs_interacting->GetBinContent( bin_interE ) + 1); 

               };
               }
               ,{"true_initKE", "true_interKE"});
      h_trueE_selAbs_interacting->Sumw2(0);
      h_trueE_selAbs_interacting->Write();
   }

   mcInteracting_selected_abs 
      .Foreach( [h_recoE_selAbs_interacting] (double init_KE, double inter_KE) { 
         //make sure incident Pion does not interact in bin it was born
         int bin_initE = (int) init_KE / bin_size_inc + 1;
         int bin_interE = (int) inter_KE / bin_size_inc + 1;
         if( checkBins(init_KE, inter_KE, bin_initE, bin_interE) ){

            h_recoE_selAbs_interacting->SetBinContent( bin_interE, h_recoE_selAbs_interacting->GetBinContent( bin_interE ) + 1); 

         };
      }
   ,{"reco_firstEntryIncident","reco_interactingKE"});

   h_recoE_selAbs_interacting->Sumw2(0);
   h_recoE_selAbs_interacting->Write();


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
   //   Unsmear Start and End Distributions for Incident Histogram
   //
   //------------------------------------------------------
   if(doEffPur){
      h_recoE_selPion_inc_initE->Multiply( h_pur_removeBG_inc_initE );
      h_recoE_selPion_inc_initE->Divide( h_eff_eventSel_inc_initE );

      //interE
      h_recoE_selPion_inc_interE->Multiply( h_pur_removeBG_inc_interE );
      h_recoE_selPion_inc_interE->Divide( h_eff_eventSel_inc_interE );


   };
   //Unsmear the start Pion Distribution with the inverse matrix
   //Loop through smearing matrix Incident
   if(doUnsmear){
      for(int i = 1; i <= nBin_int; i++){
         double help_sum_initE = 0, help_sum_interE=0;
         //rows
         for(int j = 1; j <= nBin_int; j++){
            help_sum_initE += h_recoE_selPion_inc_initE->GetBinContent(j)*h2_inverse_smearing_incident_initE->GetBinContent( i, j);
            help_sum_interE += h_recoE_selPion_inc_interE->GetBinContent(j)*h2_inverse_smearing_incident_interE->GetBinContent( i, j);
         }; 
         h_recoE_unsmeared_inc_initE->SetBinContent( i , help_sum_initE);
         h_recoE_unsmeared_inc_interE->SetBinContent( i , help_sum_interE);
      };

      h_recoE_unsmeared_inc_initE->Sumw2(0);
      h_recoE_unsmeared_inc_interE->Sumw2(0);
   }

   else {
      h_recoE_unsmeared_inc_initE = h_recoE_selPion_inc_initE;
      h_recoE_unsmeared_inc_interE = h_recoE_selPion_inc_interE;

   }
   if(doRecoEff){
      h_recoE_unsmeared_inc_initE->Divide( h_trueE_eff_reconstruction_inc_initE );
      h_recoE_unsmeared_inc_interE->Divide( h_trueE_eff_reconstruction_inc_interE );
   };

   h_recoE_unsmeared_inc_initE->Sumw2(0);
   h_recoE_unsmeared_inc_initE->Write();
   h_recoE_unsmeared_inc_interE->Sumw2(0);
   h_recoE_unsmeared_inc_interE->Write();
   //------------------------------------------------------
   //    Build recoE_unsmeared Incident Histo
   //    from the distributions for XS calculation after
   //------------------------------------------------------


   build_incidentHist( h_recoE_unsmeared_inc_initE, h_recoE_unsmeared_inc_interE, h_recoE_unsmeared_incident );
   h_recoE_unsmeared_incident->Sumw2(0);
   h_recoE_unsmeared_incident->Write();

   //------------------------------------------------------
   //   Unsmear Interacting
   //
   //------------------------------------------------------

   //Take account for purity of TrueE == RecoE and Efficiency from Smearing
   if(doEffPur){
      h_recoE_selAbs_interacting->Multiply( h_pur_removeBG_int);
      h_recoE_selAbs_interacting->Divide( h_eff_eventSel_int );

   };
   if(doUnsmear){   for(int i = 1; i <= nBin_int; i++){
      double help_sum = 0;
      //rows
      for(int j = 1; j <= nBin_int; j++){
         help_sum += h_recoE_selAbs_interacting->GetBinContent(j)*h2_inverse_smearing_interacting->GetBinContent( i, j);
      }; 
      h_recoE_unsmeared_interacting->SetBinContent( i , help_sum);
   };

   h_recoE_unsmeared_interacting->Sumw2(0);
   }

   else {
      h_recoE_unsmeared_interacting = h_recoE_selAbs_interacting;
   }
   if(doRecoEff){
      h_recoE_unsmeared_interacting->Divide( h_trueE_eff_reconstruction_int );

   };
   h_recoE_unsmeared_interacting->Sumw2(0);
   h_recoE_unsmeared_interacting->Write();
   //=====================================================
   //             Computing the XS
   //=====================================================
   //
   // xs(Ebin) = (A / (Na*density*bin_size)) * dEdX(Ebin) * hInteracting / hIncident
   //
   // More Accurate use log( Ninc / (Ninc - Nint ))

   //------------------------------------------------------
   //    Absorption, Selected Interactions Reconstrucetd Energy
   //------------------------------------------------------
   //
   string xs_abs_name, xs_abs_title, xs_totInel_name, xs_totInel_title;
   if(isMC) {
      xs_totInel_name = "h_xs_recoE_selected_totInel_mc";
      xs_totInel_title = "Total Inelastic MC";
      xs_abs_name = "h_xs_recoE_selected_abs_mc";
      xs_abs_title = "Absorption MC";
   }
   else {
      xs_totInel_name = "h_xs_recoE_selected_totInel_data";
      xs_totInel_title = "Total Inelastic Data";
      xs_abs_name = "h_xs_recoE_selected_abs_data";
      xs_abs_title = "Absorption Data";
   }


   TH1D* h_xs_recoE_selected_abs = new TH1D( xs_abs_name.c_str(), xs_abs_title.c_str(), nBin_int, eEnd, eStart);
   TH1D* h_xs_recoE_selected_totInel = new TH1D( xs_totInel_name.c_str(), xs_totInel_title.c_str(), nBin_int, eEnd, eStart);

   //Function doing XS
   do_XS_log( h_xs_recoE_selected_abs, h_recoE_unsmeared_interacting, h_recoE_unsmeared_incident, h_betheMean_muon );   
   do_XS_log_binomial_error( h_xs_recoE_selected_abs, h_recoE_unsmeared_interacting, h_recoE_unsmeared_incident, h_betheMean_muon );   

   do_XS_log( h_xs_recoE_selected_totInel, h_recoE_unsmeared_inc_interE, h_recoE_unsmeared_incident, h_betheMean_muon );   
   do_XS_log_binomial_error( h_xs_recoE_selected_totInel, h_recoE_unsmeared_inc_interE, h_recoE_unsmeared_incident, h_betheMean_muon );   

   
   h_xs_recoE_selected_abs->Write();
   h_xs_recoE_selected_totInel->Write();

   //------------------------------------------------------
   //     For TrueE  
   //------------------------------------------------------
   h_trueE_selPion_inc_initE->Multiply( h_trueE_pur_removeBG_inc_initE );
   h_trueE_selPion_inc_initE->Divide( h_trueE_eff_eventSel_inc_initE );
   //also reco efficiency
   h_trueE_selPion_inc_initE->Divide( h_trueE_eff_reconstruction_inc_initE );

   h_trueE_selPion_inc_interE->Multiply( h_trueE_pur_removeBG_inc_interE );
   h_trueE_selPion_inc_interE->Divide( h_trueE_eff_eventSel_inc_interE );
   //also reco efficiency
   h_trueE_selPion_inc_interE->Divide( h_trueE_eff_reconstruction_inc_interE );

   h_trueE_selAbs_interacting->Multiply( h_trueE_pur_removeBG_int);
   h_trueE_selAbs_interacting->Divide( h_trueE_eff_eventSel_int );
   //also reco efficiency
   h_trueE_selAbs_interacting->Divide( h_trueE_eff_reconstruction_int);


   build_incidentHist( h_trueE_selPion_inc_initE, h_trueE_selPion_inc_interE, h_trueE_selPion_incident );
   h_trueE_selPion_incident->Sumw2(0);
   h_trueE_selPion_incident->Write();
   //-----------

   TH1D* h_xs_trueE_selected_abs = new TH1D( "h_xs_trueE_selected_abs", "true Kinetic Energy (MeV)", nBin_int, eEnd, eStart );
   TH1D* h_xs_trueE_selected_totInel = new TH1D( "h_xs_trueE_selected_totInel", "true Kinetic Energy (MeV)", nBin_int, eEnd, eStart );
   //Function doing XS
   if(isMC){   

      do_XS_log( h_xs_trueE_selected_abs, h_trueE_selAbs_interacting, h_trueE_selPion_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_xs_trueE_selected_abs, h_trueE_selAbs_interacting, h_trueE_selPion_incident, h_betheMean_muon );   
      
      do_XS_log( h_xs_trueE_selected_totInel, h_trueE_selPion_inc_interE, h_trueE_selPion_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_xs_trueE_selected_totInel, h_trueE_selPion_inc_interE, h_trueE_selPion_incident, h_betheMean_muon );   

   }

   h_xs_trueE_selected_abs->Write();
   h_xs_trueE_selected_totInel->Write();

   //------------------------------------------------------

   //=====================================================*
   //            Plotting and Style
   //=====================================================
   //
   //h_xs_recoE_selected_abs->Scale(2);

   TCanvas *c_recoE_abs = new TCanvas("c_recoE_abs", "c_recoE_abs");
   gPad->SetGrid(1,1);
   h_xs_recoE_selected_abs->SetTitle( "Selected Absorption; Reco Kinetic Energy (MeV); #sigma (mb)");
   h_xs_recoE_selected_abs->GetXaxis()->SetRangeUser(300,1000);
   h_xs_recoE_selected_abs->GetXaxis()->SetNdivisions(1020);
   h_xs_recoE_selected_abs->GetYaxis()->SetNdivisions(1020);

   string c_abs_title;
   if(isMC) c_abs_title = "Absorption MC; Reco kinetic Energy (MeV); #sigma (mbarn)";
   else c_abs_title = "Absorption Data; Reco kinetic Energy (MeV); #sigma (mbarn)";

   abs_KE->SetTitle( c_abs_title.c_str());
   abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   abs_KE->SetLineColor(kBlue);
   abs_KE->SetLineWidth(3);
   abs_KE->Draw("AC");
   h_xs_recoE_selected_abs->SetMarkerSize(0.7);
   h_xs_recoE_selected_abs->Draw("PE0 SAME");

   c_recoE_abs->Write();

   if(isMC){
      TCanvas *c_trueE_abs = new TCanvas("c_trueE_abs", "c_trueE_abs");
      gPad->SetGrid(1,1);
      h_xs_trueE_selected_abs->SetTitle( "Selected Absorption;true Kinetic Energy (MeV); #sigma (mb)");
      h_xs_trueE_selected_abs->GetXaxis()->SetRangeUser(300,1000);
      h_xs_trueE_selected_abs->GetXaxis()->SetNdivisions(1020);
      h_xs_trueE_selected_abs->GetYaxis()->SetNdivisions(1020);

      abs_KE->SetTitle( "Absorption;true Kinetic Energy (MeV); #sigma (mb)");
      abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
      abs_KE->SetLineColor(kRed);
      abs_KE->SetLineWidth(3);
      abs_KE->Draw("AC");
      h_xs_trueE_selected_abs->SetMarkerSize(0.7);
      h_xs_trueE_selected_abs->Draw("PE0 SAME");

      c_trueE_abs->Write();
   };


   TCanvas *c_recoE_totInel = new TCanvas("c_recoE_totInel", "c_recoE_totInel");
   gPad->SetGrid(1,1);
   h_xs_recoE_selected_totInel->SetTitle( "Selected Total Inelastic; Reco Kinetic Energy (MeV); #sigma (mb)");
   h_xs_recoE_selected_totInel->GetXaxis()->SetRangeUser(300,1000);
   h_xs_recoE_selected_totInel->GetXaxis()->SetNdivisions(1020);
   h_xs_recoE_selected_totInel->GetYaxis()->SetNdivisions(1020);

   string c_totInel_title;
   if(isMC) c_totInel_title = "Total Inelastic MC; Reco kinetic Energy (MeV); #sigma (mbarn)";
   else c_totInel_title = "Total Inelastic Data; Reco kinetic Energy (MeV); #sigma (mbarn)";

   totInel_KE->SetTitle( c_totInel_title.c_str());
   totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   totInel_KE->SetLineColor(kRed);
   totInel_KE->SetLineWidth(3);
   totInel_KE->Draw("AC");
   h_xs_recoE_selected_totInel->SetMarkerSize(0.7);
   h_xs_recoE_selected_totInel->Draw("PE0 SAME");

   c_recoE_totInel->Write();

   if(isMC){
      TCanvas *c_trueE_totInel = new TCanvas("c_trueE_totInel", "c_trueE_totInel");
      gPad->SetGrid(1,1);
      h_xs_trueE_selected_totInel->SetTitle( "Selected Total Inelastic; true Kinetic Energy (MeV); #sigma (mb)");
      h_xs_trueE_selected_totInel->GetXaxis()->SetRangeUser(300,1000);
      h_xs_trueE_selected_totInel->GetXaxis()->SetNdivisions(1020);
      h_xs_trueE_selected_totInel->GetYaxis()->SetNdivisions(1020);


      totInel_KE->SetTitle( c_totInel_title.c_str());
      totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
      totInel_KE->SetLineColor(kRed);
      totInel_KE->SetLineWidth(3);
      totInel_KE->Draw("AC");
      h_xs_trueE_selected_totInel->SetMarkerSize(0.7);
      h_xs_trueE_selected_totInel->Draw("PE0 SAME");

      c_trueE_totInel->Write();
   }
   
   //output->Write();
   //f1.Close();
   return 0;
}

