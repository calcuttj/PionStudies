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

int eSliceMethod_trueProcess_trueE(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);

   //output file
   string outputNameMC = "output_eSliceMethod_trueProcess_trueE.root";
   string outputName;
   outputName = outputNameMC;


   TFile *output = new TFile ( outputName.c_str() , "RECREATE");

   //access Jakes GeantFile in folder
   TFile f1("exclusive_xsec.root");
   TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
   TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
   TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
   f1.Close();

   string mg_title;
   mg_title = "Cross-Section MC; True kinetic Energy (MeV); #sigma (mbarn)";

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(totInel_KE);
   mg->Add(abs_KE);
   //mg->Add(cex_KE);
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
   TH1D* h_trueE_truePion_inc_initE = new TH1D("h_trueE_truePion_inc_initE", "Incident Selected Pion true_initE", nBin_int, eEnd, eStart);
   h_trueE_truePion_inc_initE->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_trueE_truePion_inc_interE = new TH1D("h_trueE_truePion_inc_interE", "Incident Selected Pion true_interE", nBin_int, eEnd, eStart);
   h_trueE_truePion_inc_interE->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_trueE_truePionInel = new TH1D("h_trueE_truePionInel", "Incident Selected Pion true_interE", nBin_int, eEnd, eStart);
   h_trueE_truePionInel->GetXaxis()->SetTitle("True KE (MeV)");
   
   TH1D* h_trueE_truePion_incident = new TH1D("h_trueE_truePion_incident", "Incident Selected Pion trueE", nBin_inc, eEnd, eStart);
   h_trueE_truePion_incident->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_trueE_trueAbs_interacting = new TH1D("h_trueE_trueAbs_interacting", "Interacting Selected ABS True Energy", nBin_int, eEnd, eStart);
   h_trueE_trueAbs_interacting->GetXaxis()->SetTitle("True KE (MeV)");

   //Initial Filters for all events
   auto mcIncident_true_primaryPi = frame
      .Filter("true_beam_endZ > 0")
      //.Filter("selected_incidentPion")
      //.Filter("primary_isBeamType && passBeamQuality_TPCjustPosition")
      .Define("true_initKE", "true_firstEntryIncident")
      .Define("true_interKE", "true_interactingKE_fromLength")
      .Filter("true_beam_PDG == 211");
      //.Filter("true_beam_PDG == 211 && passBeamQuality_TPCjustPosition && primary_isBeamType"); //shows how reco efficiency acts ununiformly on the XS computation

   //.Range(70,100)

   auto mcInteracting_true_abs = mcIncident_true_primaryPi
      .Filter("true_absSignal && true_pion_daughter == 0");
   //========================================================
   //Build the Incident Histogram
   //---------
   mcIncident_true_primaryPi
      .Foreach( [h_trueE_truePion_inc_initE, h_trueE_truePion_inc_interE] (double true_initKE, double true_beam_interactingEnergy) { 
            //make sure incident Pion does not interact in bin it was born
            int binNum_initE = (int) true_initKE / bin_size_inc + 1;
            int binNum_interE = (int) true_beam_interactingEnergy / bin_size_inc + 1;
            if( checkBins(true_initKE, true_beam_interactingEnergy, binNum_initE, binNum_interE) ){

               h_trueE_truePion_inc_initE->SetBinContent( binNum_initE, h_trueE_truePion_inc_initE->GetBinContent( binNum_initE ) + 1); 
               h_trueE_truePion_inc_interE->SetBinContent( binNum_interE, h_trueE_truePion_inc_interE->GetBinContent( binNum_interE ) + 1); 

            };   
            }
            ,{"true_initKE", "true_interKE"});

   h_trueE_truePion_inc_initE->Write();
   h_trueE_truePion_inc_interE->Write();


   //=====================================================
   //------------------------------------------------------
   //Interacting selected samples
   //------------------------------------------------------
   mcIncident_true_primaryPi
      .Filter("true_primPionInel")
      .Foreach( [h_trueE_truePionInel] (double true_initKE, double true_beam_interactingEnergy) { 
            //make sure incident Pion does not interact in bin it was born
            int binNum_initE = (int) true_initKE / bin_size_inc + 1;
            int binNum_interE = (int) true_beam_interactingEnergy / bin_size_inc + 1;
            if( checkBins(true_initKE, true_beam_interactingEnergy, binNum_initE, binNum_interE) ){

               h_trueE_truePionInel->SetBinContent( binNum_interE, h_trueE_truePionInel->GetBinContent( binNum_interE ) + 1); 

            };   
            }
            ,{"true_initKE", "true_interKE"});

   h_trueE_truePionInel->Write();


   mcInteracting_true_abs
      .Foreach( [h_trueE_trueAbs_interacting] ( double true_initKE, double true_beam_interactingEnergy ){

            int binNum_initE = (int) true_initKE / bin_size_inc + 1;
            int binNum_interE = (int) true_beam_interactingEnergy / bin_size_inc + 1;
            if( checkBins(true_initKE, true_beam_interactingEnergy, binNum_initE, binNum_interE) ){

            h_trueE_trueAbs_interacting->SetBinContent( binNum_interE, h_trueE_trueAbs_interacting->GetBinContent( binNum_interE ) + 1); 
    
            };
            }            
            ,{"true_initKE","true_interKE"});

   h_trueE_trueAbs_interacting->Sumw2(0);
   h_trueE_trueAbs_interacting->Write();

   //=====================================================
   //           Incident
   //=====================================================
   //Function in eSlice.h to build Incident Histogram

   
   build_incidentHist( h_trueE_truePion_inc_initE, h_trueE_truePion_inc_interE, h_trueE_truePion_incident );
   h_trueE_truePion_incident->Write();

   //for(int i = nBin_int; i >=1; i--) std::cout << "Entry[" << i << "] = " << h_trueE_truePion_incident->GetBinContent(i)  << std::endl;


   //=====================================================
   //            Prepare BetheBloch Mean for each Bin 
   //            QUESTION: take betheBloch of Pion or Muon?? Comparison to data fits better muon Bethe... 
   //            at hihger momentum ~400-800 they anway are almost the same
   //=====================================================
   TH1D* h_betheMean_muon = new TH1D("h_betheMean_muon", "Mean Energy Loss", nBin_int, eEnd, eStart);

   //fill histo with Mean dEdX of bin center
   fill_betheHisto( h_betheMean_muon, mass_pion );
   h_betheMean_muon->Write();


   //=====================================================
   //             Computing the XS
   //=====================================================
   //
   // xs(Ebin) = (A / (Na*density*bin_size)) * dEdX(Ebin) * hInteracting / hIncident
   //
   // More Accurate use log( Ninc / (Ninc - Nint ))
   //
   //scale_factor is done in eSlice.h

   //------------------------------------------------------
   //    Absorption, Selected Interactions Reconstrucetd Energy
   //------------------------------------------------------


   TH1D* h_xs_trueE_trueAbs = new TH1D("h_xs_trueE_trueAbs", "Absortpion MC", nBin_int, eEnd, eStart);

   //Function to do XS with log formula

   do_XS_log( h_xs_trueE_trueAbs, h_trueE_trueAbs_interacting, h_trueE_truePion_incident, h_betheMean_muon);
   
   do_XS_log_binomial_error( h_xs_trueE_trueAbs, h_trueE_trueAbs_interacting, h_trueE_truePion_incident, h_betheMean_muon);
   
   h_trueE_trueAbs_interacting->Write();
   h_xs_trueE_trueAbs->Write();


   //------------------------------------------------------
   //    Total Inelastici
   //------------------------------------------------------
   //

   TH1D* h_xs_trueE_truePion_totInel = new TH1D("h_xs_trueE_truePion_totInel", "TotInel MC", nBin_int, eEnd, eStart);
 
   do_XS_log( h_xs_trueE_truePion_totInel , h_trueE_truePionInel, h_trueE_truePion_incident, h_betheMean_muon);
   
   do_XS_log_binomial_error( h_xs_trueE_truePion_totInel, h_trueE_truePionInel, h_trueE_truePion_incident, h_betheMean_muon);  
   
   h_trueE_truePionInel->Write();
   h_xs_trueE_truePion_totInel->Write();

   //=====================================================*
   //            Plotting and Style
   //=====================================================
   //
   TCanvas *c_trueE_abs = new TCanvas("c_trueE_abs", "c_trueE_abs");
   gPad->SetGrid(1,1);
   h_xs_trueE_trueAbs->SetTitle( "True Absorption;true Kinetic Energy (MeV); #sigma (mb)");
   h_xs_trueE_trueAbs->GetXaxis()->SetRangeUser(100,1000);
   h_xs_trueE_trueAbs->GetXaxis()->SetNdivisions(1020);
   h_xs_trueE_trueAbs->GetYaxis()->SetNdivisions(1020);

   abs_KE->SetTitle( "Absorption;true Kinetic Energy (MeV); #sigma (mb)");
   abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   abs_KE->SetLineColor(kRed);
   abs_KE->SetLineWidth(3);
   abs_KE->Draw("AC");
   h_xs_trueE_trueAbs->SetMarkerSize(0.7);
   h_xs_trueE_trueAbs->Draw("PE0 SAME");

   c_trueE_abs->Write();


   TCanvas *c_trueE_totInel = new TCanvas("c_trueE_totInel", "c_trueE_totInel");
   gPad->SetGrid(1,1);
   h_xs_trueE_truePion_totInel->SetTitle( "True Total Inelastic; true Kinetic Energy (MeV); #sigma (mb)");
   h_xs_trueE_truePion_totInel->GetXaxis()->SetRangeUser(100,1100);
   h_xs_trueE_truePion_totInel->GetXaxis()->SetNdivisions(1020);
   h_xs_trueE_truePion_totInel->GetYaxis()->SetNdivisions(1020);

   totInel_KE->SetTitle("Total Inelastic MC; true kinetic Energy (MeV); #sigma (mbarn)");
   totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   totInel_KE->SetLineColor(kRed);
   totInel_KE->SetLineWidth(3);
   totInel_KE->Draw("AC");
   h_xs_trueE_truePion_totInel->SetMarkerSize(0.7);
   h_xs_trueE_truePion_totInel->Draw("PE0 SAME");

   c_trueE_totInel->Write();
   //output->Write();
   //f1.Close();
   TCanvas *c_all = new TCanvas("c_all", "c_all");
   gPad->SetGrid(1,1);

   //cex_KE->SetLineColor(kGreen);
   abs_KE->SetLineColor(kBlue);
   mg->GetXaxis()->SetRangeUser(0,1000);
   mg->Draw("AC");
   h_xs_trueE_truePion_totInel->Draw("PE0 SAME");
   h_xs_trueE_trueAbs->Draw("PE0 SAME");

   c_all->Write();
   return 0;
}

//=====================================================
//           Debug
//=====================================================


/*std::cout << "Nbins Histogram = " << nBin_int << "  Energy from 0, 1200 MeV " << std::endl;
  std::cout << "Bin             Birth     Death     Cumulative Shifted      Ratio " << std::endl;


  int birth = 0;
  int death = 0;
  int cumulative_shifted = 0;
  double ratio = 0;
  int cum_before = 0;
  for(int i = nBin_int; i >= 1; i--){

  cumulative_shifted = birth + cum_before - death;
  cum_before = cumulative_shifted;
  birth = h_trueE_truePion_inc_initE->GetBinContent(i);
  death = h_trueE_truePion_inc_interE->GetBinContent(i);
  ratio = (double) death / cumulative_shifted;

  std::cout << i << "  " << birth << "     " << death << "   " << cumulative_shifted << "    " << ratio << std::endl;

  }*/


