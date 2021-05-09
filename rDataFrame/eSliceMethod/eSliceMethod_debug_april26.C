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

int eSliceMethod_debug_april26(){

   //gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   //ROOT::RDataFrame frame(inputTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   
   //output file
   TFile *output = new TFile ("debug_april26.root", "RECREATE");

   //access Jakes GeantFile in folder
   TFile f1("output_eSliceMethod_selectedEvents_MC_effPur_noSmear.root");
   TH1D *mc_pion_inc = (TH1D*)f1.Get("h_selected_pion_incidentRecoE");
   TH1D *mc_abs_int = (TH1D*)f1.Get("h_selected_abs_interactingRecoE");
   //TH1D *true_xs_abs = (TH1D*)f1.Get("h_xs_true_abs_true_E");

   if(mc_pion_inc) std::cout << "Accessed first file" << std::endl;

   TFile f2("output_eSliceMethod_selectedEvents_DATA_effPur_noSmear.root");
   TH1D *data_pion_inc = (TH1D*)f2.Get("h_selected_pion_incidentRecoE");
   TH1D *data_abs_int = (TH1D*)f2.Get("h_selected_abs_interactingRecoE");
   //TH1D *selected_xs_abs = (TH1D*)f2.Get("h_xs_selected_abs");

   if(data_pion_inc) std::cout << "Accessed second file" << std::endl;
   
   output->cd();

   //Scale True to number of Entries of selected
   //std::cout << "Entries true incident = " << mc_pion_inc->GetEntries() << std::endl;
   //std::cout << "Entries selected incident = " << data_pion_inc->GetEntries() << std::endl;

   mc_pion_inc->Scale( data_pion_inc->Integral() / mc_pion_inc->Integral() );
   mc_pion_inc->Sumw2(0);
   mc_abs_int->Scale( data_abs_int->Integral() / mc_abs_int->Integral() );
   mc_abs_int->Sumw2(0);

   //Normalise Incident and Interacting Histos to own integral
   //mc_pion_inc->Scale( 1 / mc_pion_inc->Integral() );
   //mc_abs_int->Scale( 1 / mc_abs_int->Integral() );
   //data_pion_inc->Scale( 1 / data_pion_inc->Integral() );
   //data_abs_int->Scale( 1 / data_abs_int->Integral() );

   mc_pion_inc->Write();
   mc_abs_int->Write();
   data_pion_inc->Write();
   data_abs_int->Write();

   //=====================================================
   TH1D* h_betheMean_muon = new TH1D("h_betheMean_muon", "Mean Energy Loss", nBin_int, eEnd, eStart);

   //fill histo with Mean dEdX of bin center
   for(int i = 1; i <= nBin_int; i++){
      h_betheMean_muon->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size_int  , mass_muon) );
   };

   double scale_factor = factor_mbarn * atomic_mass / ( density * N_avogadro * bin_size_int );

   //compute XS by taking DATA Int and MC Inc
   
   TH1D* h_xs_dataInt_mcInc = (TH1D*) mc_pion_inc->Clone("h_xs_dataInt_mcInc");
   
   TH1D* dummy_recoE_abs = (TH1D*) mc_pion_inc->Clone("");
   dummy_recoE_abs->Add( data_abs_int, -1);
   
   h_xs_dataInt_mcInc->Divide( dummy_recoE_abs );
   
   for(int i = 1; i <= nBin_int; i++) h_xs_dataInt_mcInc->SetBinContent(i, log( h_xs_dataInt_mcInc->GetBinContent(i) ));
   h_xs_dataInt_mcInc->Multiply( h_betheMean_muon );
   h_xs_dataInt_mcInc->Scale( scale_factor );
   h_xs_dataInt_mcInc->Sumw2(0);

   h_xs_dataInt_mcInc->Write();
///
//
//
   TH1D* h_xs_mcInt_dataInc = (TH1D*) data_pion_inc->Clone("h_xs_mcInt_dataInc");
   
   TH1D* stupid_recoE_abs = (TH1D*) data_pion_inc->Clone("");
   stupid_recoE_abs->Add( mc_abs_int, -1);
   
   h_xs_mcInt_dataInc->Divide( stupid_recoE_abs );
   
   for(int i = 1; i <= nBin_int; i++) h_xs_mcInt_dataInc->SetBinContent(i, log( h_xs_mcInt_dataInc->GetBinContent(i) ));
   h_xs_mcInt_dataInc->Multiply( h_betheMean_muon );
   h_xs_mcInt_dataInc->Scale( scale_factor );
   h_xs_mcInt_dataInc->Sumw2(0);

   h_xs_mcInt_dataInc->Write();

   return 0;
}

