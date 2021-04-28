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
#include "TPad.h"
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

int compare_incInt_dataMC(){

   //
   gStyle->SetNdivisions(1020);

   //output file
   TFile *output = new TFile ("compare_incInt_dataMC_effPur_francesco.root", "RECREATE");

   //access Jakes GeantFile in folder
   TFile f1("output_eSliceMethod_selectedEvents_MC_effPur_francesco.root");
   TH1D *mc_incident_pion = (TH1D*)f1.Get("h_selected_pion_incidentRecoE");
   TH1D *mc_interacting_abs = (TH1D*)f1.Get("h_selected_abs_interactingRecoE");
   //TH1D *mc_interacting_totInel = (TH1D*)f1.Get("h_selected_totInel_interactingRecoE");
   //TH1D *mc_interacting_cex = (TH1D*)f1.Get("h_selected_cex_interactingRecoE");


   TFile f2("output_eSliceMethod_selectedEvents_DATA_effPur_francesco.root");
   TH1D *data_incident_pion = (TH1D*)f2.Get("h_selected_pion_incidentRecoE");
   TH1D *data_interacting_abs = (TH1D*)f2.Get("h_selected_abs_interactingRecoE");
   //TH1D *data_interacting_totInel = (TH1D*)f2.Get("h_selected_totInel_interactingRecoE");
   //TH1D *data_interacting_cex = (TH1D*)f2.Get("h_selected_cex_interactingRecoE");

   output->cd();

   //Scale mc to number of Entries of data

   mc_incident_pion->Scale( data_incident_pion->Integral() / mc_incident_pion->Integral() );
   mc_interacting_abs->Scale( data_interacting_abs->Integral() / mc_interacting_abs->Integral() );
   //mc_interacting_cex->Scale( data_interacting_cex->Integral() / mc_interacting_cex->Integral() );
   //mc_interacting_totInel->Scale( data_interacting_totInel->Integral() / mc_interacting_totInel->Integral() );

   mc_incident_pion->Sumw2(0);
   mc_interacting_abs->Sumw2(0);
   //mc_interacting_cex->Sumw2(0);
   //mc_interacting_totInel->Sumw2(0);
   
   mc_incident_pion->SetNameTitle("mc_inc_pi", "MC incident Pion");
   mc_interacting_abs->SetNameTitle("mc_int_abs", "MC interacting Abs");
   //mc_interacting_cex->SetNameTitle("mc_int_cex", "MC interacting Cex");
   //mc_interacting_totInel->SetNameTitle("mc_int_totInel","MC interacting total Inelastic" );
   
   data_incident_pion->SetNameTitle("data_inc_pi", "data incident Pion");
   data_interacting_abs->SetNameTitle("data_int_abs", "data interacting Abs");
   //data_interacting_cex->SetNameTitle("data_int_cex", "data interacting Cex");
   //data_interacting_totInel->SetNameTitle("data_int_totInel","data interacting total Inelastic" );

   mc_incident_pion->SetLineColor(kBlue);
   mc_interacting_abs->SetLineColor(kBlue);
   //mc_interacting_cex->SetLineColor(kBlue);
   //mc_interacting_totInel->SetLineColor(kBlue);
   
   data_incident_pion->SetLineColor(kGreen + 2);
   data_interacting_abs->SetLineColor(kGreen + 2);
   //data_interacting_cex->SetLineColor(kGreen + 2);
   //data_interacting_totInel->SetLineColor(kGreen + 2);
 
   mc_incident_pion->Write();
   mc_interacting_abs->Write();
   //mc_interacting_cex->Write();
   //mc_interacting_totInel->Write();
   
   data_incident_pion->Write();
   data_interacting_abs->Write();
   //data_interacting_cex->Write();
   //data_interacting_totInel->Write();

   TH1D* div_dat_mc_inc = (TH1D*) data_incident_pion->Clone("div_dat_mc_inc");
   div_dat_mc_inc->SetTitle("Divide Data over MC Incident");
   div_dat_mc_inc->Divide( mc_incident_pion );
   div_dat_mc_inc->Write(); 
 
   TH1D* div_dat_mc_int_abs = (TH1D*) data_interacting_abs->Clone("div_dat_mc_int_abs");
   div_dat_mc_int_abs->SetTitle("Divide Data over MC Interacting Abs");
   div_dat_mc_int_abs->Divide( mc_interacting_abs );
   div_dat_mc_int_abs->Write(); 
   

   THStack *stack_inc = new THStack("stack_inc", "stacked incident");
   stack_inc->Add(mc_incident_pion);
   stack_inc->Add(data_incident_pion);

   stack_inc->Write();

   THStack *stack_int_abs = new THStack("stack_int_abs", "stacked interacting Abs");
   stack_int_abs->Add(mc_interacting_abs);
   stack_int_abs->Add(data_interacting_abs);

   stack_int_abs->Write(); 
   
   //THStack *stack_int_cex = new THStack("stack_int_cex", "stacked interacting Cex");
   //stack_int_cex->Add(mc_interacting_cex);
   //stack_int_cex->Add(data_interacting_cex);

   //stack_int_cex->Write(); 
   
   //THStack *stack_int_totInel = new THStack("stack_int_totInel", "stacked interacting totInel");
   //stack_int_totInel->Add(mc_interacting_totInel);
   //stack_int_totInel->Add(data_interacting_totInel);

   //stack_int_totInel->Write(); 

   TCanvas *c_same_inc = new TCanvas("c_same_inc", "");
   gPad->SetGrid(1,1);
   mc_incident_pion->Draw("HIST");
   data_incident_pion->Draw("HIST SAME");

   TCanvas *c_same_int_abs = new TCanvas("c_same_int_abs", "");
   gPad->SetGrid(1,1);
   mc_interacting_abs->Draw("HIST");
   data_interacting_abs->Draw("HIST SAME");
   
   //TCanvas *c_same_int_cex = new TCanvas("c_same_int_cex", "");
   //gPad->SetGrid(1,1);
   //mc_interacting_cex->Draw("HIST");
   //data_interacting_cex->Draw("HIST SAME");
   //
   //TCanvas *c_same_int_totInel = new TCanvas("c_same_int_totInel", "");
   //gPad->SetGrid(1,1);
   //mc_interacting_totInel->Draw("HIST");
   //data_interacting_totInel->Draw("HIST SAME");


   c_same_inc->Write();
   c_same_int_abs->Write();
   //c_same_int_cex->Write();
   //c_same_int_totInel->Write();
   output->Write();
   //f1.Close();
   return 0;
}

