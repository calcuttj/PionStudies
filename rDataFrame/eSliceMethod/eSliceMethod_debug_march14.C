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

int eSliceMethod_debug_march14(){

   //gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   //ROOT::RDataFrame frame(inputTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   
   //output file
   TFile *output = new TFile ("debug.root", "RECREATE");

   //access Jakes GeantFile in folder
   TFile f1("output_eSliceMethod_trueProcess_trueEnergy.root");
   TH1D *true_pion_inc = (TH1D*)f1.Get("h_true_pion_true_incidentE");
   TH1D *true_abs_int = (TH1D*)f1.Get("h_true_abs_true_interactingE");
   TH1D *true_xs_abs = (TH1D*)f1.Get("h_xs_true_abs_true_E");

   if(true_pion_inc) std::cout << "Accessed first file" << std::endl;

   TFile f2("output_eSliceMethod_selectedEvents_onlyTrueInt.root");
   TH1D *selected_pion_inc = (TH1D*)f2.Get("h_selected_pion_incidentE");
   TH1D *selected_abs_int = (TH1D*)f2.Get("h_selected_abs_interactingE");
   TH1D *selected_xs_abs = (TH1D*)f2.Get("h_xs_selected_abs");

   if(selected_pion_inc) std::cout << "Accessed second file" << std::endl;
   
   output->cd();

   //Scale True to number of Entries of selected
   std::cout << "Entries true incident = " << true_pion_inc->GetEntries() << std::endl;
   std::cout << "Entries selected incident = " << selected_pion_inc->GetEntries() << std::endl;

   true_pion_inc->Scale( selected_pion_inc->Integral() / true_pion_inc->Integral() );
   true_abs_int->Scale( selected_abs_int->Integral() / true_abs_int->Integral() );

   //Normalise Incident and Interacting Histos to own integral
   //true_pion_inc->Scale( 1 / true_pion_inc->Integral() );
   //true_abs_int->Scale( 1 / true_abs_int->Integral() );
   //selected_pion_inc->Scale( 1 / selected_pion_inc->Integral() );
   //selected_abs_int->Scale( 1 / selected_abs_int->Integral() );

   true_pion_inc->Write();
   true_abs_int->Write();
   selected_pion_inc->Write();
   selected_abs_int->Write();

   //Names
   true_pion_inc->SetNameTitle("true_inc_pi", "true incident Pion");
   true_pion_inc->SetFillColor(kBlue);

   true_abs_int->SetNameTitle("true_abs_int", "true abs interacting");
   true_abs_int->SetFillColor(kBlue);

   selected_pion_inc->SetNameTitle("selected_true_pi", "selected true incident Pion");
   selected_pion_inc->SetFillColor(kRed);

   selected_abs_int->SetNameTitle("selected_true_abs", "selected true abs interacting");
   selected_abs_int->SetFillColor(kRed);

   THStack *stack_inc = new THStack("stack_inc", "stacked incident");
   stack_inc->Add(true_pion_inc);
   stack_inc->Add(selected_pion_inc);

   stack_inc->Write();

   THStack *stack_int = new THStack("stack_int", "stacked interacting");
   stack_int->Add(true_abs_int);
   stack_int->Add(selected_abs_int);

   stack_int->Write(); 
   
   TCanvas *c_stack_inc = new TCanvas("c_stack_inc", "");
   gPad->SetGrid(1,1);
   stack_inc->Draw("nostackb HIST");
      
   TCanvas *c_stack_int = new TCanvas("c_stack_int", "");
   gPad->SetGrid(1,1);
   stack_int->Draw("nostackb HIST");
   //c_stack->Close();

   
   c_stack_inc->Write();
   c_stack_int->Write();
   output->Write();
   //f1.Close();
   return 0;
}

