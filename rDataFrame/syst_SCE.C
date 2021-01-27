#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "lambda.h"
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>

//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;


//***********************
//Main Function
//Get estimate on Energy of primary pions for MC with and without SCE

int syst_SCE(const string path1, const string path2, const string path3){

   //Int_t palette_sce[] = {kRed, kBlue, kRed-6, kBlue-10};
   //gStyle->SetPalette(4,palette_sce);

   //read in Data, MC SCE and MC without SCE
   ROOT::RDataFrame frame_mc_SCE(inputTree, path1);
   ROOT::RDataFrame frame_mc_noSCE(inputTree, path2);
   ROOT::RDataFrame frame_data(inputTree, path3);

   TFile *output = new TFile ("output_study.root", "RECREATE");

//   THStack *stack_incidentEnergy = new THStack ("stack_incidentEnergy", "incident Energy");
//   THStack *stack_interactingEnergy = new THStack ("stack_interactingEnergy", "interacting Energy");


   auto h_mcSCE_reco_incidentEnergy = frame_mc_SCE.Histo1D({"mc_sce_incidentE", "mc SCE", 150, 0, 1500}, "reco_beam_incidentEnergies");
   auto h_mcSCE_reco_interactingEnergy = frame_mc_SCE.Histo1D({"mc_sce_interactingE", "mc SCE", 150, 0, 1500}, "reco_beam_interactingEnergy");
   auto h_mcSCE_true_incidentEnergy = frame_mc_SCE.Histo1D({"mc_sce_incidentE", "true mc SCE", 150, 0, 1500}, "true_beam_incidentEnergies");
   auto h_mcSCE_true_interactingEnergy = frame_mc_SCE.Histo1D({"mc_sce_interactingE", "true mc SCE", 150, 0, 1500}, "true_beam_interactingEnergy");
   
   auto h_mcnoSCE_reco_incidentEnergy = frame_mc_noSCE.Histo1D({"mc_noSCE_incidentE", "mc noSCE", 150, 0, 1500}, "reco_beam_incidentEnergies");
   auto h_mcnoSCE_reco_interactingEnergy = frame_mc_noSCE.Histo1D({"mc_noSCE_interactingE", "mc noSCE", 150, 0, 1500}, "reco_beam_interactingEnergy");
   auto h_mcnoSCE_true_incidentEnergy = frame_mc_noSCE.Histo1D({"mc_noSCE_incidentE", "true mc noSCE", 150, 0, 1500}, "true_beam_incidentEnergies");
   auto h_mcnoSCE_true_interactingEnergy = frame_mc_noSCE.Histo1D({"mc_noSCE_interactingE", "true mc noSCE", 150, 0, 1500}, "true_beam_interactingEnergy");
   
   auto h_data_reco_incidentEnergy = frame_data.Histo1D({"data_incidentE", "data", 150, 0, 1500}, "reco_beam_incidentEnergies");
   auto h_data_reco_interactingEnergy = frame_data.Histo1D({"data_interactingE", "data", 150, 0, 1500}, "reco_beam_interactingEnergy");

   //Scaling Histograms to data

/*   double int_data_incidentEnergy = h_data_reco_incidentEnergy->Integral();
   double int_data_interactingEnergy = h_data_reco_interactingEnergy->Integral();

   double sum_incidentE = h_mcSCE_reco_incidentEnergy->Integral() + h_mcSCE_true_incidentEnergy->Integral()
           + h_mcnoSCE_reco_incidentEnergy->Integral() + h_mcnoSCE_true_incidentEnergy->Integral();
   
   double sum_interactingE = h_mcSCE_reco_interactingEnergy->Integral() + h_mcSCE_true_interactingEnergy->Integral()
           + h_mcnoSCE_reco_interactingEnergy->Integral() + h_mcnoSCE_true_interactingEnergy->Integral();


   h_mcSCE_reco_incidentEnergy->Scale( int_data_incidentEnergy / sum_incidentE );
   h_mcSCE_reco_interactingEnergy->Scale( int_data_interactingEnergy / sum_interactingE );
   
   h_mcSCE_true_incidentEnergy->Scale( int_data_incidentEnergy / sum_incidentE );
   h_mcSCE_true_interactingEnergy->Scale( int_data_interactingEnergy / sum_interactingE );
   
   h_mcnoSCE_reco_incidentEnergy->Scale( int_data_incidentEnergy / sum_incidentE );
   h_mcnoSCE_reco_interactingEnergy->Scale( int_data_interactingEnergy / sum_interactingE );
   
   h_mcnoSCE_true_incidentEnergy->Scale( int_data_incidentEnergy / sum_incidentE );
   h_mcnoSCE_true_interactingEnergy->Scale( int_data_interactingEnergy / sum_interactingE );
*/
   double int_data_incidentEnergy = h_data_reco_incidentEnergy->Integral();
   double int_data_interactingEnergy = h_data_reco_interactingEnergy->Integral();

   h_mcSCE_reco_incidentEnergy->Scale( int_data_incidentEnergy / h_mcSCE_reco_incidentEnergy->Integral());
   h_mcSCE_reco_interactingEnergy->Scale( int_data_interactingEnergy / h_mcSCE_reco_interactingEnergy->Integral() );
   
   h_mcSCE_true_incidentEnergy->Scale( int_data_incidentEnergy /h_mcSCE_true_incidentEnergy->Integral());
   h_mcSCE_true_interactingEnergy->Scale( int_data_interactingEnergy / h_mcSCE_true_interactingEnergy->Integral() );
   
   h_mcnoSCE_reco_incidentEnergy->Scale( int_data_incidentEnergy / h_mcnoSCE_reco_incidentEnergy->Integral());
   h_mcnoSCE_reco_interactingEnergy->Scale( int_data_interactingEnergy / h_mcnoSCE_reco_interactingEnergy->Integral() );
   
   h_mcnoSCE_true_incidentEnergy->Scale( int_data_incidentEnergy / h_mcnoSCE_true_incidentEnergy->Integral() );
   h_mcnoSCE_true_interactingEnergy->Scale( int_data_interactingEnergy /h_mcnoSCE_true_interactingEnergy->Integral() );


   h_mcSCE_reco_incidentEnergy->Sumw2(0);
   h_mcSCE_reco_interactingEnergy->Sumw2(0);                                    
   h_mcSCE_true_incidentEnergy->Sumw2(0);  
   h_mcSCE_true_interactingEnergy->Sumw2(0);                                
   h_mcnoSCE_reco_incidentEnergy->Sumw2(0);
   h_mcnoSCE_reco_interactingEnergy->Sumw2(0); 
   h_mcnoSCE_true_incidentEnergy->Sumw2(0);
   h_mcnoSCE_true_interactingEnergy->Sumw2(0);

   h_mcSCE_reco_incidentEnergy->SetLineColor(kRed);
   h_mcSCE_reco_interactingEnergy->SetLineColor(kRed);                                    
   h_mcSCE_true_incidentEnergy->SetLineColor(kRed-6);  
   h_mcSCE_true_interactingEnergy->SetLineColor(kRed-6);                                
   h_mcnoSCE_reco_incidentEnergy->SetLineColor(kBlue);
   h_mcnoSCE_reco_interactingEnergy->SetLineColor(kBlue); 
   h_mcnoSCE_true_incidentEnergy->SetLineColor(kBlue-10);
   h_mcnoSCE_true_interactingEnergy->SetLineColor(kBlue-10);


   /*
   stack_incidentEnergy->Add(h_mcSCE_reco_incidentEnergy.GetPtr());
   stack_incidentEnergy->Add(h_mcnoSCE_reco_incidentEnergy.GetPtr());
   stack_incidentEnergy->Add(h_mcSCE_true_incidentEnergy.GetPtr());
   stack_incidentEnergy->Add(h_mcnoSCE_true_incidentEnergy.GetPtr());

   stack_interactingEnergy->Add(h_mcSCE_reco_interactingEnergy.GetPtr());
   stack_interactingEnergy->Add(h_mcnoSCE_reco_interactingEnergy.GetPtr());
   stack_interactingEnergy->Add(h_mcSCE_true_interactingEnergy.GetPtr());
   stack_interactingEnergy->Add(h_mcnoSCE_true_interactingEnergy.GetPtr());
*/
   auto c_incidentE = new TCanvas("incidentEnergy", "Incident Energy", 1600,1200);

   //stack_incidentEnergy->Draw("PFC HIST");
   h_mcSCE_reco_incidentEnergy->Draw("HIST");
   h_mcSCE_true_incidentEnergy->Draw("HIST SAME");
   h_mcnoSCE_reco_incidentEnergy->Draw("HIST SAME");
   h_mcnoSCE_true_incidentEnergy->Draw("HIST SAME");
   h_data_reco_incidentEnergy->Draw("PSAME SAME");
   c_incidentE->BuildLegend();
   c_incidentE->Write();
   c_incidentE->Close();

   auto c_interactingE = new TCanvas("interactingEnergy", "interacting Energy", 1600,1200);

   //stack_interactingEnergy->Draw("PFC HIST");
   h_mcSCE_reco_interactingEnergy->Draw("HIST");
   h_mcSCE_true_interactingEnergy->Draw("HIST SAME");
   h_mcnoSCE_reco_interactingEnergy->Draw("HIST SAME");
   h_mcnoSCE_true_interactingEnergy->Draw("HIST SAME");
   h_data_reco_interactingEnergy->Draw("PSAME");
   c_interactingE->BuildLegend();
   c_interactingE->Write();
   c_interactingE->Close();

   h_mcSCE_reco_incidentEnergy->Write();
   h_mcSCE_reco_interactingEnergy->Write();                                    
   h_mcSCE_true_incidentEnergy->Write();  
   h_mcSCE_true_interactingEnergy->Write();                                
   h_mcnoSCE_reco_incidentEnergy->Write();
   h_mcnoSCE_reco_interactingEnergy->Write(); 
   h_mcnoSCE_true_incidentEnergy->Write();
   h_mcnoSCE_true_interactingEnergy->Write();
   h_data_reco_incidentEnergy->Write();
   h_data_reco_interactingEnergy->Write();


   output->Write();
   return 0;
}
