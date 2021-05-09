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

//compare how Pandora Shower Tag and CNN Track Id can discriminate ChEx signal from the ChExAbs signal, count the matches to the real ChEx signal from MC truth


//***********************
//Main Function

int draw_CNN_michelPlots(const string path1, const string path2){


   ROOT::RDataFrame frame(pionTree, path1);
   ROOT::RDataFrame frame_data(pionTree, path2);

   TFile *output = new TFile ("plot_CNNmichel.root", "RECREATE");

   Int_t pal_CNN[9] = { kViolet -6, kGreen + 3, kCyan +1, kCyan +1, kYellow -2, kOrange + 7, kOrange + 7, kRed -3, kRed -3};
   gStyle->SetPalette(9,pal_CNN);

   auto data_frame = frame_data.Filter("primary_isBeamType && passBeamQuality && passBeamCut");
   
   auto h_data_michelScore = data_frame.Histo1D({"h_data_michelScore","data", 20,0,1}, "reco_daughter_PFP_michelScore_collection");

   auto define_frame = frame
      .Filter("primary_isBeamType && passBeamCut && passBeamCutBI")

      .Define("gamma", pdg_gamma) //put in also nuclear gammas
      .Define("proton", pdg_proton)
      .Define("piPlus", pdg_piPlus)
      .Define("piMinus", pdg_piMinus)
      .Define("electron", pdg_electron)
      .Define("positron", pdg_positron)
      .Define("muMinus", pdg_muMinus)
      .Define("muPlus", pdg_muPlus)
      .Define("nucleus", pdg_nucleus)

      .Define("michelCNN_proton", daughter_property<std::vector<double>>, {"proton", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_michelScore_collection"})
      .Define("michelCNN_gamma", daughter_property<std::vector<double>>, {"gamma", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_michelScore_collection"})
      .Define("michelCNN_piPlus", daughter_property<std::vector<double>>, {"piPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_michelScore_collection"})
      .Define("michelCNN_piMinus", daughter_property<std::vector<double>>, {"piMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_michelScore_collection"})
      .Define("michelCNN_electron", daughter_property<std::vector<double>>, {"electron", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_michelScore_collection"})
      .Define("michelCNN_positron", daughter_property<std::vector<double>>, {"positron", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_michelScore_collection"})
      .Define("michelCNN_muMinus", daughter_property<std::vector<double>>, {"muMinus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_michelScore_collection"})
      .Define("michelCNN_muPlus", daughter_property<std::vector<double>>, {"muPlus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_michelScore_collection"})
      .Define("michelCNN_nucleus", daughter_property<std::vector<double>>, {"nucleus", "reco_daughter_PFP_true_byHits_PDG", "reco_daughter_PFP_michelScore_collection"});

   THStack *stack_cnnMichel = new THStack("cnnMichel_stack", "CNN MichelScore");


   auto h_allmc_michelScore = define_frame.Histo1D({"h_allmc_michelScore","allmc", 20,0,1}, "reco_daughter_PFP_michelScore_collection");
   double tot_data = h_data_michelScore->GetEntries();
   double tot_mc = h_allmc_michelScore->GetEntries();
   
   auto h_proton_michelScore = define_frame.Histo1D({"h_proton_michelScore","proton", 20,0,1}, "michelCNN_proton");
   auto h_gamma_michelScore = define_frame.Histo1D({"h_gamma_michelScore","gamma", 20,0,1}, "michelCNN_gamma");
   auto h_piPlus_michelScore = define_frame.Histo1D({"h_piPlus_michelScore","piPlus", 20,0,1}, "michelCNN_piPlus");
   auto h_piMinus_michelScore = define_frame.Histo1D({"h_piMinus_michelScore","piMinus", 20,0,1}, "michelCNN_piMinus");
   auto h_electron_michelScore = define_frame.Histo1D({"h_electron_michelScore","electron", 20,0,1}, "michelCNN_electron");
   auto h_positron_michelScore = define_frame.Histo1D({"h_positron_michelScore","positron", 20,0,1}, "michelCNN_positron");
   auto h_muMinus_michelScore = define_frame.Histo1D({"h_muMinus_michelScore","muMinus", 20,0,1}, "michelCNN_muMinus");
   auto h_muPlus_michelScore = define_frame.Histo1D({"h_muPlus_michelScore","muPlus", 20,0,1}, "michelCNN_muPlus");
   auto h_nucleus_michelScore = define_frame.Histo1D({"h_nucleus_michelScore","nucleus", 20,0,1}, "michelCNN_nucleus");

   h_proton_michelScore->Scale( tot_data / tot_mc);            h_proton_michelScore->Sumw2(0);
   h_gamma_michelScore->Scale( tot_data / tot_mc);             h_gamma_michelScore->Sumw2(0);
   h_piPlus_michelScore->Scale( tot_data / tot_mc);            h_piPlus_michelScore->Sumw2(0);
   h_piMinus_michelScore->Scale( tot_data / tot_mc);           h_piMinus_michelScore->Sumw2(0);
   h_electron_michelScore->Scale( tot_data / tot_mc);          h_electron_michelScore->Sumw2(0);
   h_positron_michelScore->Scale( tot_data / tot_mc);          h_positron_michelScore->Sumw2(0);
   h_muMinus_michelScore->Scale( tot_data / tot_mc);           h_muMinus_michelScore->Sumw2(0);
   h_muPlus_michelScore->Scale( tot_data / tot_mc);            h_muPlus_michelScore->Sumw2(0);
   h_nucleus_michelScore->Scale( tot_data / tot_mc);           h_nucleus_michelScore->Sumw2(0);

   h_proton_michelScore->SetFillColor( kViolet -6 );
   h_gamma_michelScore->SetFillColor( kGreen + 3); 
   h_piPlus_michelScore->SetFillColor( kCyan +1);
   h_piMinus_michelScore->SetFillColor( kCyan +1); 
   h_electron_michelScore->SetFillColor( kOrange + 7);
   h_positron_michelScore->SetFillColor( kOrange + 7);
   h_muMinus_michelScore->SetFillColor( kRed -3 ); 
   h_muPlus_michelScore->SetFillColor( kRed -3 ); 
   h_nucleus_michelScore->SetFillColor(kYellow -2 );    
   

   stack_cnnMichel->Add(h_proton_michelScore.GetPtr());
   stack_cnnMichel->Add(h_gamma_michelScore.GetPtr());
   stack_cnnMichel->Add(h_piPlus_michelScore.GetPtr());
   stack_cnnMichel->Add(h_piMinus_michelScore.GetPtr());
   stack_cnnMichel->Add(h_nucleus_michelScore.GetPtr());
   stack_cnnMichel->Add(h_muMinus_michelScore.GetPtr());
   stack_cnnMichel->Add(h_muPlus_michelScore.GetPtr()); 
   stack_cnnMichel->Add(h_electron_michelScore.GetPtr());
   stack_cnnMichel->Add(h_positron_michelScore.GetPtr());

   auto c_CNNmichel = new TCanvas("c_CNNmichel", "", 1600,2000);
   stack_cnnMichel->Draw("PFC");
   h_data_michelScore->Draw("PE SAME");
   c_CNNmichel->BuildLegend();
   c_CNNmichel->Write();








   output->Write();
   return 0;
}

