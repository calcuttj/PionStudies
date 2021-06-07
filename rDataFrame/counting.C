#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "lambda.h"
#include "eventSelection.h"
#include <ROOT/RDataFrame.hxx>

#include "TGraphAsymmErrors.h"

#include "backgrounds.h"
#include "selection_defs.h"

#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>

const std::string mc_all_file = "eventSelection_mc_all.root";
const std::string data_all_file = "eventSelection_data_all.root";
const string mc_beamType_file = "eventSelection_mc_beamType.root";
const string mc_endAPA3_file = "eventSelection_mc_PRIMARYPION.root";
const std::string mc_noPionDaughter_file = "eventSelection_mc_COMBINED.root";
const std::string data_noPionDaughter_file = "eventSelection_data_COMBINED.root";
const string mc_signal_cex_file = "eventSelection_mc_CEX.root";
const string data_signal_cex_file = "eventSelection_data_CEX.root";
const std::string mc_signal_abs_file = "eventSelection_mc_ABS.root";
const std::string data_signal_abs_file = "eventSelection_data_ABS.root";
             

int counting(bool doCounting = false) {

  gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
  gROOT->SetBatch(1);

  TFile *output = new TFile ("output_eventSelection.root", "RECREATE");
  THStack *stack_cutFlow = new THStack("cutFlow", "Cut Flow MC and Data");

  //************** CUT Flow Histogram ************//
  TH1I *h_mc_total = new TH1I("mc_total", "MC selected", 7 , 0., 7.);
  TH1I *h_data_total = new TH1I("data_total", "Data selected", 7 , 0., 7.);
  TH1I *h_true_abs = new TH1I("true_abs", "True Abs Signal", 7 , 0., 7.);
  TH1I *h_true_cex = new TH1I("true_cex", "True Cex Signal", 7 , 0., 7.);
  TH1I *h_true_nPi0 = new TH1I("true_nPi0", "True N-Pi0 Signal", 7 , 0., 7.);
  TH1I *h_true_BG = new TH1I("true_BG", "Background", 7 , 0., 7.);


  //Before filling, split it up MC by the true signal/BG
  //This will later be used in the counting as well
  
  ROOT::RDataFrame mc_all(inputTree, mc_all_file);
  ROOT::RDataFrame mcCUT_beamType(inputTree, mc_beamType_file);
  auto mcCUT_beamCut = mcCUT_beamType.Filter("passBeamCut");
  ROOT::RDataFrame mcCUT_endAPA3(inputTree, mc_endAPA3_file);
  ROOT::RDataFrame mcCUT_noPionDaughter(inputTree, mc_noPionDaughter_file);
  ROOT::RDataFrame mcSIGNAL_cex(inputTree, mc_signal_cex_file);
  ROOT::RDataFrame mcSIGNAL_abs(inputTree, mc_signal_abs_file);

  ROOT::RDataFrame data_all(inputTree, data_all_file);
  auto dataCUT_beamType = data_all.Filter("primary_isBeamType");
  auto dataCUT_beamCut = dataCUT_beamType.Filter("passBeamQuality && passBeamCut");
  auto dataCUT_endAPA3 = dataCUT_beamCut.Filter("primary_ends_inAPA3"); 
  ROOT::RDataFrame dataCUT_noPionDaughter(inputTree, data_noPionDaughter_file);
  ROOT::RDataFrame dataSIGNAL_cex(inputTree, data_signal_cex_file);
  ROOT::RDataFrame dataSIGNAL_abs(inputTree, data_signal_abs_file);
  
  auto beamType_absSignal = mcCUT_beamType.Filter("true_absSignal").Count();
  auto beamType_chexSignal = mcCUT_beamType.Filter("true_chexSignal").Count();
  auto beamType_nPi0Signal = mcCUT_beamType.Filter("true_nPi0Signal").Count();
  auto beamType_background = mcCUT_beamType.Filter("true_backGround").Count();

  auto beamCut_absSignal = mcCUT_beamCut.Filter("true_absSignal").Count();
  auto beamCut_chexSignal = mcCUT_beamCut.Filter("true_chexSignal").Count();
  auto beamCut_nPi0Signal = mcCUT_beamCut.Filter("true_nPi0Signal").Count();
  auto beamCut_background = mcCUT_beamCut.Filter("true_backGround").Count();

  auto endAPA3_absSignal_filter = mcCUT_endAPA3.Filter("true_absSignal")/*.Count()*/;
  auto endAPA3_chexSignal_filter = mcCUT_endAPA3.Filter("true_chexSignal")/*.Count()*/;
  auto endAPA3_nPi0Signal_filter = mcCUT_endAPA3.Filter("true_nPi0Signal")/*.Count()*/;
  auto endAPA3_background_filter = mcCUT_endAPA3.Filter("true_backGround")/*.Count()*/;

  auto endAPA3_absSignal =  endAPA3_absSignal_filter.Count();
  auto endAPA3_chexSignal = endAPA3_chexSignal_filter.Count();
  auto endAPA3_nPi0Signal = endAPA3_nPi0Signal_filter.Count();
  auto endAPA3_background = endAPA3_background_filter.Count();

  auto noPionDaughter_absSignal = mcCUT_noPionDaughter.Filter("true_absSignal").Count();
  auto noPionDaughter_chexSignal = mcCUT_noPionDaughter.Filter("true_chexSignal").Count();
  auto noPionDaughter_nPi0Signal = mcCUT_noPionDaughter.Filter("true_nPi0Signal").Count();
  auto noPionDaughter_background = mcCUT_noPionDaughter.Filter("true_backGround").Count();

  auto chex_absSignal = mcSIGNAL_cex.Filter("true_absSignal").Count();
  auto chex_chexSignal = mcSIGNAL_cex.Filter("true_chexSignal").Count();
  auto chex_nPi0Signal = mcSIGNAL_cex.Filter("true_nPi0Signal").Count();
  auto chex_background = mcSIGNAL_cex.Filter("true_backGround").Count();


  auto abs_absSignal_filter = mcSIGNAL_abs.Filter("true_absSignal");
  auto abs_chexSignal_filter = mcSIGNAL_abs.Filter("true_chexSignal");
  auto abs_nPi0Signal_filter = mcSIGNAL_abs.Filter("true_nPi0Signal");
  auto abs_background_filter = mcSIGNAL_abs.Filter("true_backGround");

  auto abs_absSignal =  abs_absSignal_filter.Count();
  auto abs_chexSignal = abs_chexSignal_filter.Count();
  auto abs_nPi0Signal = abs_nPi0Signal_filter.Count();
  auto abs_background = abs_background_filter.Count();
  
  auto n_mc_all = mc_all.Count();
  auto n_data_all = data_all.Count();
  auto n_true_combinedSignal = mc_all.Filter("true_combinedSignal").Count();
  auto n_true_absSignal = mc_all.Filter("true_absSignal").Count();
  auto n_true_cexSignal = mc_all.Filter("true_chexSignal").Count();
  auto n_true_nPi0Signal = mc_all.Filter("true_nPi0Signal").Count();
  auto n_true_backGround = mc_all.Filter("true_backGround").Count();

  h_mc_total->SetBinContent(1 , *n_mc_all );
  h_data_total->SetBinContent(1, *n_data_all );
  h_true_abs->SetBinContent(1, *n_true_absSignal );
  h_true_cex->SetBinContent(1, *n_true_cexSignal );
  h_true_nPi0->SetBinContent(1, *n_true_nPi0Signal );
  h_true_BG->SetBinContent(1, *n_true_backGround );

  auto N_mcCUT_beamType = mcCUT_beamType.Count();
  auto N_dataCUT_beamType = dataCUT_beamType.Count();
  //auto beamType_absSignal = mcCUT_beamType.Filter("true_absSignal").Count();
  //auto beamType_chexSignal = mcCUT_beamType.Filter("true_chexSignal").Count();
  //auto beamType_nPi0Signal = mcCUT_beamType.Filter("true_nPi0Signal").Count();
  //auto beamType_background = mcCUT_beamType.Filter("true_backGround").Count();

  h_mc_total->SetBinContent(2 , *N_mcCUT_beamType );
  h_data_total->SetBinContent(2, *N_dataCUT_beamType);
  h_true_abs->SetBinContent(2, *beamType_absSignal);
  h_true_cex->SetBinContent(2, *beamType_chexSignal);
  h_true_nPi0->SetBinContent(2, *beamType_nPi0Signal);
  h_true_BG->SetBinContent(2, *beamType_background);

  auto N_mcCUT_beamCut = mcCUT_beamCut.Count();
  auto N_dataCUT_beamCut = dataCUT_beamCut.Count();
  //auto beamCut_absSignal = mcCUT_beamCut.Filter("true_absSignal").Count();
  //auto beamCut_chexSignal = mcCUT_beamCut.Filter("true_chexSignal").Count();
  //auto beamCut_nPi0Signal = mcCUT_beamCut.Filter("true_nPi0Signal").Count();
  //auto beamCut_background = mcCUT_beamCut.Filter("true_backGround").Count();

  h_mc_total->SetBinContent(3 , *N_mcCUT_beamCut );
  h_data_total->SetBinContent(3, *N_dataCUT_beamCut);
  h_true_abs->SetBinContent(3, *beamCut_absSignal);
  h_true_cex->SetBinContent(3, *beamCut_chexSignal);
  h_true_nPi0->SetBinContent(3, *beamCut_nPi0Signal);
  h_true_BG->SetBinContent(3, *beamCut_background);

  auto N_mcCUT_endAPA3 = mcCUT_endAPA3.Count();
  auto N_dataCUT_endAPA3 = dataCUT_endAPA3.Count();
  //auto endAPA3_absSignal =  mcCUT_endAPA3.Filter("true_absSignal");
  //auto endAPA3_chexSignal = mcCUT_endAPA3.Filter("true_chexSignal");
  //auto endAPA3_nPi0Signal = mcCUT_endAPA3.Filter("true_nPi0Signal");
  //auto endAPA3_background = mcCUT_endAPA3.Filter("true_backGround");

  h_mc_total->SetBinContent(4 , *N_mcCUT_endAPA3 );
  h_data_total->SetBinContent(4, *N_dataCUT_endAPA3);
  h_true_abs->SetBinContent(4, *endAPA3_absSignal );
  h_true_cex->SetBinContent(4, *endAPA3_chexSignal );
  h_true_nPi0->SetBinContent(4, *endAPA3_nPi0Signal );
  h_true_BG->SetBinContent(4, *endAPA3_background );

  auto N_mcCUT_noPionDaughter = mcCUT_noPionDaughter.Count();
  auto N_dataCUT_noPionDaughter = dataCUT_noPionDaughter.Count();
  /*
  auto noPionDaughter_absSignal =
      mcCUT_noPionDaughter.Filter("true_absSignal").Count();
  auto noPionDaughter_chexSignal =
      mcCUT_noPionDaughter.Filter("true_chexSignal").Count();
  auto noPionDaughter_nPi0Signal =
      mcCUT_noPionDaughter.Filter("true_nPi0Signal").Count();
  auto noPionDaughter_background =
      mcCUT_noPionDaughter.Filter("true_backGround").Count();
      */


  h_mc_total->SetBinContent(5 , *N_mcCUT_noPionDaughter );
  h_data_total->SetBinContent(5, *N_dataCUT_noPionDaughter);
  h_true_abs->SetBinContent(5, *noPionDaughter_absSignal );
  h_true_cex->SetBinContent(5, *noPionDaughter_chexSignal );
  h_true_nPi0->SetBinContent(5, *noPionDaughter_nPi0Signal );
  h_true_BG->SetBinContent(5, *noPionDaughter_background );

  auto N_mcSIGNAL_cex = mcSIGNAL_cex.Count();
  auto N_dataSIGNAL_cex = dataSIGNAL_cex.Count();
  /*
  auto chex_absSignal = mcSIGNAL_cex.Filter("true_absSignal").Count();
  auto chex_chexSignal = mcSIGNAL_cex.Filter("true_chexSignal").Count();
  auto chex_nPi0Signal = mcSIGNAL_cex.Filter("true_nPi0Signal").Count();
  auto chex_background = mcSIGNAL_cex.Filter("true_backGround").Count();
  */


  h_mc_total->SetBinContent(6 , *N_mcSIGNAL_cex );
  h_data_total->SetBinContent(6, *N_dataSIGNAL_cex);
  h_true_abs->SetBinContent(6,  *chex_absSignal );
  h_true_cex->SetBinContent(6,  *chex_chexSignal );
  h_true_nPi0->SetBinContent(6, *chex_nPi0Signal );
  h_true_BG->SetBinContent(6,   *chex_background );

  auto N_mcSIGNAL_abs = mcSIGNAL_abs.Count();
  auto N_dataSIGNAL_abs = dataSIGNAL_abs.Count();

  h_mc_total->SetBinContent(7 , *N_mcSIGNAL_abs );
  h_data_total->SetBinContent(7, *N_dataSIGNAL_abs);
  h_true_abs->SetBinContent(7,  *abs_absSignal );
  h_true_cex->SetBinContent(7,  *abs_chexSignal );
  h_true_nPi0->SetBinContent(7, *abs_nPi0Signal );
  h_true_BG->SetBinContent(7,   *abs_background );

 
  h_mc_total->SetFillColor(kBlue - 9);
  h_true_abs->SetFillColor(kGreen);
  h_true_cex->SetFillColor(kMagenta);
  h_true_nPi0->SetFillColor(kTeal);
  h_true_BG->SetFillColor(kBlue);

  //Scaling MC to DATA
  double total = h_data_total->GetBinContent(3);
  double total_2 = h_true_abs->GetBinContent(3) + h_true_cex->GetBinContent(3) +
                   h_true_nPi0->GetBinContent(3) + h_true_BG->GetBinContent(3);

  h_true_BG->Scale(total / total_2);
  h_true_abs->Scale(total / total_2);
  h_true_cex->Scale(total / total_2);
  h_true_nPi0->Scale(total / total_2);

  h_true_BG->Sumw2(0);
  h_true_abs->Sumw2(0);
  h_true_cex->Sumw2(0);
  h_true_nPi0->Sumw2(0);

  stack_cutFlow->Add(h_true_BG);
  stack_cutFlow->Add(h_true_abs);
  stack_cutFlow->Add(h_true_cex);
  stack_cutFlow->Add(h_true_nPi0);   



  //****Scaling Data to MC --> shouldn't do that***/

  //double total = h_true_abs->GetBinContent(5) + h_true_cex->GetBinContent(5) +
  //               h_true_nPi0->GetBinContent(5) + h_true_BG->GetBinContent(5);
  //h_data_total->Scale(total / h_data_total->GetBinContent(5));


  auto c1 = new TCanvas("Event Selection Flow", "",1600,1200);
  stack_cutFlow->Draw();
  h_data_total->SetMarkerStyle(20);
  h_data_total->Draw("P e1 SAME");
  c1->BuildLegend();
  c1->Write();
  c1->Close();

  auto c2 = new TCanvas("zoomed_flow", "", 1600, 1200);
  stack_cutFlow->Draw();
  stack_cutFlow->GetHistogram()->GetXaxis()->SetRangeUser(3, 7);
  stack_cutFlow->GetXaxis()->SetBinLabel(4, "Primary Pion");
  stack_cutFlow->GetXaxis()->SetBinLabel(5, "Abs + Cex + n#pi^{0}");
  stack_cutFlow->GetXaxis()->SetBinLabel(6, "Cex + n#pi^{0}");
  stack_cutFlow->GetXaxis()->SetBinLabel(7, "Abs");
  stack_cutFlow->GetXaxis()->SetLabelSize(.04);
  stack_cutFlow->Draw();
  h_data_total->Draw(" P e1 same");
  c2->BuildLegend();
  c2->Write();
  
  

  //*******************************
  //Efficiencies and Purity available Events are after endAPA3 CUT
  //*******************************

  auto comb_signal_help = (double)*mcCUT_noPionDaughter.Filter("true_combinedSignal == 1").Count();
  auto abs_signal_help = (double)*abs_absSignal;
  auto cex_signal_help = (double)(*chex_chexSignal + *chex_nPi0Signal) ;


  auto CUTendAPA3_combined_eff = 100* comb_signal_help / (double)*mcCUT_endAPA3.Filter("true_combinedSignal == 1").Count();
  auto CUTendAPA3_combined_pur = 100* comb_signal_help/ (double)*N_mcCUT_noPionDaughter;

  auto CUTendAPA3_abs_eff = 100* abs_signal_help / (double)*endAPA3_absSignal;
  auto CUTendAPA3_abs_pur = 100* abs_signal_help / (double)*N_mcSIGNAL_abs;

  auto CUTendAPA3_cex_eff = 100* cex_signal_help / (double)(*endAPA3_chexSignal + *endAPA3_nPi0Signal);
  auto CUTendAPA3_cex_pur = 100* cex_signal_help / (double)*N_mcSIGNAL_cex;


  //*******************
  //Efficiency and Purities after each Cut for Abs
  //******************
  //
  //TGraph for Purity Efficiency Development
  //
  int n_cuts = 3;
  double x[] = {1,2,3};
  double eff_abs[n_cuts];
  double pur_abs[n_cuts];
  double eff_times_pur_abs[n_cuts];


  pur_abs[1] = (double)*endAPA3_absSignal / (double)*N_mcCUT_endAPA3;
  eff_abs[1] = (double)*endAPA3_absSignal / (double)*endAPA3_absSignal;

  pur_abs[2] = (double)*noPionDaughter_absSignal / (double)*N_mcCUT_noPionDaughter;
  eff_abs[2] = (double)*noPionDaughter_absSignal / (double)*endAPA3_absSignal;

  pur_abs[3] = (double)*abs_absSignal / (double)*N_mcSIGNAL_abs;
  eff_abs[3] = (double)*abs_absSignal / (double)*endAPA3_absSignal;

  for(int i=0; i < n_cuts; i++){
     eff_times_pur_abs[i] = pur_abs[i]*eff_abs[i];
  };

 
  auto c_eff_pur_abs = new TCanvas("eff_pur_cuts", "",1600,2000 );
  c_eff_pur_abs->SetGrid();
  auto multi_graph = new TMultiGraph();
  multi_graph->SetTitle("Efficiency and Purity");
    
  auto graph_eff_abs = new TGraph(n_cuts,x,eff_abs);
  graph_eff_abs->SetTitle("Efficiency Absorption");
  graph_eff_abs->SetLineColor(kBlue);
  graph_eff_abs->SetLineWidth(4);
  graph_eff_abs->SetMarkerStyle(21);

  auto graph_pur_abs = new TGraph(n_cuts,x,pur_abs);
  graph_pur_abs->SetTitle("Purity Absorption");
  graph_pur_abs->SetLineColor(6);
  graph_pur_abs->SetLineWidth(4);
  graph_pur_abs->SetMarkerStyle(21);

  auto graph_eff_times_pur = new TGraph(n_cuts,x,eff_times_pur_abs);
  graph_eff_times_pur->SetTitle("Efficiency x Purity");
  graph_eff_times_pur->SetLineColor(8);
  graph_eff_times_pur->SetLineWidth(4);
  graph_eff_times_pur->SetMarkerStyle(29);

  multi_graph->Add(graph_eff_abs);
  multi_graph->Add(graph_pur_abs);
  multi_graph->Add(graph_eff_times_pur);

  multi_graph->GetYaxis()->SetRangeUser(0.,1.);
  multi_graph->GetXaxis()->SetTitle("Cuts");

  multi_graph->Draw("ALP");
    
  c_eff_pur_abs->BuildLegend();
  c_eff_pur_abs->Write();
  c_eff_pur_abs->Close();


  if( doCounting ){
    //Output Counting  //
    //

    std::cout << "Event Selection Cuts " << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "Total MC Events = " << *n_mc_all << std::endl;
    std::cout << "Total Data Events = " << *n_data_all << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "TRUE Event types" << std::endl;

    auto n_true_primPionInel = mc_all.Filter("true_primPionInel").Count();
    auto n_true_primPionInel_withElastic = mc_all.Filter(
        "true_primPionInel_withElastic").Count();
    std::cout << "Primary Pion Inelastic with Elastic before = " << *n_true_primPionInel_withElastic << std::endl;
    std::cout << "Events without primary Pion = " << *n_mc_all - *n_true_primPionInel << std::endl;
    std::cout << "Events with true Pion Daughters = " << *mc_all.Filter("true_pion_daughter > 0").Count() << std::endl;
    std::cout << "Combined Signal = " << *n_true_combinedSignal << std::endl;
    std::cout << "True Abs Signal = " << *n_true_absSignal << std::endl;
    std::cout << "True Cex Signal = " << *n_true_cexSignal << std::endl;
    std::cout << "True N-Pi0 Signal = " << *n_true_nPi0Signal << std::endl;
    std::cout << "True BackGround = " << *n_true_backGround << std::endl;
    std::cout << "********************************" << std::endl;

    std::cout << "Starting the CUTs on Total MC Events = " << *n_mc_all << std::endl;
    std::cout << "Starting the CUTs on Total DATA Events = " << *n_data_all << std::endl;
    std::cout << "---------PIMARY PARTICLE----" << std::endl;
    std::cout << "CUT 1 = Beam Particle Track Like  = " << *N_mcCUT_beamType << std::endl;
    std::cout << "DATA CUT 1 = " << *N_dataCUT_beamType << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcCUT_beamType.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *beamType_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << *beamType_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *beamType_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << *beamType_background << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_beamType.Filter("true_backGround").Filter("true_primPionInel == 0 ").Count() << std::endl;

    std::cout << std::endl;
    std::cout << "CUT 2 = Beam Position  = " << *N_mcCUT_beamCut << std::endl;
    std::cout << "DATA CUT 2 = " << *N_dataCUT_beamCut << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcCUT_beamCut.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *beamCut_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << *beamCut_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *beamCut_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << *beamCut_background << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_beamCut.Filter("true_backGround").Filter("true_primPionInel == 0").Count() << std::endl;

    std::cout << std::endl;
    std::cout << "CUT 3 = Primary in APA 3  = " << *N_mcCUT_endAPA3 << std::endl;
    std::cout << "DATA CUT 3 = " << *N_dataCUT_endAPA3 << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcCUT_endAPA3.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *endAPA3_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << *endAPA3_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *endAPA3_nPi0Signal << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_endAPA3.Filter("true_primPionInel == 0 ").Count() << std::endl;
    std::cout << "--------- True BackGround = " << *endAPA3_background << std::endl;
    std::cout << "--------- Events with true Pion Daughters = " << *mcCUT_endAPA3.Filter("true_backGround").Filter("true_pion_daughter > 0").Count() << std::endl;


    std::cout << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "RECO COMBINED EVENTS " << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "CUT 4 = Reject Pion Like Track Objects Chi2  = " << *N_mcCUT_noPionDaughter << std::endl;
    std::cout << "DATA CUT 4 = " << *N_dataCUT_noPionDaughter << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcCUT_noPionDaughter.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *mcCUT_noPionDaughter.Filter("true_absSignal == 1").Count() << std::endl;
    std::cout << "--------- True Chex  Signal = " << *mcCUT_noPionDaughter.Filter("true_chexSignal == 1").Count() << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *mcCUT_noPionDaughter.Filter("true_nPi0Signal == 1").Count() << std::endl;
    std::cout << "--------- True BackGround = " << *mcCUT_noPionDaughter.Filter("true_backGround == 1").Count() << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcCUT_noPionDaughter.Filter("true_primPionInel == 0").Count() << std::endl;
    std::cout << "--------- Contamination of Events with Pion Daughter = " << *mcCUT_noPionDaughter.Filter("true_backGround").Filter("true_pion_daughter > 0").Count() << std::endl;

    std::cout << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "RECO Charge Exchange EVENTS " << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "Events with Pi0 Shower  = " << *N_mcSIGNAL_cex << std::endl;
    std::cout << "DATA  = " << *N_dataSIGNAL_cex << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcSIGNAL_cex.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *chex_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << *chex_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *chex_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << *chex_background << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcSIGNAL_cex.Filter("true_primPionInel == 0 ").Count() << std::endl;
    std::cout << "--------- Contamination of Events with Pion Daughter = " << *mcSIGNAL_cex.Filter("true_backGround").Filter("true_pion_daughter > 0").Count() << std::endl;

    std::cout << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "RECO Absorption EVENTS " << std::endl;
    std::cout << "********************************" << std::endl;
    std::cout << "Events without Pi0 Shower  = " << *N_mcSIGNAL_abs << std::endl;
    std::cout << "DATA = " << *N_dataSIGNAL_abs << std::endl;
    std::cout << "--------- True Combined Signal = " << *mcSIGNAL_abs.Filter("true_combinedSignal == 1").Count() << std::endl;
    std::cout << "--------- True Absorption  Signal = " << *abs_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << *abs_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << *abs_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << *abs_background << std::endl;
    std::cout << "--------- Contamination of primary NON-pions = " << *mcSIGNAL_abs.Filter("true_primPionInel == 0").Count() << std::endl;
    std::cout << "--------- Contamination of Events with Pion Daughter = " << *mcSIGNAL_abs.Filter("true_backGround").Filter("true_pion_daughter > 0").Count() << std::endl;

  

  }

  std::cout << std::endl;
  std::cout << "********************************" << std::endl;
  std::cout << "EFFICIENCY & PURITY " << std::endl;
  std::cout << "********************************" << std::endl;
  std::cout << "**SELECTION********EFFICIENCY to Chi2 on primary*****************PURITY***********" << std::endl;
  std::cout << "*******************************************************" << std::endl;
  std::cout << "--COMBINED------" <<  CUTendAPA3_combined_eff << "--------------------------" << CUTendAPA3_combined_pur << std::endl;
  std::cout << "-----CEX*--------" <<  CUTendAPA3_cex_eff << "--------------------------" << CUTendAPA3_cex_pur << std::endl;
  std::cout << "-----ABS--------" <<  CUTendAPA3_abs_eff << "---------------------------" << CUTendAPA3_abs_pur << std::endl;
  std::cout << "********************************" << std::endl;

  c1->Write();
  stack_cutFlow->Write();
  h_data_total->Write();
  output->Write();

  return 0;
}
