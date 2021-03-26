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

const std::string mc_all_file = "test055_eventSelection_mc_all.root";
const std::string data_all_file = "test055_eventSelection_data_all.root";

int eventSelection_cutEffPur() {

  //gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");

  TFile *output = new TFile ("eventSelection_cutEffPur_test04MuRemove.root", "RECREATE");
  THStack *stack_cutFlow = new THStack("cutFlow", "Cut Flow MC and Data");

  //************** CUT Flow Histogram ************//
  TH1I *h_mc_total = new TH1I("mc_total", "MC selected", 8 , 0., 8.);
  TH1I *h_data_total = new TH1I("data_total", "Data selected", 8 , 0., 8.);
  TH1I *h_true_abs = new TH1I("true_abs", "True Abs Signal", 8 , 0., 8.);
  TH1I *h_true_cex = new TH1I("true_cex", "True Cex Signal", 8 , 0., 8.);
  TH1I *h_true_nPi0 = new TH1I("true_nPi0", "True N-Pi0 Signal", 8 , 0., 8.);
  TH1I *h_true_BG = new TH1I("true_BG", "Background", 8 , 0., 8.);
  TH1I *h_true_primMu = new TH1I("true_primMu", "primaryMuons", 8 , 0., 8.);


  //Before filling, split it up MC by the true signal/BG
  //This will later be used in the counting as well
  
  ROOT::RDataFrame mc_all(pionTree, mc_all_file);
  auto mc_trackLike = mc_all.Filter("primary_isBeamType");
  auto mc_beamCut = mc_trackLike.Filter("passBeamCut && passBeamCutBI");
  auto mc_removeMu = mc_beamCut.Filter("!isPrimaryMuonCandidate");
  //auto mc_removeMu = mc_beamCut;
  //until here is inicident Pion sample
  auto mc_endAPA3 = mc_removeMu.Filter("primary_ends_inAPA3");
  auto mc_combined = mc_endAPA3.Filter("has_noPion_daughter");
  
  auto mc_abs = mc_combined.Filter("!has_shower_nHits_distance");
  auto mc_cex = mc_combined.Filter("has_shower_nHits_distance");


  ROOT::RDataFrame data_all(pionTree, data_all_file);
  auto data_trackLike = data_all.Filter("primary_isBeamType");
  auto data_beamCut = data_trackLike.Filter("passBeamQuality && passBeamCut");
  auto data_removeMu = data_beamCut.Filter("!isPrimaryMuonCandidate");
  //not removingMu
  //auto data_removeMu = data_beamCut;
  //until here is inicident Pion sample
  auto data_endAPA3 = data_removeMu.Filter("primary_ends_inAPA3");
  auto data_combined = data_endAPA3.Filter("has_noPion_daughter");
  
  auto data_abs = data_combined.Filter("!has_shower_nHits_distance");
  auto data_cex = data_combined.Filter("has_shower_nHits_distance");
  
  double n_mc_all = (double) *mc_all.Count();
  double n_data_all = (double) *data_all.Count();
  double n_true_combinedSignal = (double) *mc_all.Filter("true_combinedSignal").Count();
  double n_true_absSignal = (double) *mc_all.Filter("true_absSignal").Count();
  double n_true_cexSignal = (double) *mc_all.Filter("true_chexSignal").Count();
  double n_true_nPi0Signal = (double) *mc_all.Filter("true_nPi0Signal").Count();
  double n_true_backGround = (double) *mc_all.Filter("true_backGround && !primaryMuon").Count();
  double n_true_primPi = (double) *mc_all.Filter("true_primPionInel").Count();
  double n_primaryMuon = (double) *mc_all.Filter("primaryMuon").Count();

  h_mc_total->SetBinContent(1 , n_mc_all );
  h_data_total->SetBinContent(1, n_data_all );
  h_true_abs->SetBinContent(1, n_true_absSignal );
  h_true_cex->SetBinContent(1, n_true_cexSignal );
  h_true_nPi0->SetBinContent(1, n_true_nPi0Signal );
  h_true_BG->SetBinContent(1, n_true_backGround );
  h_true_primMu->SetBinContent(1, n_primaryMuon );

  double N_mc_trackLike = (double) *mc_trackLike.Count();
  double N_data_trackLike = (double) *data_trackLike.Count();
  double trackLike_absSignal = (double) *mc_trackLike.Filter("true_absSignal").Count();
  double trackLike_chexSignal = (double) *mc_trackLike.Filter("true_chexSignal").Count();
  double trackLike_nPi0Signal = (double) *mc_trackLike.Filter("true_nPi0Signal").Count();
  double trackLike_background = (double) *mc_trackLike.Filter("true_backGround && !primaryMuon").Count();
  double trackLike_primPi = (double) *mc_trackLike.Filter("true_primPionInel").Count();
  double trackLike_primaryMuon = (double) *mc_trackLike.Filter("primaryMuon").Count();

  h_mc_total->SetBinContent(2 , N_mc_trackLike );
  h_data_total->SetBinContent(2, N_data_trackLike);
  h_true_abs->SetBinContent(2, trackLike_absSignal);
  h_true_cex->SetBinContent(2, trackLike_chexSignal);
  h_true_nPi0->SetBinContent(2, trackLike_nPi0Signal);
  h_true_BG->SetBinContent(2, trackLike_background);
  h_true_primMu->SetBinContent(2, trackLike_primaryMuon );
  

  double N_mc_beamCut = (double) *mc_beamCut.Count();
  double N_data_beamCut = (double) *data_beamCut.Count();
  double beamCut_absSignal = (double) *mc_beamCut.Filter("true_absSignal").Count();
  double beamCut_chexSignal = (double) *mc_beamCut.Filter("true_chexSignal").Count();
  double beamCut_nPi0Signal = (double) *mc_beamCut.Filter("true_nPi0Signal").Count();
  double beamCut_background = (double) *mc_beamCut.Filter("true_backGround && !primaryMuon").Count();
  double beamCut_primPi = (double) *mc_beamCut.Filter("true_primPionInel").Count();
  double beamCut_primaryMuon = (double) *mc_beamCut.Filter("primaryMuon").Count();

  h_mc_total->SetBinContent(3 , N_mc_beamCut );
  h_data_total->SetBinContent(3, N_data_beamCut);
  h_true_abs->SetBinContent(3, beamCut_absSignal);
  h_true_cex->SetBinContent(3, beamCut_chexSignal);
  h_true_nPi0->SetBinContent(3, beamCut_nPi0Signal);
  h_true_BG->SetBinContent(3, beamCut_background);
  h_true_primMu->SetBinContent(3, beamCut_primaryMuon );

  double N_mc_removeMu = (double) *mc_removeMu.Count();
  double N_data_removeMu = (double) *data_removeMu.Count();
  double removeMu_absSignal = (double) *mc_removeMu.Filter("true_absSignal").Count();
  double removeMu_chexSignal = (double) *mc_removeMu.Filter("true_chexSignal").Count();
  double removeMu_nPi0Signal = (double) *mc_removeMu.Filter("true_nPi0Signal").Count();
  double removeMu_background = (double) *mc_removeMu.Filter("true_backGround && !primaryMuon").Count();
  double removeMu_primPi = (double) *mc_removeMu.Filter("true_primPionInel").Count();
  double removeMu_primaryMuon = (double) *mc_removeMu.Filter("primaryMuon").Count();

  h_mc_total->SetBinContent(4 , N_mc_removeMu );
  h_data_total->SetBinContent(4, N_data_removeMu);
  h_true_abs->SetBinContent(4, removeMu_absSignal);
  h_true_cex->SetBinContent(4, removeMu_chexSignal);
  h_true_nPi0->SetBinContent(4, removeMu_nPi0Signal);
  h_true_BG->SetBinContent(4, removeMu_background);
  h_true_primMu->SetBinContent(4, removeMu_primaryMuon );

  double N_mc_endAPA3 = (double) *mc_endAPA3.Count();
  double N_data_endAPA3 = (double) *data_endAPA3.Count();
  double endAPA3_absSignal = (double) * mc_endAPA3.Filter("true_absSignal").Count();
  double endAPA3_chexSignal = (double) *mc_endAPA3.Filter("true_chexSignal").Count();
  double endAPA3_nPi0Signal = (double) *mc_endAPA3.Filter("true_nPi0Signal").Count();
  double endAPA3_background = (double) *mc_endAPA3.Filter("true_backGround && !primaryMuon").Count();
  double endAPA3_primaryMuon = (double) *mc_endAPA3.Filter("primaryMuon").Count();

  h_mc_total->SetBinContent(5 , N_mc_endAPA3 );
  h_data_total->SetBinContent(5, N_data_endAPA3);
  h_true_abs->SetBinContent(5, endAPA3_absSignal );
  h_true_cex->SetBinContent(5, endAPA3_chexSignal );
  h_true_nPi0->SetBinContent(5, endAPA3_nPi0Signal );
  h_true_BG->SetBinContent(5, endAPA3_background );
  h_true_primMu->SetBinContent(5, endAPA3_primaryMuon );

  double N_mc_combined = (double) *mc_combined.Count();
  double N_data_combined = (double) *data_combined.Count();
  double combined_absSignal = (double) *mc_combined.Filter("true_absSignal").Count();
  double combined_chexSignal = (double) *mc_combined.Filter("true_chexSignal").Count();
  double combined_nPi0Signal = (double) *mc_combined.Filter("true_nPi0Signal").Count();
  double combined_background = (double) *mc_combined.Filter("true_backGround && !primaryMuon").Count();
  double combined_primaryMuon = (double) *mc_combined.Filter("primaryMuon").Count();

  h_mc_total->SetBinContent(6 , N_mc_combined );
  h_data_total->SetBinContent(6, N_data_combined);
  h_true_abs->SetBinContent(6, combined_absSignal );
  h_true_cex->SetBinContent(6, combined_chexSignal );
  h_true_nPi0->SetBinContent(6, combined_nPi0Signal );
  h_true_BG->SetBinContent(6, combined_background );
  h_true_primMu->SetBinContent(6, combined_primaryMuon );

  double N_mc_abs = (double) *mc_abs.Count();
  double N_data_abs = (double) *data_abs.Count();
  double abs_absSignal = (double) *mc_abs.Filter("true_absSignal").Count();
  double abs_chexSignal = (double) *mc_abs.Filter("true_chexSignal").Count();
  double abs_nPi0Signal = (double) *mc_abs.Filter("true_nPi0Signal").Count();
  double abs_background = (double) *mc_abs.Filter("true_backGround && !primaryMuon").Count();
  double abs_primaryMuon = (double) *mc_abs.Filter("primaryMuon").Count();

  h_mc_total->SetBinContent(7 , N_mc_abs );
  h_data_total->SetBinContent(7, N_data_abs);
  h_true_abs->SetBinContent(7,  abs_absSignal );
  h_true_cex->SetBinContent(7,  abs_chexSignal );
  h_true_nPi0->SetBinContent(7, abs_nPi0Signal );
  h_true_BG->SetBinContent(7,   abs_background );
  h_true_primMu->SetBinContent(7, abs_primaryMuon );

  double N_mc_cex = (double) *mc_cex.Count();
  double N_data_cex = (double) *data_cex.Count();
  double cex_absSignal = (double) *mc_cex.Filter("true_absSignal").Count();
  double cex_chexSignal = (double) *mc_cex.Filter("true_chexSignal").Count();
  double cex_nPi0Signal = (double) *mc_cex.Filter("true_nPi0Signal").Count();
  double cex_background = (double) *mc_cex.Filter("true_backGround && !primaryMuon").Count();
  double cex_primaryMuon = (double) *mc_cex.Filter("primaryMuon").Count();

  h_mc_total->SetBinContent(8 , N_mc_cex );
  h_data_total->SetBinContent(8, N_data_cex);
  h_true_abs->SetBinContent(8,  cex_absSignal );
  h_true_cex->SetBinContent(8,  cex_chexSignal );
  h_true_nPi0->SetBinContent(8, cex_nPi0Signal );
  h_true_BG->SetBinContent(8,   cex_background );
  h_true_primMu->SetBinContent(8, cex_primaryMuon );
 
  h_mc_total->SetFillColor(kBlue - 9);
  h_true_abs->SetFillColor(kGreen);
  h_true_cex->SetFillColor(kMagenta);
  h_true_nPi0->SetFillColor(kTeal);
  h_true_BG->SetFillColor(kBlue);
  h_true_primMu->SetFillColor(kAzure+1);


  //Scaling MC to DATA after the beamQuality Cuts
  double total = h_data_total->GetBinContent(3);
  double total_2 = h_true_abs->GetBinContent(3) + h_true_cex->GetBinContent(3) +
                   h_true_nPi0->GetBinContent(3) + h_true_BG->GetBinContent(3) + h_true_primMu->GetBinContent(3);

  h_true_primMu->Scale(total / total_2);
  h_true_BG->Scale(total / total_2);
  h_true_abs->Scale(total / total_2);
  h_true_cex->Scale(total / total_2);
  h_true_nPi0->Scale(total / total_2);

  h_true_primMu->Sumw2(0);
  h_true_BG->Sumw2(0);
  h_true_abs->Sumw2(0);
  h_true_cex->Sumw2(0);
  h_true_nPi0->Sumw2(0);

  
  stack_cutFlow->Add(h_true_primMu);
  stack_cutFlow->Add(h_true_BG);
  stack_cutFlow->Add(h_true_abs);
  stack_cutFlow->Add(h_true_cex);
  stack_cutFlow->Add(h_true_nPi0);   


  auto c1 = new TCanvas("Event Selection Flow", "",1600,1200);
  //stack_cutFlow->GetXaxis()->SetRangeUser(3,8);
  //h_data_total->GetXaxis()->SetRangeUser(3,8);
  stack_cutFlow->Draw();
  h_data_total->Draw("PSAME");
  c1->BuildLegend();
  c1->Write();


  //c1->Close();

  

  
  //------------------------------
  //Efficiencies and Purity 
  //------------------------------
  //efficiency cut after cut 
  //and also overall efficiency wrt to available events after beamCuts
  //
  double eff_incidentPi_trackLike = 100 * ( trackLike_primPi / n_true_primPi );
  double pur_incidentPi_trackLike = 100 * ( trackLike_primPi / N_mc_trackLike );
  
  double eff_incidentPi_beamCut = 100 * ( beamCut_primPi / trackLike_primPi );
  double pur_incidentPi_beamCut = 100 * ( beamCut_primPi / N_mc_beamCut );
  
  double eff_incidentPi_removeMu = 100 * ( removeMu_primPi / beamCut_primPi );
  double pur_incidentPi_removeMu = 100 * ( removeMu_primPi / N_mc_removeMu );
  double effGeneral_incidentPi_removeMu = 100 * ( removeMu_primPi / beamCut_primPi );

  //this is now the incident Pion Sample

  double eff_abs_removeMu = 100 * ( removeMu_absSignal / beamCut_absSignal );
  double pur_abs_removeMu = 100 * ( removeMu_absSignal / N_mc_removeMu );
  double effGeneral_abs_removeMu = 100 * ( removeMu_absSignal / beamCut_absSignal );
  
  double eff_cex_removeMu = 100 * ( removeMu_chexSignal / beamCut_chexSignal );
  double pur_cex_removeMu = 100 * ( removeMu_chexSignal / N_mc_removeMu );
  double effGeneral_cex_removeMu = 100 * ( removeMu_chexSignal / beamCut_chexSignal );
  
  double eff_abs_endAPA3 = 100 * ( endAPA3_absSignal / removeMu_absSignal );
  double pur_abs_endAPA3 = 100 * ( endAPA3_absSignal / N_mc_endAPA3 );
  double effGeneral_abs_endAPA3 = 100 * ( endAPA3_absSignal / beamCut_absSignal );
  
  double eff_cex_endAPA3 = 100 * ( endAPA3_chexSignal / removeMu_chexSignal );
  double pur_cex_endAPA3 = 100 * ( endAPA3_chexSignal / N_mc_endAPA3 );
  double effGeneral_cex_endAPA3 = 100 * ( endAPA3_chexSignal / beamCut_chexSignal );
  
  double eff_abs_combined = 100 * ( combined_absSignal / endAPA3_absSignal );
  double pur_abs_combined = 100 * ( combined_absSignal / N_mc_combined );
  double effGeneral_abs_combined = 100 * ( combined_absSignal / beamCut_absSignal );
  
  double eff_cex_combined = 100 * ( combined_chexSignal / endAPA3_chexSignal );
  double pur_cex_combined = 100 * ( combined_chexSignal / N_mc_combined );
  double effGeneral_cex_combined = 100 * ( combined_chexSignal / beamCut_chexSignal );
  
  double eff_abs_abs = 100 * ( abs_absSignal / combined_absSignal );
  double pur_abs_abs = 100 * ( abs_absSignal / N_mc_abs );
  double effGeneral_abs_abs = 100 * ( abs_absSignal / beamCut_absSignal );
  
  double eff_cex_cex = 100 * ( cex_chexSignal / combined_chexSignal );
  double pur_cex_cex = 100 * ( cex_chexSignal / N_mc_cex );
  double effGeneral_cex_cex = 100 * ( cex_chexSignal / beamCut_chexSignal );

  //------------------------------
  //Efficiency and Purities after each Cut for Abs
  //------------------------------
  //
  //TGraph for Purity Efficiency Development
  //
  int n_cuts = 4;
  double x[] = {1,2,3,4};
  double eff_abs[n_cuts];
  double effGeneral_abs[n_cuts];
  double pur_abs[n_cuts];
  double eff_times_pur_abs[n_cuts];

  pur_abs[0] = pur_abs_removeMu / 100.;
  eff_abs[0] = eff_abs_removeMu / 100.;
  effGeneral_abs[0] = effGeneral_abs_removeMu / 100.;
  
  pur_abs[1] = pur_abs_endAPA3 / 100.;
  eff_abs[1] = eff_abs_endAPA3 / 100.;
  effGeneral_abs[1] = effGeneral_abs_endAPA3 / 100.;
  //std::cout << "effGeneral_abs = " << effGeneral_abs[] << std::endl;

  pur_abs[2] = pur_abs_combined / 100.;
  eff_abs[2] = eff_abs_combined / 100.;
  effGeneral_abs[2] = effGeneral_abs_combined / 100.;

  pur_abs[3] = pur_abs_abs / 100.;
  eff_abs[3] = eff_abs_abs / 100.;
  effGeneral_abs[3] = effGeneral_abs_abs / 100.;
  
 
  for(int i=0; i < n_cuts; i++){
     eff_times_pur_abs[i] = pur_abs[i]*effGeneral_abs[i];
  };


  auto c_eff_pur_abs = new TCanvas("eff_pur_cuts", "",1600,2000 );
  c_eff_pur_abs->SetGrid();
  auto multi_graph = new TMultiGraph();
  multi_graph->SetTitle("Efficiency and Purity");
    
  auto graph_eff_abs = new TGraph(n_cuts,x,eff_abs);
  graph_eff_abs->SetTitle("Efficiency2 = (NtrueAbsAtCut / NtrueAbsBeforeCut)");
  graph_eff_abs->SetLineColor(kBlue);
  graph_eff_abs->SetLineWidth(4);
  graph_eff_abs->SetMarkerStyle(21);
  
  auto graph_effGeneral_abs = new TGraph(n_cuts,x,effGeneral_abs);
  graph_effGeneral_abs->SetTitle("Efficiency1 = (NtrueAbsAtCut / NtrueAbsAfterBeamCut)");
  graph_effGeneral_abs->SetLineColor(kRed);
  graph_effGeneral_abs->SetLineWidth(4);
  graph_effGeneral_abs->SetMarkerStyle(21);

  auto graph_pur_abs = new TGraph(n_cuts,x,pur_abs);
  graph_pur_abs->SetTitle("Purity Absorption");
  graph_pur_abs->SetLineColor(6);
  graph_pur_abs->SetLineWidth(4);
  graph_pur_abs->SetMarkerStyle(21);

  auto graph_eff_times_pur = new TGraph(n_cuts,x,eff_times_pur_abs);
  graph_eff_times_pur->SetTitle("Efficiency1 x Purity");
  graph_eff_times_pur->SetLineColor(8);
  graph_eff_times_pur->SetLineWidth(4);
  graph_eff_times_pur->SetMarkerStyle(29);

  multi_graph->Add(graph_effGeneral_abs);
  multi_graph->Add(graph_eff_abs);
  multi_graph->Add(graph_pur_abs);
  multi_graph->Add(graph_eff_times_pur);

  multi_graph->GetYaxis()->SetRangeUser(0.,1.);
  multi_graph->GetXaxis()->SetTitle("Cuts");

  multi_graph->Draw("ALP");
    
  c_eff_pur_abs->BuildLegend();
  c_eff_pur_abs->Write();
  c_eff_pur_abs->Close();


      //Output Counting  //
    //

    std::cout << "Event Selection Cuts " << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "Total MC Events = " << n_mc_all << std::endl;
    std::cout << "Total Data Events = " << n_data_all << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "TRUE Event types" << std::endl;

    std::cout << "Primary Pion Inelastic = " << n_true_primPi << std::endl;
    std::cout << "Events with PrimaryMuon = " << n_primaryMuon << std::endl;
    std::cout << "Combined Signal = " << n_true_combinedSignal << std::endl;
    std::cout << "True Abs Signal = " << n_true_absSignal << std::endl;
    std::cout << "True Cex Signal = " << n_true_cexSignal << std::endl;
    std::cout << "True N-Pi0 Signal = " << n_true_nPi0Signal << std::endl;
    std::cout << "True BackGround = " << n_true_backGround << std::endl;
    std::cout << "------------------------------" << std::endl;

    std::cout << "Starting the CUTs on Total MC Events = " << n_mc_all << std::endl;
    std::cout << "Starting the CUTs on Total DATA Events = " << n_data_all << std::endl;
    std::cout << "---------PIMARY PARTICLE----" << std::endl;
    std::cout << "CUT 1 = Beam Particle Track Like  = " << N_mc_trackLike << std::endl;
    std::cout << "DATA CUT 1 = " << N_data_trackLike << std::endl;
    std::cout << "--------- True Absorption  Signal = " << trackLike_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << trackLike_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << trackLike_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << trackLike_background << std::endl;
    std::cout << "--------- Contamination of primaryMuons= " << trackLike_primaryMuon << std::endl;

    std::cout << std::endl;
    std::cout << "CUT 2 = Beam Window Cut  = " << N_mc_beamCut << std::endl;
    std::cout << "DATA CUT 2 = " << N_data_beamCut << std::endl;
    std::cout << "--------- True Absorption  Signal = " << beamCut_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << beamCut_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << beamCut_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << beamCut_background << std::endl;
    std::cout << "--------- Contamination of primaryMuons= " << beamCut_primaryMuon << std::endl;
    
    std::cout << std::endl;
    std::cout << "CUT 3 = Remove Primary Muons  = " << N_mc_removeMu << std::endl;
    std::cout << "DATA CUT 3 = " << N_data_removeMu << std::endl;
    std::cout << "--------- True Absorption  Signal = " << removeMu_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << removeMu_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << removeMu_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << removeMu_background << std::endl;
    std::cout << "--------- Contamination of primaryMuons= " << removeMu_primaryMuon << std::endl;
    
    std::cout << std::endl;
    std::cout << "CUT 4 = end APA3 = " << N_mc_endAPA3 << std::endl;
    std::cout << "DATA CUT 4 = " << N_data_endAPA3 << std::endl;
    std::cout << "--------- True Absorption  Signal = " << endAPA3_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << endAPA3_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << endAPA3_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << endAPA3_background << std::endl;
    std::cout << "--------- Contamination of primaryMuons= " << endAPA3_primaryMuon << std::endl;

    std::cout << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "RECO COMBINED EVENTS " << std::endl;
    std::cout << "------------------------------" << std::endl;

    std::cout << "CUT 5 = remove Pion Daughters = " << N_mc_combined << std::endl;
    std::cout << "DATA CUT 5 = " << N_data_combined << std::endl;
    std::cout << "--------- True Absorption  Signal = " << combined_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << combined_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << combined_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << combined_background << std::endl;
    std::cout << "--------- Contamination of primaryMuons= " << combined_primaryMuon << std::endl;
    
    std::cout << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "RECO Charge Exchange EVENTS " << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "CUT  = find Showers = " << N_mc_cex << std::endl;
    std::cout << "DATA CUT  = " << N_data_cex << std::endl;
    std::cout << "--------- True Absorption  Signal = " << cex_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << cex_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << cex_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << cex_background << std::endl;
    std::cout << "--------- Contamination of primaryMuons= " << cex_primaryMuon << std::endl;

    std::cout << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "RECO Absorption EVENTS " << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "CUT  = reject Showers = " << N_mc_abs << std::endl;
    std::cout << "DATA CUT  = " << N_data_abs << std::endl;
    std::cout << "--------- True Absorption  Signal = " << abs_absSignal << std::endl;
    std::cout << "--------- True Chex  Signal = " << abs_chexSignal << std::endl;
    std::cout << "--------- True N-Pi0  Signal = " << abs_nPi0Signal << std::endl;
    std::cout << "--------- True BackGround = " << abs_background << std::endl;
    std::cout << "--------- Contamination of primaryMuons= " << abs_primaryMuon << std::endl;

  

  std::cout << std::endl;
  std::cout << "------------------------------" << std::endl;
  std::cout << "EFFICIENCY & PURITY " << std::endl;
  std::cout << "------------------------------" << std::endl;
  std::cout << "------SELECTION--------------EFFICIENCY--------------PURITY-------" << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "--Incident Pions----" <<  effGeneral_incidentPi_removeMu << "-------------------" << pur_incidentPi_removeMu << std::endl;

  std::cout << "--CEX--------" <<  effGeneral_cex_cex<< "--------------------------" << pur_cex_cex << std::endl;
  std::cout << "--ABS--------" <<  effGeneral_abs_abs << "---------------------------" << pur_abs_abs << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  //c1->Write();
  //stack_cutFlow->Write();
  //h_data_total->Write();
 
  output->Write();

  return 0;
}
