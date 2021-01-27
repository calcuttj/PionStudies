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

using namespace ROOT::VecOps;

int leading_proton_study(const string selected_abs_file, const string APA3_cut_file,
                         bool doBatch = false) {
  gROOT->SetBatch(doBatch);
  gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");

  TFile *output = new TFile ("leading_proton.root", "RECREATE");

  ROOT::RDataFrame mcSIGNAL_abs(inputTree, selected_abs_file);
  ROOT::RDataFrame mcCUT_endAPA3(inputTree, APA3_cut_file); 
  auto endAPA3_absSignal_filter = mcCUT_endAPA3.Filter("true_absSignal");
  auto abs_absSignal_filter = mcSIGNAL_abs.Filter("true_absSignal");

  abs_absSignal_filter = abs_absSignal_filter.Define(
      "leading_proton_momentum", leading_proton_momentum,
      {"true_beam_daughter_startP", "true_beam_daughter_PDG"})
      .Define("leading_proton_det_theta", leading_proton_det_theta,
              {"true_beam_daughter_startP", "true_beam_daughter_PDG",
              "true_beam_daughter_startPz"})
      .Define("leading_proton_det_phi", leading_proton_det_phi,
              {"true_beam_daughter_startP", "true_beam_daughter_PDG",
              "true_beam_daughter_startPx", "true_beam_daughter_startPy"});

  endAPA3_absSignal_filter = endAPA3_absSignal_filter.Define(
      "leading_proton_momentum", leading_proton_momentum,
      {"true_beam_daughter_startP", "true_beam_daughter_PDG"})
      .Define("leading_proton_det_theta", leading_proton_det_theta,
              {"true_beam_daughter_startP", "true_beam_daughter_PDG",
              "true_beam_daughter_startPz"})
      .Define("leading_proton_det_phi", leading_proton_det_phi,
              {"true_beam_daughter_startP", "true_beam_daughter_PDG",
              "true_beam_daughter_startPx", "true_beam_daughter_startPy"});

  TH2D total_abs_nProton_leadingP = *(endAPA3_absSignal_filter.Histo2D(
      {"total_abs_nProton_leadingP", "", 20, 0, 20, 50, 0., 2.},
       "true_daughter_nProton", "leading_proton_momentum"));
  total_abs_nProton_leadingP.Write();

  TH2D selected_abs_nProton_leadingP = *(abs_absSignal_filter.Histo2D(
      {"selected_abs_nProton_leadingP", "", 20, 0, 20, 50, 0., 2.},
       "true_daughter_nProton", "leading_proton_momentum"));
  selected_abs_nProton_leadingP.Write();

  TH1D * abs_selected_nProton =
      (TH1D*)selected_abs_nProton_leadingP.ProjectionX()->Clone();
  TH1D * abs_total_nProton =
      (TH1D*)total_abs_nProton_leadingP.ProjectionX()->Clone();
  TGraphAsymmErrors abs_eff_nProton(abs_selected_nProton, abs_total_nProton);
  abs_total_nProton->Write();
  abs_selected_nProton->Write();
  abs_eff_nProton.Write("abs_eff_nProton");
  
  TH2D selected_abs_nProton_det_theta = *(abs_absSignal_filter.Histo2D(
      {"selected_abs_nProton_nProton_det_theta", "", 20, 0, 20, 10, -1., 1.},
       "true_daughter_nProton", "leading_proton_det_theta"));
  TH2D total_abs_nProton_det_theta = *(endAPA3_absSignal_filter.Histo2D(
      {"total_abs_nProton_nProton_det_theta", "", 20, 0, 20, 10, -1., 1.},
       "true_daughter_nProton", "leading_proton_det_theta"));
  
  selected_abs_nProton_det_theta.Write();
  total_abs_nProton_det_theta.Write();
  TH2D * eff_abs_nProton_det_theta = (TH2D*)selected_abs_nProton_det_theta.Clone();
  eff_abs_nProton_det_theta->Divide(&total_abs_nProton_det_theta);
  eff_abs_nProton_det_theta->Write("eff_abs_nProton_det_theta");

  for (int i = 1; i <= total_abs_nProton_det_theta.GetNbinsX(); ++i) {
    TH1D * r = (TH1D*)selected_abs_nProton_det_theta.ProjectionY(
        "", i, i)->Clone();
    r->Divide(total_abs_nProton_det_theta.ProjectionY("", i, i));
    std::string name = "eff_nProton_det_theta_" + std::to_string(i);
    r->Write(name.c_str());
  }

  TH2D selected_abs_nProton_det_phi = *(abs_absSignal_filter.Histo2D(
      {"selected_abs_nProton_nProton_det_phi", "", 20, 0, 20, 10, -90., 90.},
       "true_daughter_nProton", "leading_proton_det_phi"));
  TH2D total_abs_nProton_det_phi = *(endAPA3_absSignal_filter.Histo2D(
      {"total_abs_nProton_nProton_det_phi", "", 20, 0, 20, 10, -90., 90.},
       "true_daughter_nProton", "leading_proton_det_phi"));
  
  selected_abs_nProton_det_phi.Write();
  total_abs_nProton_det_phi.Write();
  TH2D * eff_abs_nProton_det_phi = (TH2D*)selected_abs_nProton_det_phi.Clone();
  eff_abs_nProton_det_phi->Divide(&total_abs_nProton_det_phi);
  eff_abs_nProton_det_phi->Write("eff_abs_nProton_det_phi");

  return 0;
}


