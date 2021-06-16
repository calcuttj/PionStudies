#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TF1.h"
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
//Look at dE vs wire without SCE corrections
//Francescos Idea

//int dEdX_correction_try(const string path1, const string path2, const string path3, const string path4){
int dEdX_correction_try(const string path1, const string path2 ){

   auto dEdX_trkPitch = [](std::vector<double> &dEdX, std::vector<double> &pitch){

      std::vector<double> product;

      for(size_t i=0; i< pitch.size(); i++){

         product.push_back( dEdX[i] * pitch[i]);
      };

      product.resize(pitch.size());
      return product;
   };

   auto cut_dEdX = [](std::vector<double> &dEdX, std::vector<double> &wire){
      std::vector<double> cut = dEdX;
      cut.resize(wire.size());
      return cut;
   };



   //Int_t palette_sce[] = {kRed+2, kRed-3,  kRed-9, kBlue -9, kBlue-3 , kBlue+2};
   //gStyle->SetPalette(6,palette_sce);

   //gROOT->ForceStyle();
   //read in Data, MC SCE and MC without SCE
   ROOT::RDataFrame frame_mc_SCE(pionTree, path1);
   ROOT::RDataFrame frame_data5387(pionTree, path2);

   //ROOT::RDataFrame frame_data5809(pionTree, path3);
   //ROOT::RDataFrame frame_data5770(pionTree, path4);


   TFile *output = new TFile ("output_dEdX_correction.root", "RECREATE");
   int nBin_wire, nBin_dE, nBin_dEdX, nBin_pitch;
   double binLow_wire, binHigh_wire,  binLow_dE, binHigh_dE, 
          binLow_dEdX, binHigh_dEdX,  binLow_pitch, binHigh_pitch;

   binLow_wire = 0; binHigh_wire = 750; nBin_wire = 750;
   binLow_dE = 0.6; binHigh_dE = 2.2; nBin_dE = 100;
   binLow_dEdX = 0.6; binHigh_dEdX = 3.2; nBin_dEdX = 200;
   binLow_pitch = 0.47; binHigh_pitch = 0.75; nBin_pitch = 100;

   auto filtered_SCE = frame_mc_SCE
      .Filter("primary_isBeamType")
      .Filter("passBeamQuality_TPCjustPosition")
      .Filter("!isPrimaryMuonCandidate")
      .Define("uncorrected_dE", dEdX_trkPitch, {"reco_beam_calibrated_dEdX_NoSCE", "reco_beam_TrkPitch_NoSCE"})
      .Define("corrected_dE", dEdX_trkPitch, {"reco_beam_calibrated_dEdX_SCE", "reco_beam_TrkPitch_SCE"})
      .Define("cut_uncorrected_dEdX", cut_dEdX, {"reco_beam_calibrated_dEdX_NoSCE", "reco_beam_calo_wire"})
      .Define("cut_corrected_dEdX", cut_dEdX, {"reco_beam_calibrated_dEdX_SCE", "reco_beam_calo_wire"});

   auto filtered_data5387 = frame_data5387
      .Filter("primary_isBeamType")
      .Filter("passBeamQuality_TPCjustPosition")
      .Filter("!isPrimaryMuonCandidate")
      .Define("uncorrected_dE", dEdX_trkPitch, {"reco_beam_calibrated_dEdX_NoSCE", "reco_beam_TrkPitch_NoSCE"})
      .Define("cut_uncorrected_dEdX", cut_dEdX, {"reco_beam_calibrated_dEdX_NoSCE", "reco_beam_calo_wire_NoSCE"})
      .Define("corrected_dE", dEdX_trkPitch, {"reco_beam_calibrated_dEdX_SCE", "reco_beam_TrkPitch_SCE"})
      .Define("cut_corrected_dEdX", cut_dEdX, {"reco_beam_calibrated_dEdX_SCE", "reco_beam_calo_wire"});

/*   
   auto filtered_data5809 = frame_data5809
      .Filter("primary_isBeamType")
      .Filter("passBeamCut")
      .Filter("!isPrimaryMuonCandidate")
      .Alias("reco_beam_calibrated_dEdX_SCE", "reco_beam_calibrated_dEdX")
      .Alias("reco_beam_TrkPitch_SCE", "reco_beam_TrkPitch")
      .Alias("reco_beam_calibrated_dEdX_NoSCE", "reco_beam_calibrated_dEdX_no_SCE")
      .Alias("reco_beam_TrkPitch_NoSCE", "reco_beam_TrkPitch_no_SCE")
      .Alias("reco_beam_calo_wire_NoSCE", "reco_beam_calo_wire_no_SCE")
      .Define("uncorrected_dE", dEdX_trkPitch, {"reco_beam_calibrated_dEdX_NoSCE", "reco_beam_TrkPitch_NoSCE"})
      .Define("cut_uncorrected_dEdX", cut_dEdX, {"reco_beam_calibrated_dEdX_NoSCE", "reco_beam_calo_wire_NoSCE"})
      .Define("corrected_dE", dEdX_trkPitch, {"reco_beam_calibrated_dEdX_SCE", "reco_beam_TrkPitch_SCE"})
      .Define("cut_corrected_dEdX", cut_dEdX, {"reco_beam_calibrated_dEdX_SCE", "reco_beam_calo_wire"});

   auto filtered_data5770 = frame_data5770
      .Filter("primary_isBeamType")
      .Filter("passBeamCut")
      .Filter("!isPrimaryMuonCandidate")
      .Alias("reco_beam_calibrated_dEdX_SCE", "reco_beam_calibrated_dEdX")
      .Alias("reco_beam_TrkPitch_SCE", "reco_beam_TrkPitch")
      .Alias("reco_beam_calibrated_dEdX_NoSCE", "reco_beam_calibrated_dEdX_no_SCE")
      .Alias("reco_beam_TrkPitch_NoSCE", "reco_beam_TrkPitch_no_SCE")
      .Alias("reco_beam_calo_wire_NoSCE", "reco_beam_calo_wire_no_SCE")
      .Define("uncorrected_dE", dEdX_trkPitch, {"reco_beam_calibrated_dEdX_NoSCE", "reco_beam_TrkPitch_NoSCE"})
      .Define("cut_uncorrected_dEdX", cut_dEdX, {"reco_beam_calibrated_dEdX_NoSCE", "reco_beam_calo_wire_NoSCE"})
      .Define("corrected_dE", dEdX_trkPitch, {"reco_beam_calibrated_dEdX_SCE", "reco_beam_TrkPitch_SCE"})
      .Define("cut_corrected_dEdX", cut_dEdX, {"reco_beam_calibrated_dEdX_SCE", "reco_beam_calo_wire"});

*/
   auto h2_mc_dE_wire_uncorrected = filtered_SCE
      .Histo2D({"h2_mc_dE_wire_uncorrected","MC wire vs dE not corrected for SCE", nBin_wire ,binLow_wire ,binHigh_wire ,nBin_dE, binLow_dE, binHigh_dE}, "reco_beam_calo_wire_NoSCE","uncorrected_dE");

   auto h2_mc_dE_wire_corrected = filtered_SCE
      .Histo2D({"h2_mc_dE_wire_corrected","MC wire vs dE corrected for SCE", nBin_wire ,binLow_wire ,binHigh_wire ,nBin_dE, binLow_dE, binHigh_dE}, "reco_beam_calo_wire","corrected_dE");

   auto h2_data5387_dE_wire_uncorrected = filtered_data5387
      .Histo2D({"h2_data5387_dE_wire_uncorrected","DATA 5387 wire vs dE uncorrected", nBin_wire ,binLow_wire ,binHigh_wire ,nBin_dE, binLow_dE, binHigh_dE}, "reco_beam_calo_wire_NoSCE","uncorrected_dE");

   auto h2_data5387_dE_wire_corrected = filtered_data5387
      .Histo2D({"h2_data5387_dE_wire_corrected","DATA 5387 wire vs dE corrected", nBin_wire ,binLow_wire ,binHigh_wire ,nBin_dE, binLow_dE, binHigh_dE}, "reco_beam_calo_wire","corrected_dE");

/*   
   auto h2_data5809_dE_wire_uncorrected = filtered_data5809
      .Histo2D({"h2_data5809_dE_wire_uncorrected","DATA 5809 wire vs dE uncorrected", nBin_wire ,binLow_wire ,binHigh_wire ,nBin_dE, binLow_dE, binHigh_dE}, "reco_beam_calo_wire_NoSCE","uncorrected_dE");


   auto h2_data5809_dE_wire_corrected = filtered_data5809
      .Histo2D({"h2_data5809_dE_wire_corrected","DATA 5809 wire vs dE corrected", nBin_wire ,binLow_wire ,binHigh_wire ,nBin_dE, binLow_dE, binHigh_dE}, "reco_beam_calo_wire","corrected_dE");

   auto h2_data5770_dE_wire_uncorrected = filtered_data5770
      .Histo2D({"h2_data5770_dE_wire_uncorrected","DATA 5770 wire vs dE uncorrected", nBin_wire ,binLow_wire ,binHigh_wire ,nBin_dE, binLow_dE, binHigh_dE}, "reco_beam_calo_wire_NoSCE","uncorrected_dE");

   auto h2_data5770_dE_wire_corrected = filtered_data5770
      .Histo2D({"h2_data5770_dE_wire_corrected","DATA 5770 wire vs dE corrected", nBin_wire ,binLow_wire ,binHigh_wire ,nBin_dE, binLow_dE, binHigh_dE}, "reco_beam_calo_wire","corrected_dE");
*/


   h2_mc_dE_wire_uncorrected->GetYaxis()->SetTitle("dE");
   h2_mc_dE_wire_corrected->GetYaxis()->SetTitle("dE");
   h2_data5387_dE_wire_uncorrected->GetYaxis()->SetTitle("dE");
   h2_data5387_dE_wire_corrected->GetYaxis()->SetTitle("dE");

   h2_mc_dE_wire_uncorrected->GetXaxis()->SetTitle("wire");
   h2_mc_dE_wire_corrected->GetXaxis()->SetTitle("wire");
   h2_data5387_dE_wire_uncorrected->GetXaxis()->SetTitle("wire");
   h2_data5387_dE_wire_corrected->GetXaxis()->SetTitle("wire");

   h2_mc_dE_wire_uncorrected->Write();
   h2_mc_dE_wire_corrected->Write();
   h2_data5387_dE_wire_uncorrected->Write();
   h2_data5387_dE_wire_corrected->Write();

/*
   h2_data5809_dE_wire_uncorrected->GetYaxis()->SetTitle("dE");
   h2_data5809_dE_wire_corrected->GetYaxis()->SetTitle("dE");

   h2_data5770_dE_wire_uncorrected->GetYaxis()->SetTitle("dE");
   h2_data5770_dE_wire_corrected->GetYaxis()->SetTitle("dE");

   h2_data5809_dE_wire_uncorrected->GetXaxis()->SetTitle("wire");
   h2_data5809_dE_wire_corrected->GetXaxis()->SetTitle("wire");
   h2_data5770_dE_wire_uncorrected->GetXaxis()->SetTitle("wire");
   h2_data5770_dE_wire_corrected->GetXaxis()->SetTitle("wire");


   h2_data5809_dE_wire_uncorrected->Write();
   h2_data5809_dE_wire_corrected->Write();

   h2_data5770_dE_wire_uncorrected->Write();

   h2_data5770_dE_wire_corrected->Write();
*/
   auto h2_mc_dEdX_wire_uncorrected = filtered_SCE
      .Histo2D({"h2_mc_dEdX_wire_uncorrected","MC wire vs dEdX not corrected for SCE", nBin_wire, binLow_wire, binHigh_wire, nBin_dEdX, binLow_dEdX, binHigh_dEdX}, "reco_beam_calo_wire_NoSCE","cut_uncorrected_dEdX");

   auto h2_mc_dEdX_wire_corrected = filtered_SCE
      .Histo2D({"h2_mc_dEdX_wire_corrected","MC wire vs dEdX corrected for SCE", nBin_wire, binLow_wire, binHigh_wire, nBin_dEdX, binLow_dEdX, binHigh_dEdX}, "reco_beam_calo_wire","cut_corrected_dEdX");

   auto h2_data5387_dEdX_wire_uncorrected = filtered_data5387
      .Histo2D({"h2_data5387_dEdX_wire_uncorrected","DATA 5387 wire vs dEdX uncorrected", nBin_wire, binLow_wire, binHigh_wire, nBin_dEdX, binLow_dEdX, binHigh_dEdX}, "reco_beam_calo_wire_NoSCE","cut_uncorrected_dEdX");

   auto h2_data5387_dEdX_wire_corrected = filtered_data5387
      .Histo2D({"h2_data5387_dEdX_wire_corrected","DATA 5387 wire vs dEdX corrected", nBin_wire, binLow_wire, binHigh_wire, nBin_dEdX, binLow_dEdX, binHigh_dEdX}, "reco_beam_calo_wire","cut_corrected_dEdX");

/*   
   auto h2_data5809_dEdX_wire_uncorrected = filtered_data5809
      .Histo2D({"h2_data5809_dEdX_wire_uncorrected","DATA 5809 wire vs dEdX uncorrected", nBin_wire, binLow_wire, binHigh_wire, nBin_dEdX, binLow_dEdX, binHigh_dEdX}, "reco_beam_calo_wire_NoSCE","cut_uncorrected_dEdX");

   auto h2_data5809_dEdX_wire_corrected = filtered_data5809
      .Histo2D({"h2_data5809_dEdX_wire_corrected","DATA 5809 wire vs dEdX corrected", nBin_wire, binLow_wire, binHigh_wire, nBin_dEdX, binLow_dEdX, binHigh_dEdX}, "reco_beam_calo_wire","cut_corrected_dEdX");

   auto h2_data5770_dEdX_wire_uncorrected = filtered_data5770
      .Histo2D({"h2_data5770_dEdX_wire_uncorrected","DATA 5770 wire vs dEdX uncorrected", nBin_wire, binLow_wire, binHigh_wire, nBin_dEdX, binLow_dEdX, binHigh_dEdX}, "reco_beam_calo_wire_NoSCE","cut_uncorrected_dEdX");
   auto h2_data5770_dEdX_wire_corrected = filtered_data5770
      .Histo2D({"h2_data5770_dEdX_wire_corrected","DATA 5770 wire vs dEdX corrected", nBin_wire, binLow_wire, binHigh_wire, nBin_dEdX, binLow_dEdX, binHigh_dEdX}, "reco_beam_calo_wire","cut_corrected_dEdX");
*/

   h2_mc_dEdX_wire_uncorrected->GetYaxis()->SetTitle("dEdX");
   h2_mc_dEdX_wire_corrected->GetYaxis()->SetTitle("dEdX");
   h2_data5387_dEdX_wire_uncorrected->GetYaxis()->SetTitle("dEdX");

   h2_mc_dEdX_wire_uncorrected->GetXaxis()->SetTitle("wire");
   h2_mc_dEdX_wire_corrected->GetXaxis()->SetTitle("wire");
   h2_data5387_dEdX_wire_uncorrected->GetXaxis()->SetTitle("wire");


   h2_mc_dEdX_wire_uncorrected->Write();
   h2_mc_dEdX_wire_corrected->Write();
   h2_data5387_dEdX_wire_uncorrected->Write();
   h2_data5387_dEdX_wire_corrected->Write();
/*
   h2_data5809_dEdX_wire_uncorrected->GetXaxis()->SetTitle("wire");
   h2_data5809_dEdX_wire_uncorrected->GetYaxis()->SetTitle("dEdX");
   h2_data5809_dEdX_wire_uncorrected->Write();
   h2_data5770_dEdX_wire_uncorrected->Write();

   h2_data5809_dEdX_wire_corrected->Write();
   h2_data5770_dEdX_wire_corrected->Write();
*/


   auto h2_mc_pitch_wire_uncorrected = filtered_SCE
      .Histo2D({"h2_mc_pitch_wire_uncorrected","MC wire vs pitch not corrected for SCE", nBin_wire, binLow_wire, binHigh_wire, nBin_pitch, binLow_pitch, binHigh_pitch}, "reco_beam_calo_wire_NoSCE","reco_beam_TrkPitch_NoSCE");

   auto h2_mc_pitch_wire_corrected = filtered_SCE
      .Histo2D({"h2_mc_pitch_wire_corrected","MC wire vs pitch corrected for SCE", nBin_wire, binLow_wire, binHigh_wire, nBin_pitch, binLow_pitch, binHigh_pitch}, "reco_beam_calo_wire","reco_beam_TrkPitch_SCE");

   auto h2_data5387_pitch_wire_uncorrected = filtered_data5387
      .Histo2D({"h2_data5387_pitch_wire_uncorrected","DATA 5387 wire vs pitch uncorrected", nBin_wire, binLow_wire, binHigh_wire, nBin_pitch, binLow_pitch, binHigh_pitch}, "reco_beam_calo_wire_NoSCE","reco_beam_TrkPitch_NoSCE");

   auto h2_data5387_pitch_wire_corrected = filtered_data5387
      .Histo2D({"h2_data5387_pitch_wire_corrected","DATA 5387 wire vs pitch corrected", nBin_wire, binLow_wire, binHigh_wire, nBin_pitch, binLow_pitch, binHigh_pitch}, "reco_beam_calo_wire","reco_beam_TrkPitch_SCE");

/*   
   auto h2_data5809_pitch_wire_uncorrected = filtered_data5809
      .Histo2D({"h2_data5809_pitch_wire_uncorrected","DATA 5809 wire vs pitch uncorrected", nBin_wire, binLow_wire, binHigh_wire, nBin_pitch, binLow_pitch, binHigh_pitch}, "reco_beam_calo_wire_NoSCE","reco_beam_TrkPitch_NoSCE");

   auto h2_data5809_pitch_wire_corrected = filtered_data5809
      .Histo2D({"h2_data5809_pitch_wire_corrected","DATA 5809 wire vs pitch corrected", nBin_wire, binLow_wire, binHigh_wire, nBin_pitch, binLow_pitch, binHigh_pitch}, "reco_beam_calo_wire","reco_beam_TrkPitch_SCE");

   auto h2_data5770_pitch_wire_uncorrected = filtered_data5770
      .Histo2D({"h2_data5770_pitch_wire_uncorrected","DATA 5770 wire vs pitch uncorrected", nBin_wire, binLow_wire, binHigh_wire, nBin_pitch, binLow_pitch, binHigh_pitch}, "reco_beam_calo_wire_NoSCE","reco_beam_TrkPitch_NoSCE");
   auto h2_data5770_pitch_wire_corrected = filtered_data5770
      .Histo2D({"h2_data5770_pitch_wire_corrected","DATA 5770 wire vs pitch corrected", nBin_wire, binLow_wire, binHigh_wire, nBin_pitch, binLow_pitch, binHigh_pitch}, "reco_beam_calo_wire","reco_beam_TrkPitch_SCE");
*/

   h2_mc_pitch_wire_uncorrected->GetYaxis()->SetTitle("pitch");
   h2_mc_pitch_wire_corrected->GetYaxis()->SetTitle("pitch");
   h2_data5387_pitch_wire_uncorrected->GetYaxis()->SetTitle("pitch");
   h2_data5387_pitch_wire_corrected->GetYaxis()->SetTitle("pitch");

   h2_mc_pitch_wire_uncorrected->GetXaxis()->SetTitle("wire");
   h2_mc_pitch_wire_corrected->GetXaxis()->SetTitle("wire");
   h2_data5387_pitch_wire_uncorrected->GetXaxis()->SetTitle("wire");
   h2_data5387_pitch_wire_corrected->GetXaxis()->SetTitle("wire");

   h2_mc_pitch_wire_uncorrected->Write();
   h2_mc_pitch_wire_corrected->Write();

   h2_data5387_pitch_wire_uncorrected->Write();
   h2_data5387_pitch_wire_corrected->Write();

/*   
   h2_data5809_pitch_wire_uncorrected->GetYaxis()->SetTitle("pitch");
   h2_data5809_pitch_wire_corrected->GetYaxis()->SetTitle("pitch");
   h2_data5809_pitch_wire_uncorrected->GetXaxis()->SetTitle("wire");
   h2_data5809_pitch_wire_corrected->GetXaxis()->SetTitle("wire");

   h2_data5809_pitch_wire_uncorrected->Write();
   h2_data5770_pitch_wire_uncorrected->Write();
   h2_data5809_pitch_wire_corrected->Write();
   h2_data5770_pitch_wire_corrected->Write();
*/
   auto h_data_cosThetaXZ = filtered_data5387
      .Define("data_cosThetaXZ", "reco_beam_trackDirX * beam_inst_dirX + reco_beam_trackDirY * beam_inst_dirY + reco_beam_trackDirZ * beam_inst_dirZ")
      .Histo1D({"h_data5387_cosThetaXZ", "DATA 5387 cosine Theta XZ", 200,0.5,1.2}, "data_cosThetaXZ");

   auto h_mc_cosThetaXZ = filtered_SCE
      .Define("mc_cosThetaXZ", "reco_beam_trackDirX * beam_inst_dirX + reco_beam_trackDirY * beam_inst_dirY + reco_beam_trackDirZ * beam_inst_dirZ")
      .Histo1D({"h_mc_cosThetaXZ", "mc 5387 cosine Theta XZ", 200,0.5,1.2}, "mc_cosThetaXZ");

   h_data_cosThetaXZ->GetXaxis()->SetTitle("cos");
   h_mc_cosThetaXZ->GetXaxis()->SetTitle("cos");

   h_data_cosThetaXZ->Write();
   h_mc_cosThetaXZ->Write();

   

   output->Write();
   return 0;
}
