#include "TCanvas.h"
#include "TFrame.h"
#include "TBenchmark.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TROOT.h"
#include "TPaveLabel.h"
#include "TError.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "lambda.h"
#include "fitfunc.h"
#include "betheBloch.h"
#include <ROOT/RDataFrame.hxx>
#include "eSliceMethod/eSlice.h"


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>

void hist_bethe(double E_init, double mass_particle, TH1D* fit_mean, TH1D* h_bethe ){
   //initial beam energy
   int cnt =0;
   for(int i=1; i <= fit_mean->GetNbinsX(); i++){
      auto temp = fit_mean->GetBinContent(i);
      //should be from wire 68 on
      //first wire with hit
      //if(temp > 0) {
         h_bethe->SetBinContent( i, 2*betheBloch_mpv(E_init, mass_particle)); 
         h_bethe->SetBinError(i, 0.001 );
         E_init = E_init - betheBloch(E_init, mass_particle)*0.5; //should put pitch apparent at that point
         //energy at each passage is reduced by mean value of bethe bloch
      //}
      if(E_init <= 0) return;
   };
};

void bethe_mean(double E_init, double mass_particle, TH1D* fit_mean, TH1D* h_bethe ){
   //initial beam energy
   int cnt =0;
   for(int i=1; i <= fit_mean->GetNbinsX(); i++){
      auto temp = fit_mean->GetBinContent(i);
      //should be from wire 68 on
      //first wire with hit
      //if(temp > 0) {
         h_bethe->SetBinContent( i, betheBloch(E_init, mass_particle)); 
         h_bethe->SetBinError(i, 0.001 );
         E_init = E_init - betheBloch(E_init, mass_particle)*0.5; //should put pitch apparent at that point
         
         //energy at each passage is reduced by mean value of bethe bloch
      //}
      
      if(E_init <= 0) return;

   };
};

int fit_dEdX_mc_validation(const string file_path) {

   //string to const char*
   const char *c_file_path = file_path.c_str();
   TFile *file =new TFile( c_file_path , "READ");
   if (!file) {
      std::cout << "Couldn't open File, wrong path? " << std::endl;
      return 1;
   }

   //get TH2D histos from file
   TH2D* h2_mc_pitch_wire = (TH2D*)file->Get("h2_mc_pitch_wire_uncorrected");
   TH2D* h2_mc_dEdX_wire = (TH2D*)file->Get("h2_mc_dEdX_wire_uncorrected");

   TH2D* h2_mc_SCEcorr_pitch_wire = (TH2D*)file->Get("h2_mc_pitch_wire_corrected");
   TH2D* h2_mc_SCEcorr_dEdX_wire = (TH2D*)file->Get("h2_mc_dEdX_wire_corrected");

   if (!h2_mc_pitch_wire || !h2_mc_dEdX_wire || !h2_mc_SCEcorr_dEdX_wire || !h2_mc_SCEcorr_pitch_wire){
      std::cout << "Couldnt find object " << h2_mc_pitch_wire << std::endl;
      return 2;
   }

   TFile *output = new TFile("output_fit_mc_Prod4.root", "RECREATE");
   int minEntries = 100;

   int nbinWire = h2_mc_pitch_wire->GetNbinsX();

   //histos of the Fit
   TH1D *fit_pitch_mean = new TH1D("fit_mc_pitch_mean", "", h2_mc_pitch_wire->GetNbinsX(), 1, h2_mc_pitch_wire->GetNbinsX() );
   fit_pitch_mean->GetXaxis()->SetTitle("wire");
   fit_pitch_mean->GetYaxis()->SetTitle("pitch (cm)");
   TH1D *fit_pitch_std = new TH1D("fit_mc_pitch_std", "", h2_mc_pitch_wire->GetNbinsX(), 1, h2_mc_pitch_wire->GetNbinsX() );
   TH1D *fit_pitch_chi2 = new TH1D("fit_mc_pitch_chi2", "", h2_mc_pitch_wire->GetNbinsX(), 1, h2_mc_pitch_wire->GetNbinsX() );
   TH1D *fit_pitch_ndf = new TH1D("fit_mc_pitch_ndf", "", h2_mc_pitch_wire->GetNbinsX(), 1, h2_mc_pitch_wire->GetNbinsX() );

   TH1D *fit_pitch_SCEcorr_mean = new TH1D("fit_mc_pitch_SCEcorr_mean", "", h2_mc_SCEcorr_pitch_wire->GetNbinsX(), 1, h2_mc_SCEcorr_pitch_wire->GetNbinsX() );
   fit_pitch_SCEcorr_mean->GetXaxis()->SetTitle("wire");
   fit_pitch_SCEcorr_mean->GetYaxis()->SetTitle("pitch (cm)");
   TH1D *fit_pitch_SCEcorr_std = new TH1D("fit_mc_pitch_SCEcorr_std", "", h2_mc_SCEcorr_pitch_wire->GetNbinsX(), 1, h2_mc_SCEcorr_pitch_wire->GetNbinsX() );
   TH1D *fit_pitch_SCEcorr_chi2 = new TH1D("fit_mc_pitch_SCEcorr_chi2", "", h2_mc_SCEcorr_pitch_wire->GetNbinsX(), 1, h2_mc_SCEcorr_pitch_wire->GetNbinsX() );
   TH1D *fit_pitch_SCEcorr_ndf = new TH1D("fit_mc_pitch_SCEcorr_ndf", "", h2_mc_SCEcorr_pitch_wire->GetNbinsX(), 1, h2_mc_SCEcorr_pitch_wire->GetNbinsX() );

   //Call fit function for Gaus fit
   gaus_fit( minEntries, h2_mc_pitch_wire, fit_pitch_mean, fit_pitch_std, fit_pitch_chi2, fit_pitch_ndf);
   gaus_fit( minEntries, h2_mc_SCEcorr_pitch_wire, fit_pitch_SCEcorr_mean, fit_pitch_SCEcorr_std, fit_pitch_SCEcorr_chi2, fit_pitch_SCEcorr_ndf);

   TH1D *fit_dEdX_mpv = new TH1D("fit_mc_dEdX_mpv", "", h2_mc_dEdX_wire->GetNbinsX(), 1, h2_mc_dEdX_wire->GetNbinsX() );
   fit_dEdX_mpv->GetXaxis()->SetTitle("wire");
   fit_dEdX_mpv->GetYaxis()->SetTitle("dEdX (MeV/cm)"); 
  TH1D *fit_dEdX_mean = new TH1D("fit_dEdX_mean", "", h2_mc_dEdX_wire->GetNbinsX(), 1, h2_mc_dEdX_wire->GetNbinsX() );
   fit_dEdX_mean->GetXaxis()->SetTitle("wire");
   fit_dEdX_mean->GetYaxis()->SetTitle("dEdX (MeV/cm)");    TH1D *fit_dEdX_std = new TH1D("fit_mc_dEdX_std", "", h2_mc_dEdX_wire->GetNbinsX(), 1, h2_mc_dEdX_wire->GetNbinsX() );
   TH1D *fit_dEdX_chi2 = new TH1D("fit_mc_dEdX_chi2", "", h2_mc_dEdX_wire->GetNbinsX(), 1, h2_mc_dEdX_wire->GetNbinsX() );
   TH1D *fit_dEdX_ndf = new TH1D("fit_mc_dEdX_ndf", "", h2_mc_dEdX_wire->GetNbinsX(), 1, h2_mc_dEdX_wire->GetNbinsX() );

   TH1D *fit_dEdX_SCEcorr_mpv = new TH1D("fit_mc_dEdX_SCEcorr_mpv", "", h2_mc_SCEcorr_dEdX_wire->GetNbinsX(), 1, h2_mc_SCEcorr_dEdX_wire->GetNbinsX() );
   fit_dEdX_SCEcorr_mpv->GetXaxis()->SetTitle("wire");
   fit_dEdX_SCEcorr_mpv->GetYaxis()->SetTitle("dEdX (MeV/cm)"); 
   TH1D *fit_dEdX_SCEcorr_mean = new TH1D("fit_dEdX_SCEcorr_mean", "", h2_mc_SCEcorr_dEdX_wire->GetNbinsX(), 1, h2_mc_SCEcorr_dEdX_wire->GetNbinsX() );
   fit_dEdX_SCEcorr_mean->GetXaxis()->SetTitle("wire");
   TH1D *fit_dEdX_SCEcorr_std = new TH1D("fit_mc_dEdX_SCEcorr_std", "", h2_mc_SCEcorr_dEdX_wire->GetNbinsX(), 1, h2_mc_SCEcorr_dEdX_wire->GetNbinsX() );
   TH1D *fit_dEdX_SCEcorr_chi2 = new TH1D("fit_mc_dEdX_SCEcorr_chi2", "", h2_mc_SCEcorr_dEdX_wire->GetNbinsX(), 1, h2_mc_SCEcorr_dEdX_wire->GetNbinsX() );
   TH1D *fit_dEdX_SCEcorr_ndf = new TH1D("fit_mc_dEdX_SCEcorr_ndf", "", h2_mc_SCEcorr_dEdX_wire->GetNbinsX(), 1, h2_mc_SCEcorr_dEdX_wire->GetNbinsX() );

   //Call fit function for Landau
   landau_fit( minEntries, h2_mc_dEdX_wire, fit_dEdX_mpv, fit_dEdX_mean, fit_dEdX_std, fit_dEdX_chi2, fit_dEdX_ndf );
   landau_fit( minEntries, h2_mc_SCEcorr_dEdX_wire, fit_dEdX_SCEcorr_mpv, fit_dEdX_mean, fit_dEdX_SCEcorr_std, fit_dEdX_SCEcorr_chi2, fit_dEdX_SCEcorr_ndf );

   //Lifetime correction histogram --> This has to be applied to the non corrected pandoracalonosce object
   //MC has 35ms
   //Build a correction histo to multiply dEdx with. exp(tDrift(wire) / lifetime)

   double x_drift = 360, x_beam_start = 28, v_drift = 0.1648, lifetime = 35000; //cm and us lifetime is 35ms
   //first bin with content in dEdX bin (x axis)
   Int_t first_wire = fit_dEdX_mpv->FindFirstBinAbove(0,1, 1,80);
   //std::cout << "first wire = " << first_wire << std::endl;

   TH1D* h_lifetime_correction = new TH1D("lifetimecorrection", "", h2_mc_dEdX_wire->GetNbinsX(), 1, h2_mc_dEdX_wire->GetNbinsX());

   //fill hist with correction factors
   for(int i=1; i <= h_lifetime_correction->GetNbinsX(); i++){
      //for(int i=1; i <= 100; i++){

      if(i >= first_wire){

         double dist_z = (i - first_wire )*0.48;
         //std::cout << "Dist z = " << dist_z << std::endl;
         double t_drift = ( x_drift - x_beam_start - dist_z*TMath::Tan( TMath::ACos(0.98) ))/ v_drift ;
         //std::cout << "T drift = " << t_drift << std::endl;
         double correction_factor = 1 / exp( - t_drift / lifetime ); //inverse bc charge gets attenuated more due to lover lifetime in detector so have to correct for it
         //std::cout << "correction factor = " << correction_factor << std::endl;
         h_lifetime_correction->SetBinContent(i, correction_factor);
         h_lifetime_correction->SetBinError(i, 0.01);
      }


   };
   h_lifetime_correction->Write();

   //correct uncorrected pandoracalonosce for lifetime
   TH1D* h_dEdX_mpv_lifetime = new TH1D("dEdX_mpv_lifetime", "", h2_mc_dEdX_wire->GetNbinsX(), 1, h2_mc_dEdX_wire->GetNbinsX());
   h_dEdX_mpv_lifetime->Multiply(h_lifetime_correction, fit_dEdX_mpv);
   h_dEdX_mpv_lifetime->Write();

   TH1D* h_dEdX_mean_lifetime = new TH1D("dEdX_mean_lifetime", "", h2_mc_dEdX_wire->GetNbinsX(), 1, h2_mc_dEdX_wire->GetNbinsX());
   h_dEdX_mean_lifetime->Multiply(h_lifetime_correction, fit_dEdX_mean);
   h_dEdX_mean_lifetime->Write();

   //Bethe Bloch
   TH1D* h_betheMPV_pion = new TH1D("betheMPV_pion", "", h2_mc_pitch_wire->GetNbinsX(), 1, h2_mc_pitch_wire->GetNbinsX() );
   h_betheMPV_pion->GetXaxis()->SetTitle("wire");
   h_betheMPV_pion->GetYaxis()->SetTitle("dEdX (MeV/cm)");
   TH1D* h_betheMPV_muon = new TH1D("betheMPV_muon", "", h2_mc_pitch_wire->GetNbinsX(), 1, h2_mc_pitch_wire->GetNbinsX() );
   h_betheMPV_muon->GetXaxis()->SetTitle("wire");
   h_betheMPV_muon->GetYaxis()->SetTitle("dEdX (MeV/cm)");

   //double mass_pion = 139, mass_muon = 105.6;

   //maybe get energy of pions and muons from distribution of MC? better
   double E_in_pion =  KE_in_pion, E_in_muon = KE_in_muon ;

   hist_bethe( E_in_pion, mass_pion, fit_pitch_mean, h_betheMPV_pion);
   hist_bethe( E_in_muon, mass_muon, fit_pitch_mean, h_betheMPV_muon);

   h_betheMPV_pion->Write();
   h_betheMPV_muon->Write();

   TH1D* h_betheMean_pion = new TH1D("betheMean_pion", "", h2_mc_pitch_wire->GetNbinsX(), 1, h2_mc_pitch_wire->GetNbinsX() );
   h_betheMean_pion->GetXaxis()->SetTitle("wire");
   h_betheMean_pion->GetYaxis()->SetTitle("dEdX (MeV/cm)");
   TH1D* h_betheMean_muon = new TH1D("betheMean_muon", "", h2_mc_pitch_wire->GetNbinsX(), 1, h2_mc_pitch_wire->GetNbinsX() );
   h_betheMean_muon->GetXaxis()->SetTitle("wire");
   h_betheMean_muon->GetYaxis()->SetTitle("dEdX (MeV/cm)");
   
   E_in_pion =  KE_in_pion, E_in_muon = KE_in_muon;
   
   bethe_mean( E_in_pion, mass_pion, fit_pitch_mean, h_betheMean_pion);
   bethe_mean( E_in_muon, mass_muon, fit_pitch_mean, h_betheMean_muon);

   h_betheMean_pion->Write();
   h_betheMean_muon->Write();


   TH1D* pitch_mean_old = new TH1D("true_pitch_mean_old", "", h2_mc_pitch_wire->GetNbinsX(), 1, h2_mc_pitch_wire->GetNbinsX() );

   // use for calculation of pitch true mean the uncalibrated but corrected dEdx
   pitch_mean_old->Multiply(fit_pitch_mean,fit_dEdX_mpv);
   pitch_mean_old->Divide(h_betheMPV_muon);
   pitch_mean_old->Write();

   TH1D* pitch_calc_lifetime = new TH1D("true_pitch_mean_lifetime", "", h2_mc_pitch_wire->GetNbinsX(), 1, h2_mc_pitch_wire->GetNbinsX() );

   // use for calculation of pitch true mean the uncalibrated but corrected dEdx
   pitch_calc_lifetime->Multiply(fit_pitch_mean,h_dEdX_mpv_lifetime);
   pitch_calc_lifetime->Divide(h_betheMPV_muon);
   pitch_calc_lifetime->Write();

   double mc_cosThetaXZ = 0.98;

   //Graphs comparing pitch sum pandoracalonosce and uncorrected
   double sum_pitch_pandoracalonosce[nbinWire], sum_pitch_true[nbinWire], sum_pitch_SCEcorr[nbinWire], wire[nbinWire];
   double pandoracalonosce_err[nbinWire], true_err[nbinWire], SCEcorr_err[nbinWire];
   double diff_true_pandoracalonosce[nbinWire], diff_pandoracalonosce_SCEcorr[nbinWire];

   double sum1 = 0 , sum2 = 0, sum3 = 0;   
   for(int i=1; i <= h2_mc_pitch_wire->GetNbinsX(); i++){
      //for(int i=1; i <= 80; i++){

      wire[i-1] = i;
      //pandoracalonosce
      sum1 = sum1 + fit_pitch_mean->GetBinContent(i) * mc_cosThetaXZ ;
      //std::cout << "pandoracalonosce pitch " << sum1 << std::endl;
      sum_pitch_pandoracalonosce[i-1] = sum1;
      pandoracalonosce_err[i-1] = fit_pitch_mean->GetBinError(i);
      //std::cout << "pandoracalonosce filled = " << sum_pitch_pandoracalonosce[i] << std::endl;

      //effective
      sum2 = sum2 + pitch_calc_lifetime->GetBinContent(i) * mc_cosThetaXZ;
      sum_pitch_true[i-1] = sum2 ;
      true_err[i-1] = pitch_calc_lifetime->GetBinError(i);

      //SCE corr
      sum3 = sum3 + fit_pitch_SCEcorr_mean->GetBinContent(i) * mc_cosThetaXZ;
      sum_pitch_SCEcorr[i-1] = sum3 ;
      SCEcorr_err[i-1] = fit_pitch_SCEcorr_mean->GetBinError(i);

      //diff pandoracalonosce true
      if(sum1 ==0) {
         diff_true_pandoracalonosce[i-1] = 0;
      }
      else{
         diff_true_pandoracalonosce[i-1] = sum2 - sum1;
      }
      //diff SCE corr true
      if(sum3 ==0) {
         diff_pandoracalonosce_SCEcorr[i-1] = 0;
      }
      else{
         diff_pandoracalonosce_SCEcorr[i-1] = sum3 - sum1;
      }
   };


   TMultiGraph *mg1 = new TMultiGraph();
   TMultiGraph *mg2 = new TMultiGraph();

   TGraphErrors* gr_pitch_pandoracalonosce_true = new TGraphErrors( nbinWire, sum_pitch_pandoracalonosce, sum_pitch_true, pandoracalonosce_err, true_err );

   TGraph* gr_diff_pandoracalonosce_true = new TGraph( nbinWire, wire, diff_true_pandoracalonosce);
   gr_pitch_pandoracalonosce_true->SetLineColor(kBlue);
   gr_pitch_pandoracalonosce_true->SetTitle("L calculated");
   gr_pitch_pandoracalonosce_true->GetXaxis()->SetTitle("L uncalibrated");
   gr_pitch_pandoracalonosce_true->GetYaxis()->SetTitle("L calculated");

   gr_diff_pandoracalonosce_true->SetLineColor(kBlue);
   gr_diff_pandoracalonosce_true->SetTitle("Difference L calculated - L uncalibrated");
   gr_diff_pandoracalonosce_true->GetXaxis()->SetTitle("L uncalibrated");
   gr_diff_pandoracalonosce_true->GetYaxis()->SetTitle("Difference");

   TCanvas *c_true_pandoracalonosce = new TCanvas("Lapp_vs_Lreal", "", 1600, 1800);
   c_true_pandoracalonosce->SetGrid();
   c_true_pandoracalonosce->Divide(2,1);
   c_true_pandoracalonosce->cd(1);
   gr_pitch_pandoracalonosce_true->Draw("ALP");
   c_true_pandoracalonosce->cd(2);
   gr_diff_pandoracalonosce_true->Draw("ALP");
   c_true_pandoracalonosce->Write();
   c_true_pandoracalonosce->Close();

   //true vs SCEcorr

   TGraphErrors* gr_pitch_SCEcorr_pandoracalonosce = new TGraphErrors( nbinWire, sum_pitch_SCEcorr, sum_pitch_pandoracalonosce, SCEcorr_err, pandoracalonosce_err );

   TGraph* gr_diff_SCEcorr_pandoracalonosce = new TGraph( nbinWire, wire, diff_pandoracalonosce_SCEcorr);
   gr_pitch_SCEcorr_pandoracalonosce->SetLineColor(kBlue);
   gr_pitch_SCEcorr_pandoracalonosce->SetTitle("L uncalibrated vs L SCEcorr");
   gr_pitch_SCEcorr_pandoracalonosce->GetXaxis()->SetTitle("L SCE corrected");
   gr_pitch_SCEcorr_pandoracalonosce->GetYaxis()->SetTitle("L uncalibrated");

   gr_diff_SCEcorr_pandoracalonosce->SetLineColor(kBlue);
   gr_diff_SCEcorr_pandoracalonosce->SetTitle("Difference L uncalibrated - L SCEcorr");
   gr_diff_SCEcorr_pandoracalonosce->GetXaxis()->SetTitle("L uncalibrated");
   gr_diff_SCEcorr_pandoracalonosce->GetYaxis()->SetTitle("Difference");
 
   TCanvas *c_pandoracalonosce_SCEcorr = new TCanvas("L_SCEcorr_vs_Lapp", "", 1600, 1800);
   c_pandoracalonosce_SCEcorr->SetGrid();
   c_pandoracalonosce_SCEcorr->Divide(2,1);
   c_pandoracalonosce_SCEcorr->cd(1);
   gr_pitch_SCEcorr_pandoracalonosce->Draw("ALP");
   c_pandoracalonosce_SCEcorr->cd(2);
   gr_diff_SCEcorr_pandoracalonosce->Draw("ALP");
   c_pandoracalonosce_SCEcorr->Write();
   c_pandoracalonosce_SCEcorr->Close();

   h2_mc_pitch_wire->Write();
   h2_mc_dEdX_wire->Write();
   h2_mc_SCEcorr_pitch_wire->Write();
   h2_mc_SCEcorr_dEdX_wire->Write();
   fit_pitch_mean->Write();
   fit_pitch_std->Write();
   fit_pitch_chi2->Write();
   fit_pitch_ndf->Write();
   fit_pitch_SCEcorr_mean->Write();
   fit_pitch_SCEcorr_std->Write();
   fit_pitch_SCEcorr_chi2->Write();
   fit_pitch_SCEcorr_ndf->Write();

   fit_dEdX_mpv->Write();
   fit_dEdX_std->Write();
   fit_dEdX_chi2->Write();
   fit_dEdX_ndf->Write();

   fit_dEdX_SCEcorr_mpv->Write();
   fit_dEdX_SCEcorr_std->Write();
   fit_dEdX_SCEcorr_chi2->Write();
   fit_dEdX_SCEcorr_ndf->Write();
   file->Close();

   TCanvas *canv_dEdx_mc = new TCanvas("canv_dEdx_mc", "canv_dEdx_mc",260,102,1241,777);
   //canv_dEdx_mc->Range(-177.2384,0.8946241,947.5077,3.175632);
   canv_dEdx_mc->SetFillColor(0);
   canv_dEdx_mc->SetBorderMode(0);
   canv_dEdx_mc->SetBorderSize(1);
   canv_dEdx_mc->SetGridx();
   canv_dEdx_mc->SetGridy();
   canv_dEdx_mc->SetTickx(1);
   canv_dEdx_mc->SetTicky(1);
   canv_dEdx_mc->SetLeftMargin(0.1584699);
   canv_dEdx_mc->SetRightMargin(0.1311475);
   canv_dEdx_mc->SetBottomMargin(0.1509147);
   canv_dEdx_mc->SetFrameBorderMode(0);
   canv_dEdx_mc->SetFrameBorderSize(0);
   canv_dEdx_mc->SetFrameBorderMode(0);
   canv_dEdx_mc->SetFrameBorderSize(0);

   fit_dEdX_mpv->SetMinimum(0.35);
   fit_dEdX_mpv->SetMaximum(2.9);
   fit_dEdX_mpv->SetLineColor(8);
   fit_dEdX_mpv->SetLineWidth(2);
   fit_dEdX_mpv->SetMarkerColor(8);
   fit_dEdX_mpv->SetMarkerStyle(20);
   fit_dEdX_mpv->SetMarkerSize(0.8);
   fit_dEdX_mpv->GetXaxis()->SetTitle("wire");
   fit_dEdX_mpv->GetXaxis()->SetNdivisions(1020);
   fit_dEdX_mpv->GetYaxis()->SetTitle("dEdX (MeV/cm)");
   fit_dEdX_mpv->GetYaxis()->SetNdivisions(1020);
   fit_dEdX_mpv->Draw("");

   h_dEdX_mpv_lifetime->SetMinimum(0.35);
   h_dEdX_mpv_lifetime->SetMaximum(2.9);
   h_dEdX_mpv_lifetime->SetLineColor(9);
   h_dEdX_mpv_lifetime->SetLineWidth(2);
   h_dEdX_mpv_lifetime->SetMarkerColor(9);
   h_dEdX_mpv_lifetime->SetMarkerStyle(20);
   h_dEdX_mpv_lifetime->SetMarkerSize(0.8);
   h_dEdX_mpv_lifetime->Draw("same");

   fit_dEdX_SCEcorr_mpv->SetMinimum(0.35);
   fit_dEdX_SCEcorr_mpv->SetMaximum(2.9);
   fit_dEdX_SCEcorr_mpv->SetLineColorAlpha(6, 0.6);
   fit_dEdX_SCEcorr_mpv->SetLineWidth(2);
   fit_dEdX_SCEcorr_mpv->SetMarkerColorAlpha(6, 0.6);
   fit_dEdX_SCEcorr_mpv->SetMarkerStyle(20);
   fit_dEdX_SCEcorr_mpv->SetMarkerSize(0.8);
   fit_dEdX_SCEcorr_mpv->Draw("same");

   h_betheMPV_pion->SetMinimum(0.35);
   h_betheMPV_pion->SetMaximum(2.9);
   h_betheMPV_pion->SetLineColor(44);
   h_betheMPV_pion->SetLineWidth(3);
   h_betheMPV_pion->SetMarkerColor(44);
   h_betheMPV_pion->SetMarkerStyle(3);
   h_betheMPV_pion->SetMarkerSize(1);
   h_betheMPV_pion->Draw("same");

   h_betheMPV_muon->SetMinimum(0.35);
   h_betheMPV_muon->SetMaximum(2.9);
   h_betheMPV_muon->SetLineColor(49);
   h_betheMPV_muon->SetLineWidth(3);
   h_betheMPV_muon->SetMarkerColor(49);
   h_betheMPV_muon->SetMarkerStyle(3);
   h_betheMPV_muon->SetMarkerSize(1);
   h_betheMPV_muon->Draw("same");
   
   //canv_dEdx_mc->BuildLegend();
   auto legend = new TLegend(0.25,0.76,0.6, 0.9,NULL,"brNDC");
   legend->AddEntry(fit_dEdX_mpv, "pandoracalonosce");
   legend->AddEntry(h_dEdX_mpv_lifetime, "pandoracalonosce + lifetime");
   legend->AddEntry(fit_dEdX_SCEcorr_mpv, "pandoracalinoxyzt");
   legend->AddEntry(h_betheMPV_pion, "betheBloch MPV Pion");
   legend->AddEntry(h_betheMPV_muon, "betheBloch MPV Muon 900 -50 MeV ekin");
   legend->Draw();
   
   TPaveLabel *pl = new TPaveLabel(187.8666,3.010961,565.6985,3.142088,"MC dEdX","br");
   pl->Draw();
   
   canv_dEdx_mc->Modified();
   canv_dEdx_mc->Update();
   canv_dEdx_mc->cd();
   canv_dEdx_mc->SetSelected(canv_dEdx_mc);
   canv_dEdx_mc->ToggleToolBar();
   canv_dEdx_mc->Write();
   //canv_dEdx_mc->Close();

   TCanvas *canv_pitch_data = new TCanvas("canv_pitch_data", "canv_pitch_data",260,102,1241,777);
   //canv_pitch_data->Range(-177.2384,0.8946241,947.5077,3.175632);
   canv_pitch_data->SetFillColor(0);
   canv_pitch_data->SetBorderMode(0);
   canv_pitch_data->SetBorderSize(1);
   canv_pitch_data->SetGridx();
   canv_pitch_data->SetGridy();
   canv_pitch_data->SetTickx(1);
   canv_pitch_data->SetTicky(1);
   canv_pitch_data->SetLeftMargin(0.1584699);
   canv_pitch_data->SetRightMargin(0.1311475);
   canv_pitch_data->SetBottomMargin(0.1509147);
   canv_pitch_data->SetFrameBorderMode(0);
   canv_pitch_data->SetFrameBorderSize(0);
   canv_pitch_data->SetFrameBorderMode(0);
   canv_pitch_data->SetFrameBorderSize(0);

   fit_pitch_mean->SetMinimum(0.35);
   fit_pitch_mean->SetMaximum(0.8);
   fit_pitch_mean->SetLineColor(8);
   fit_pitch_mean->SetLineWidth(2);
   fit_pitch_mean->SetMarkerColor(8);
   fit_pitch_mean->SetMarkerStyle(20);
   fit_pitch_mean->SetMarkerSize(0.8);
   fit_pitch_mean->GetXaxis()->SetTitle("wire");
   fit_pitch_mean->GetXaxis()->SetNdivisions(1020);
   fit_pitch_mean->GetYaxis()->SetTitle("pitch (cm)");
   fit_pitch_mean->GetYaxis()->SetNdivisions(1020);
   fit_pitch_mean->Draw("");

   pitch_calc_lifetime->SetMinimum(0.35);
   pitch_calc_lifetime->SetMaximum(0.8);
   pitch_calc_lifetime->SetLineColor(9);
   pitch_calc_lifetime->SetLineWidth(2);
   pitch_calc_lifetime->SetMarkerColor(9);
   pitch_calc_lifetime->SetMarkerStyle(20);
   pitch_calc_lifetime->SetMarkerSize(0.8);
   pitch_calc_lifetime->Draw("same");

   fit_pitch_SCEcorr_mean->SetMinimum(0.35);
   fit_pitch_SCEcorr_mean->SetMaximum(0.8);
   fit_pitch_SCEcorr_mean->SetLineColorAlpha(6, 0.6);
   fit_pitch_SCEcorr_mean->SetLineWidth(2);
   fit_pitch_SCEcorr_mean->SetMarkerColorAlpha(6, 0.6);
   fit_pitch_SCEcorr_mean->SetMarkerStyle(20);
   fit_pitch_SCEcorr_mean->SetMarkerSize(0.8);
   fit_pitch_SCEcorr_mean->Draw("same");

   
   //canv_pitch_data->BuildLegend();
   auto legend_pitch = new TLegend(0.25,0.76,0.6, 0.9,NULL,"brNDC");
   legend_pitch->AddEntry(fit_pitch_mean, "pandoracalonosce");
   legend_pitch->AddEntry(pitch_calc_lifetime ,"calculated pitch");
   legend_pitch->AddEntry(fit_pitch_SCEcorr_mean, "pandoracalinoxyzt");
   legend_pitch->Draw();
   
   TPaveLabel *pla = new TPaveLabel(187.8666,3.010961,565.6985,3.142088,"Data mc pitch","br");
   pla->Draw();
   
   canv_pitch_data->Modified();
   canv_pitch_data->Update();
   canv_pitch_data->cd();
   canv_pitch_data->SetSelected(canv_pitch_data);
   canv_pitch_data->ToggleToolBar();
   canv_pitch_data->Write();
   //canv_pitch_data->Close();

   output->Write();
   //output->Close();

   //Save for eSlice Method values
   TFile *output_fitval = new TFile("output_fit_mc_Prod4_eSlice.root", "RECREATE");
   fit_pitch_SCEcorr_mean->Write();
   fit_pitch_SCEcorr_std->Write();
   fit_pitch_SCEcorr_chi2->Write();
   fit_pitch_ndf->Write();
   fit_pitch_mean->Write();
   fit_pitch_std->Write();
   fit_pitch_chi2->Write();
   fit_pitch_ndf->Write();
   fit_dEdX_mpv->Write();
   fit_dEdX_std->Write();
   fit_dEdX_chi2->Write();
   fit_dEdX_ndf->Write();
   h_dEdX_mpv_lifetime->Write();
   fit_dEdX_SCEcorr_mpv->Write();
   fit_dEdX_SCEcorr_std->Write();
   fit_dEdX_SCEcorr_chi2->Write();
   fit_dEdX_SCEcorr_ndf->Write();

   output_fitval->Write();
   output_fitval->Close();


   

   return 0;
   }


