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
#include "betheBloch.h"
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>

int fit_dEdX_correction_try(const string file_path) {
   // TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
   // c1->SetGridx();
   // c1->SetGridy();
   // c1->GetFrame()->SetFillColor(21);
   // c1->GetFrame()->SetBorderMode(-1);
   // c1->GetFrame()->SetBorderSize(5);

   //string to const char*
   const char *c_file_path = file_path.c_str();
   TFile *file =new TFile( c_file_path , "READ");
   if (!file) {
      std::cout << "Couldn't open File, wrong path? " << std::endl;
      return 1;
   }

   //get TH2D histos from file
   TH2D* h2_5387_pitch_wire = (TH2D*)file->Get("h2_data5387_pitch_wire_uncorrected");
   TH2D* h2_5387_dEdX_wire = (TH2D*)file->Get("h2_data5387_dEdX_wire_uncorrected");
   
   TH2D* h2_5387_SCEcorr_pitch_wire = (TH2D*)file->Get("h2_data5387_pitch_wire_corrected");
   TH2D* h2_5387_SCEcorr_dEdX_wire = (TH2D*)file->Get("h2_data5387_dEdX_wire_corrected");
   
   //file->GetObject("h2_data5387_pitch_wire_uncorrected", h2_5387_pitch_wire);
   if (!h2_5387_pitch_wire || !h2_5387_dEdX_wire || !h2_5387_SCEcorr_dEdX_wire || !h2_5387_SCEcorr_pitch_wire){
      std::cout << "Couldnt find object " << h2_5387_pitch_wire << std::endl;
      return 2;
   }

   TFile *output = new TFile("output_fit_5387.root", "RECREATE");
   //h2_5387_pitch_wire->DrawClone("COLZ");
   //h2_5387_pitch_wire->Write();

   double temp_mean = 0.;
   double temp_std = 0.;
   double temp_maxEntry = 0.;
   int temp_maxLeft = 0.;
   int temp_maxRight = 0.;
   int minEntries = 100;

   int nbinWire = h2_5387_pitch_wire->GetNbinsX();

   //histos of the Fit
   //--- somehow are saved to file after exiting the loop, no need to save again
   TH1D *fit_pitch_mean = new TH1D("fit_pitch_mean", "", h2_5387_pitch_wire->GetNbinsX(), 1, h2_5387_pitch_wire->GetNbinsX() );
   TH1D *fit_pitch_std = new TH1D("fit_pitch_std", "", h2_5387_pitch_wire->GetNbinsX(), 1, h2_5387_pitch_wire->GetNbinsX() );
   TH1D *fit_pitch_chi2 = new TH1D("fit_pitch_chi2", "", h2_5387_pitch_wire->GetNbinsX(), 1, h2_5387_pitch_wire->GetNbinsX() );
   TH1D *fit_pitch_ndf = new TH1D("fit_pitch_ndf", "", h2_5387_pitch_wire->GetNbinsX(), 1, h2_5387_pitch_wire->GetNbinsX() );

   TH1D *fit_pitch_SCEcorr_mean = new TH1D("fit_pitch_SCEcorr_mean", "", h2_5387_SCEcorr_pitch_wire->GetNbinsX(), 1, h2_5387_SCEcorr_pitch_wire->GetNbinsX() );
   TH1D *fit_pitch_SCEcorr_std = new TH1D("fit_pitch_SCEcorr_std", "", h2_5387_SCEcorr_pitch_wire->GetNbinsX(), 1, h2_5387_SCEcorr_pitch_wire->GetNbinsX() );
   TH1D *fit_pitch_SCEcorr_chi2 = new TH1D("fit_pitch_SCEcorr_chi2", "", h2_5387_SCEcorr_pitch_wire->GetNbinsX(), 1, h2_5387_SCEcorr_pitch_wire->GetNbinsX() );
   TH1D *fit_pitch_SCEcorr_ndf = new TH1D("fit_pitch_SCEcorr_ndf", "", h2_5387_SCEcorr_pitch_wire->GetNbinsX(), 1, h2_5387_SCEcorr_pitch_wire->GetNbinsX() );

   //Loop through 2D histo get projection of each bin and perform Gaus Fit

   for(int i = 1; i <= h2_5387_pitch_wire->GetNbinsX(); i++){
      //for(int i = 383; i <= 390; i++){

      TH1D *temp_projection = h2_5387_pitch_wire->ProjectionY("py",i,i);
  
      //fit only if Entries > 100
      if( temp_projection->GetEntries() > minEntries){
         
         //bin with highest entry WATCHOUT hard coded Bins..
         int max_bin = temp_projection->GetMaximumBin();
      
         temp_maxEntry = temp_projection->GetBinCenter( max_bin );
         temp_std = temp_projection->GetStdDev(1);
         TF1* f1 = new TF1("f1", "gaus", temp_maxEntry - temp_std , temp_maxEntry + temp_std);
         //TF1* f1 = new TF1("f1", "gaus", temp_mean - temp_std/2 , temp_mean + temp_std/2);
         temp_projection->Fit("f1", "RS QN");
         //gPad->WaitPrimitive();
         //save Parameters of Gaus Fit each to a histo
         fit_pitch_mean->SetBinContent( i , f1->GetParameter(1) );
         fit_pitch_mean->SetBinError( i , f1->GetParError(1) );
         fit_pitch_std->SetBinContent( i , f1->GetParameter(2) );
         fit_pitch_std->SetBinError( i , f1->GetParError(2) );
         fit_pitch_chi2->SetBinContent( i , f1->GetChisquare() );
         fit_pitch_ndf->SetBinContent( i , f1->GetNDF() );

         //Delete fitted function and histo
         delete f1;
         delete temp_projection;
      }
   };

   for(int i = 1; i <= h2_5387_SCEcorr_pitch_wire->GetNbinsX(); i++){
      //for(int i = 65; i <= 70; i++){

      TH1D *temp_projection = h2_5387_SCEcorr_pitch_wire->ProjectionY("py",i,i);

      //fit only if Entries > 100
      if( temp_projection->GetEntries() > minEntries){

         int max_bin = temp_projection->GetMaximumBin();
      
         temp_maxEntry = temp_projection->GetBinCenter( max_bin );
         temp_std = temp_projection->GetStdDev(1);
         TF1* f1 = new TF1("f1", "gaus", temp_maxEntry - temp_std , temp_maxEntry + temp_std);
         //TF1* f1 = new TF1("f1", "gaus", temp_mean - temp_std/2 , temp_mean + temp_std/2);
         temp_projection->Fit("f1", "RS QN");
          //gPad->WaitPrimitive();
         //save Parameters of Gaus Fit each to a histo
         fit_pitch_SCEcorr_mean->SetBinContent( i , f1->GetParameter(1) );
         fit_pitch_SCEcorr_mean->SetBinError( i , f1->GetParError(1) );
         fit_pitch_SCEcorr_std->SetBinContent( i , f1->GetParameter(2) );
         fit_pitch_SCEcorr_std->SetBinError( i , f1->GetParError(2) );
         fit_pitch_SCEcorr_chi2->SetBinContent( i , f1->GetChisquare() );
         fit_pitch_SCEcorr_ndf->SetBinContent( i , f1->GetNDF() );

         //Delete fitted function and histo
         delete f1;
         delete temp_projection;
      }
   };


   
   TH1D *fit_dEdX_mpv = new TH1D("fit_dEdX_mpv", "", h2_5387_dEdX_wire->GetNbinsX(), 1, h2_5387_dEdX_wire->GetNbinsX() );
   TH1D *fit_dEdX_std = new TH1D("fit_dEdX_std", "", h2_5387_dEdX_wire->GetNbinsX(), 1, h2_5387_dEdX_wire->GetNbinsX() );
   TH1D *fit_dEdX_chi2 = new TH1D("fit_dEdX_chi2", "", h2_5387_dEdX_wire->GetNbinsX(), 1, h2_5387_dEdX_wire->GetNbinsX() );
   TH1D *fit_dEdX_ndf = new TH1D("fit_dEdX_ndf", "", h2_5387_dEdX_wire->GetNbinsX(), 1, h2_5387_dEdX_wire->GetNbinsX() );

   TH1D *fit_dEdX_SCEcorr_mpv = new TH1D("fit_dEdX_SCEcorr_mpv", "", h2_5387_SCEcorr_dEdX_wire->GetNbinsX(), 1, h2_5387_SCEcorr_dEdX_wire->GetNbinsX() );
   TH1D *fit_dEdX_SCEcorr_std = new TH1D("fit_dEdX_SCEcorr_std", "", h2_5387_SCEcorr_dEdX_wire->GetNbinsX(), 1, h2_5387_SCEcorr_dEdX_wire->GetNbinsX() );
   TH1D *fit_dEdX_SCEcorr_chi2 = new TH1D("fit_dEdX_SCEcorr_chi2", "", h2_5387_SCEcorr_dEdX_wire->GetNbinsX(), 1, h2_5387_SCEcorr_dEdX_wire->GetNbinsX() );
   TH1D *fit_dEdX_SCEcorr_ndf = new TH1D("fit_dEdX_SCEcorr_ndf", "", h2_5387_SCEcorr_dEdX_wire->GetNbinsX(), 1, h2_5387_SCEcorr_dEdX_wire->GetNbinsX() );

   //Loop through 2D histo get projection of each bin and perform Gaus Fit
   for(int i = 1; i <= h2_5387_dEdX_wire->GetNbinsX(); i++){
   //   for(int i = 100; i <= 120; i++){

      TH1D *temp_projection = h2_5387_dEdX_wire->ProjectionY("py",i,i);
      //temp_projection->Draw();
      //gPad->WaitPrimitive();
      //temp_projection->Close();

      //fit only if Entries > 100
      if( temp_projection->GetEntries() > minEntries){

         //bin with highest entry WATCHOUT hard coded Bins..
         int max_bin = temp_projection->GetMaximumBin();
         int widthLeft = max_bin - 100;
         int widthRight = max_bin + 100;

         temp_maxEntry = temp_projection->GetBinContent( max_bin );
         widthLeft = temp_projection->FindFirstBinAbove( 0.2*temp_maxEntry, 1, widthLeft, max_bin );
         widthRight = temp_projection->FindLastBinAbove( 0.2*temp_maxEntry, 1, max_bin , widthRight );

         TF1* f2 = new TF1("f2", "landau", temp_projection->GetBinCenter(widthLeft) , temp_projection->GetBinCenter(widthRight));
         //temp_projection->Fit("f2", "RS ");
         temp_projection->Fit("f2", "RS QN");
         //gPad->WaitPrimitive();
         //save Parameters of Landau Fit each to a histo
         //Landa Par1 MPV is not real MPV, do GetMaximumX for landau most probable value
         fit_dEdX_mpv->SetBinContent( i , f2->GetMaximumX() );
         fit_dEdX_mpv->SetBinError( i , f2->GetParError(1) );
         fit_dEdX_std->SetBinContent( i , f2->GetParameter(2) );
         fit_dEdX_std->SetBinError( i , f2->GetParError(2) );
         fit_dEdX_chi2->SetBinContent( i , f2->GetChisquare() );
         fit_dEdX_ndf->SetBinContent( i , f2->GetNDF() );

         //Delete fitted function and histo
         delete f2;
         delete temp_projection;
      }
   };

   for(int i = 1; i <= h2_5387_SCEcorr_dEdX_wire->GetNbinsX(); i++){
   //   for(int i = 77; i <= 87; i++){

      TH1D *temp_projection = h2_5387_SCEcorr_dEdX_wire->ProjectionY("py",i,i);
      //temp_projection->Draw();
      //gPad->WaitPrimitive();
      //temp_projection->Close();

      //fit only if Entries > 100
      if( temp_projection->GetEntries() > minEntries){

         //bin with highest entry WATCHOUT hard coded Bins..
         int max_bin = temp_projection->GetMaximumBin();
         int widthLeft = max_bin - 100;
         int widthRight = max_bin + 100;

         temp_maxEntry = temp_projection->GetBinContent( max_bin );
         widthLeft = temp_projection->FindFirstBinAbove( 0.2*temp_maxEntry, 1, widthLeft, max_bin );
         widthRight = temp_projection->FindLastBinAbove( 0.2*temp_maxEntry, 1, max_bin , widthRight );

         TF1* f2 = new TF1("f2", "landau", temp_projection->GetBinCenter(widthLeft) , temp_projection->GetBinCenter(widthRight));
         temp_projection->Fit("f2", "RS QN");
         //gPad->WaitPrimitive();
         //save Parameters of Landau Fit each to a histo
         //Landa Par1 MPV is not real MPV, do GetMaximumX for landau most probable value
         fit_dEdX_SCEcorr_mpv->SetBinContent( i , f2->GetMaximumX() );
         fit_dEdX_SCEcorr_mpv->SetBinError( i , f2->GetParError(1) );
         fit_dEdX_SCEcorr_std->SetBinContent( i , f2->GetParameter(2) );
         fit_dEdX_SCEcorr_std->SetBinError( i , f2->GetParError(2) );
         fit_dEdX_SCEcorr_chi2->SetBinContent( i , f2->GetChisquare() );
         fit_dEdX_SCEcorr_ndf->SetBinContent( i , f2->GetNDF() );

         //Delete fitted function and histo
         delete f2;
         delete temp_projection;
      }
   };

   //prepare for correction get pitch_true = pitch_reco * dEdX_reco (mpv) / dEdX_true (mpv)
   
   /* std::cout << "Mean Energy Loss Bethe Bloch minEntriesMeV = " << betheBloch(1000) << std::endl;
   std::cout << "MPV Bethe bloch minEntriesMeV = " << betheBloch_mpv(1000) << std::endl;
   std::cout << "Mean Energy Loss Bethe Bloch 860MeV = " << betheBloch(860) << std::endl;
   std::cout << "MPV Bethe bloch 860MeV = " << betheBloch_mpv(860) << std::endl;
   std::cout << "Mean Energy Loss Bethe Bloch 500MeV = " << betheBloch(500) << std::endl;
   std::cout << "MPV Bethe bloch 500MeV = " << betheBloch_mpv(500) << std::endl;
   */

   TH1D* hist_betheBloch_mpv = new TH1D("betheBloch_mpv", "", h2_5387_pitch_wire->GetNbinsX(), 1, h2_5387_pitch_wire->GetNbinsX() );

   //initial beam energy
   double Etemp = 890;
   int cnt =0;
   for(int i=1; i <= fit_pitch_mean->GetNbinsX(); i++){
      auto temp = fit_pitch_mean->GetBinContent(i);
      //should be from wire 68 on
      //first wire with hit
      if(temp > 0) {
         hist_betheBloch_mpv->SetBinContent( i, 2*betheBloch_mpv(Etemp)); 
         hist_betheBloch_mpv->SetBinError(i, 0.001 );
         Etemp = Etemp - betheBloch(Etemp)*0.5; //should put pitch apparent at that point
         //energy at each passage is reduced by mean value of bethe bloch
      }
      //if(temp > 0) {
      //   hist_betheBloch_mpv->SetBinContent( i, 2*betheBloch_mpv(860 - cnt)); 
      //   hist_betheBloch_mpv->SetBinError(i, 0.001 );
      //   cnt++; //energy at each passage is reduced by 1 MeV / bc of wire pitch
      //}
   };

   hist_betheBloch_mpv->Write();


   TH1D* pitch_true_mean = new TH1D("true_pitch_mean", "", h2_5387_pitch_wire->GetNbinsX(), 1, h2_5387_pitch_wire->GetNbinsX() );
   //Test scale
   //fit_pitch_mean->Scale(1.05);
   //fit_dEdX_mpv->Scale(0.95);

   pitch_true_mean->Multiply(fit_pitch_mean,fit_dEdX_mpv);
   pitch_true_mean->Divide(hist_betheBloch_mpv);
   //pitch_true_mean->Write();

   //TH1D* hypothesis_pitch_mean = new TH1D("hypothesis_pitch_mean", "", h2_5387_pitch_wire->GetNbinsX(), 1, h2_5387_pitch_wire->GetNbinsX() );

   //hypothesis_pitch_mean->Multiply(fit_pitch_SCEcorr_mean, fit_dEdX_mpv);
   //hypothesis_pitch_mean->Divide(hist_betheBloch_mpv);
   //hypothesis_pitch_mean->Write();

   //Graphs comparing pitch sum pandoracalo and uncorrected
   double sum_pitch_pandoracalo[nbinWire], sum_pitch_true[nbinWire], wire[nbinWire];
   double pandoracalo_err[nbinWire], true_err[nbinWire];
   double diff_true_pandoracalo[nbinWire];

   double sum1 = 0 , sum2 = 0;
   for(int i=1; i <= h2_5387_pitch_wire->GetNbinsX(); i++){
   //for(int i=1; i <= 80; i++){
      
      wire[i-1] = i;
      //pandoracalo
      sum1 = sum1 + fit_pitch_mean->GetBinContent(i);
      //std::cout << "Pandoracalo pitch " << sum1 << std::endl;
      sum_pitch_pandoracalo[i-1] = sum1;
      pandoracalo_err[i-1] = fit_pitch_mean->GetBinError(i);
      //std::cout << "Pandoracalo filled = " << sum_pitch_pandoracalo[i] << std::endl;

      //effective
      sum2 = sum2 + pitch_true_mean->GetBinContent(i);
      sum_pitch_true[i-1] = sum2 ;
      true_err[i-1] = pitch_true_mean->GetBinError(i);

         //diffs
      if(sum1 ==0) {
         diff_true_pandoracalo[i-1] = 0;
      }
      else{
         diff_true_pandoracalo[i-1] = sum2 - sum1;
      }
   };

   TMultiGraph *mg1 = new TMultiGraph();
   TMultiGraph *mg2 = new TMultiGraph();
   
   TGraphErrors* gr_pitch_pandoracalo_true = new TGraphErrors( nbinWire, sum_pitch_pandoracalo, sum_pitch_true, pandoracalo_err, true_err );

   TGraph* gr_diff_pandoracalo_true = new TGraph( nbinWire, sum_pitch_pandoracalo, diff_true_pandoracalo);
   gr_pitch_pandoracalo_true->SetLineColor(kBlue);
   gr_pitch_pandoracalo_true->SetTitle("L real");
   
   gr_diff_pandoracalo_true->SetLineColor(kBlue);
   gr_diff_pandoracalo_true->SetTitle("Difference L real - L apparent");
    
   TCanvas *c_pitch_pandoracalo_true = new TCanvas("c_pitch_pandoracalo_true", "", 1600, 1800);
   gr_pitch_pandoracalo_true->Draw("ALP");
   c_pitch_pandoracalo_true->Write();
   c_pitch_pandoracalo_true->Close();

   mg1->Add(gr_pitch_pandoracalo_true);

   mg2->Add(gr_diff_pandoracalo_true);

   TCanvas *c_combined = new TCanvas("graphs", "", 1600, 1800);
   c_combined->Divide(2,1);
   c_combined->cd(1);
   mg1->Draw("ALP");
   c_combined->cd(2);
   mg2->Draw("ALP");
   c_combined->Write();
   c_combined->Close();

   h2_5387_pitch_wire->Write();
   h2_5387_dEdX_wire->Write();
   h2_5387_SCEcorr_pitch_wire->Write();
   h2_5387_SCEcorr_dEdX_wire->Write();
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

   return 0;
   }

//   TH1D* pitch_sum = new TH1D("sum_true_pitch_wire", "", h2_5387_pitch_wire->GetNbinsX(), 1*0.5, h2_5387_pitch_wire->GetNbinsX()*0.5 );
//   
//   double sum = 0;
//   for(int i = 1; i <= pitch_sum->GetNbinsX(); i++){
//
//      //help = help + energy_true_mean->GetBinContent(i);
//      //energy_sum->SetBinContent(i, help);
//      sum = sum + pitch_true_mean->GetBinContent(i);
//      //pitch_sum->SetBinContent(i , sum);
//      pitch_sum->SetBinContent(i , sum / (i*0.5) );
//
//   };
//
//   //energy_sum->Write();
//   //
//   //TEST multiply by corrected pitch measured dEdX
//
//   TH1D* test_pitch_sum = new TH1D("test_sum_true_pitch_wire", "", h2_5387_pitch_wire->GetNbinsX(), 1*0.5, h2_5387_pitch_wire->GetNbinsX()*0.5 );
//   
//   sum = 0;
//   for(int i = 1; i <= test_pitch_sum->GetNbinsX(); i++){
//
//      sum = sum + test_multiply->GetBinContent(i);
//      //test_pitch_sum->SetBinContent(i , sum);
//      test_pitch_sum->SetBinContent(i , sum /(i*0.5));
//
//   };
//
//
//   test_pitch_sum->Write();
//
//   pitch_sum->Write();
//
