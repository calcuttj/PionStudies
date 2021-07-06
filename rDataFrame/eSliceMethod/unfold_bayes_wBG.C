#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;
using namespace std;
using namespace ROOT::VecOps;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
//#include "RooUnfoldIds.h"
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
#include "TRandom3.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "TMatrixDBase.h"
#include "TArray.h"
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#endif

#include "../lambda.h"
#include "../betheBloch.h"
#include "eSlice.h"

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;
const string path = "eSliceMethod_Prod4a_mc_1GeV_all_06_11_21.root";
const string dataPath = "eSliceMethod_Prod4a_5387_1GeV_all_06_11_21.root";


//==============================================================================
// Example Unfolding
//==============================================================================

void unfold_bayes_wBG(const string mcFilepath = path, bool doMC = true, bool doXS = true)
{
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   //ROOT::RDataFrame data_pre(pionTree, dataPath);
   
   string output_name;
   if(doMC) output_name = "unfold_wBG_mc_" + std::to_string((int) bin_size_int) + "MeV.root";
   else output_name = "unfold_data_" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile( output_name.c_str() , "RECREATE"); //maybe save with binning?

   //no need to filter for only reconstructed events as with the Miss function one can take into account the not-reconstructed events
   auto frame = inputFrame
      .Define("true_equalBin", equalBin, {"true_initKE", "true_interKE"})
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"})
      //Definitions for Response Matrix generation
      
      //Pion InitialE and interE distributions
      .Define("truePion_response", 
            " true_beam_PDG == 211 && selected_incidentPion"
            " && !true_equalBin && !reco_equalBin"
            " && reco_initKE != -999. && reco_interKE != -999.")
      
      .Define("truePion_miss", 
            " true_beam_PDG == 211 && !true_equalBin"
             " && ( reco_initKE == -999. || reco_interKE == -999. || reco_equalBin || !selected_incidentPion )")
      
      .Define("truePion_fake", 
            " selected_incidentPion && (true_beam_PDG != 211 || true_equalBin)")

      //Pion Total Inelastic Distribution
      .Define("truePion_totInel_response", 
            " true_primPionInel && selected_incidentPion"
            " && !true_equalBin && !reco_equalBin"
            " && reco_initKE != -999. && reco_interKE != -999.")

      .Define("truePion_totInel_miss", 
            " true_primPionInel && !true_equalBin"
            " && ( reco_initKE == -999. || reco_interKE == -999. || reco_equalBin || !selected_incidentPion )")

      .Define("truePion_totInel_fake", 
            " selected_incidentPion && ( !true_primPionInel || true_equalBin ) ")

      //Pion Total Inelastic Distribution
      .Define("truePion_abs_response", 
            " true_absSignal && selected_abs"
            " && !true_equalBin && !reco_equalBin"
            " && reco_initKE != -999. && reco_interKE != -999.")

      .Define("truePion_abs_miss", 
            " true_absSignal && !true_equalBin"
            " && ( reco_initKE == -999. || reco_interKE == -999. || reco_equalBin || !selected_abs )")

      .Define("truePion_abs_fake", 
            " selected_abs && ( !true_absSignal || true_equalBin )")
      
      .Filter("true_beam_endZ > 0");
      //.Filter("true_beam_PDG == 211");

// Checking wether categorisation is not giving one event two contributions
   auto h_test = frame
      .Define("sum", "truePion_response + truePion_miss + truePion_fake")
      .Histo1D("sum");
   h_test->Write();


   cout << "==================================== TRAIN ====================================" << endl;
   //events that have reco_initBin == reco_interBin go into the Miss function as they are treated like a reco-ineff
   //events with true_initBin == true_interBin havre been filtered out already
   RooUnfoldResponse response_interE = RooUnfoldResponse(nBin_int, eEnd, eStart);
   RooUnfoldResponse response_initE = RooUnfoldResponse(nBin_int, eEnd, eStart);
   RooUnfoldResponse response_totInel = RooUnfoldResponse(nBin_int, eEnd, eStart);
   RooUnfoldResponse response_abs = RooUnfoldResponse(nBin_int, eEnd, eStart);

   frame
      //.Range(50000) //last 50k ev
      .Foreach( [ &response_initE, &response_interE](double true_initE, double true_interE, 
                                                     double reco_initE, double reco_interE, 
                                                     bool ev_response, bool ev_miss, bool ev_fake){

            //Signal
            if(ev_response){
            response_initE.Fill( reco_initE, true_initE);
            response_interE.Fill( reco_interE, true_interE);
            }
            //MisReconstruction
            else if( ev_miss ) {
            response_initE.Miss (true_initE);
            response_interE.Miss (true_interE);
            }
            //BackGround will be subtracted before unfolding
            else if( ev_fake ){

            response_initE.Fake(reco_initE);
            response_interE.Fake(reco_interE);
            }

            }
            ,{"true_initKE", "true_interKE", "reco_initKE", "reco_interKE", 
              "truePion_response", "truePion_miss", "truePion_fake"});

   //Total Inelastic
   frame
      //.Range(50000)
      //.Filter("true_primPionInel")//last 50k ev
      .Foreach( [ &response_totInel](double true_initE, double true_interE, 
                                     double reco_initE, double reco_interE, 
                                     bool ev_response, bool ev_miss, bool ev_fake){

            //Signal
            if(ev_response){
            response_totInel.Fill(reco_interE, true_interE);
            }
            //MisReconstruction
            else if( ev_miss ) {
            response_totInel.Miss (true_interE);
            }
            //BackGround will be subtracted before unfolding
            else if( ev_fake ){
            response_totInel.Fake(reco_interE);
            }

            }
            ,{"true_initKE", "true_interKE", "reco_initKE", "reco_interKE", 
              "truePion_totInel_response", "truePion_totInel_miss", "truePion_totInel_fake"});

   //Total Inelastic
   frame
      //.Range(50000)
      //.Filter("true_absSignal")//last 50k ev
      .Foreach( [ &response_abs](double true_initE, double true_interE, 
                                     double reco_initE, double reco_interE, 
                                     bool ev_response, bool ev_miss, bool ev_fake){

            //Signal
            if(ev_response){
            response_abs.Fill(reco_interE, true_interE);
            }
            //MisReconstruction
            else if( ev_miss ) {
            response_abs.Miss (true_interE);
            }
            //BackGround will be subtracted before unfolding
            else if( ev_fake ){
            response_abs.Fake(reco_interE);
            }

            }
            ,{"true_initKE", "true_interKE", "reco_initKE", "reco_interKE", 
            "truePion_abs_response", "truePion_abs_miss", "truePion_abs_fake"});


   cout << "==================================== TEST =====================================" << endl;


   TH1D* hTrue_initE= new TH1D ("trueMC_initE", "MC Truth initE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hMeas_initE= new TH1D ("meas_initE", "Measured initE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   TH1D* hTrue_interE= new TH1D ("trueMC_interE", "MC Truth interE; Energy [MeV]; ev/bin",    nBin_int, eEnd, eStart);
   TH1D* hMeas_interE= new TH1D ("meas_interE", "Measured interE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   
   TH1D* hTrue_totInel= new TH1D ("trueMC_totInel", "MC Truth totInel; Energy [MeV]; ev/bin",    nBin_int, eEnd, eStart);
   TH1D* hMeas_totInel= new TH1D ("meas_totInel", "Measured totInel; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   TH1D* hTrue_abs= new TH1D ("trueMC_abs", "MC Truth abs; Energy [MeV]; ev/bin",    nBin_int, eEnd, eStart);
   TH1D* hMeas_abs= new TH1D ("meas_abs", "Measured abs; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   //Build hMeas from same prefiltered sample that was supplied to response and fake (the Miss are misReco so won't catch them in reco)
   //Build hTrue from what was part of response and Miss (the fake are not part of hTrue
   
   cout << "==================================== Building the Measured Distribution =====================================" << endl;

   frame
      .Filter("selected_incidentPion")
      .Foreach( [ hMeas_initE, hMeas_interE, hMeas_totInel ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){

               hMeas_initE->Fill( reco_initE );
               hMeas_interE->Fill( reco_interE );
               hMeas_totInel->Fill( reco_interE );
            }

            }
            ,{"reco_initKE", "reco_interKE", "reco_equalBin"}
            );
   
   frame
      .Filter("selected_abs")
      .Foreach( [ hMeas_abs ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){
               hMeas_abs->Fill( reco_interE );
            }

            }
            ,{ "reco_initKE", "reco_interKE", "reco_equalBin"}
            );

   hMeas_initE->Write();
   hMeas_interE->Write();
   hMeas_totInel->Write();
   hMeas_abs->Write();

   cout << "==================================== Building the True Distribution =====================================" << endl;

   frame
      //.Filter("selected_incidentPion")
      .Foreach( [ hTrue_initE, hTrue_interE ]( double true_initE, double true_interE, 
                                               bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){

               hTrue_initE->Fill( true_initE );
               hTrue_interE->Fill( true_interE );
            }

            }
            ,{"true_initKE", "true_interKE", "truePion_response", "truePion_miss"}
            );
    frame
      //.Filter("selected_incidentPion")
      .Foreach( [ hTrue_totInel ]( double true_initE, double true_interE, 
                                   bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){
               hTrue_totInel->Fill( true_interE );
            }

            }
            ,{ "true_initKE", "true_interKE", "truePion_totInel_response", "truePion_totInel_miss"}
            );
  
   frame
      //.Filter("selected_abs")
      .Foreach( [ hTrue_abs ]( double true_initE, double true_interE, 
                               bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){
               hTrue_abs->Fill( true_interE );
            }

            }
            ,{ "true_initKE", "true_interKE", "truePion_abs_response", "truePion_abs_miss"}
            );

   hTrue_initE->Write();
   hTrue_interE->Write();
   hTrue_totInel->Write();
   hTrue_abs->Write();


   cout << "==================================== UNFOLD ===================================" << endl;
   RooUnfoldBayes unfold_initE = RooUnfoldBayes(&response_initE, hMeas_initE, 4);    // OR
   RooUnfoldBayes unfold_interE = RooUnfoldBayes(&response_interE, hMeas_interE, 4);    // OR
   RooUnfoldBayes unfold_totInel = RooUnfoldBayes(&response_totInel, hMeas_totInel, 4);    // OR
   RooUnfoldBayes unfold_abs = RooUnfoldBayes(&response_abs, hMeas_abs, 4);    // OR

   cout << "------------------------------------ Init E -----------------------------------" << endl;
   TH1D* hUnfold_initE = (TH1D*) unfold_initE.Hreco();
   cout << "------------------------------------ Inter E -----------------------------------" << endl;
   TH1D* hUnfold_interE = (TH1D*) unfold_interE.Hreco();
   cout << "------------------------------------ TotInel -----------------------------------" << endl;
   TH1D* hUnfold_totInel = (TH1D*) unfold_totInel.Hreco();
   cout << "------------------------------------ abs -----------------------------------" << endl;
   TH1D* hUnfold_abs = (TH1D*) unfold_abs.Hreco();

   hUnfold_initE->SetNameTitle("hUnfold_initE", "Unfold initE; Energy [MeV]; ev/bin");
   hUnfold_interE->SetNameTitle("hUnfold_interE", "Unfold interE; Energy [MeV]; ev/bin");
   hUnfold_totInel->SetNameTitle("hUnfold_totInel", "Unfold totInel; Energy [MeV]; ev/bin");
   hUnfold_abs->SetNameTitle("hUnfold_abs", "Unfold abs; Energy [MeV]; ev/bin");

   hUnfold_initE->Write();
   hUnfold_interE->Write();
   hUnfold_totInel->Write();
   hUnfold_abs->Write();
   //print infos with PrintTable command
   cout << "==================================== UNFOLD INIT E===================================" << endl;
   unfold_initE.PrintTable (cout, hTrue_initE);
   cout << "==================================== UNFOLD INTER E===================================" << endl;
   unfold_interE.PrintTable (cout, hTrue_interE);
   cout << "==================================== UNFOLD TOTINEL E===================================" << endl;
   unfold_totInel.PrintTable (cout, hTrue_totInel);
   cout << "==================================== UNFOLD abs E===================================" << endl;
   unfold_abs.PrintTable (cout, hTrue_abs);
   
   cout << "==================================== CHECK INTEGRALS===================================" << endl;
   cout << "Integral of Measured initE = " << hMeas_initE->Integral() << endl;
   cout << "Integral of Measured interE = " << hMeas_interE->Integral() << endl;
   cout << "Integral of True initE = " << hTrue_initE->Integral() << endl;
   cout << "Integral of True interE = " << hTrue_interE->Integral() << endl;
   cout << "Integral of Unfold initE = " << hUnfold_initE->Integral() << endl;
   cout << "Integral of Unfold interE = " << hUnfold_interE->Integral() << endl;

   //InitE has not same Integral as InterE
   //Scale initialE to number of interacting Pions for incident, bc interE for XS also is in the nominator  
   hUnfold_initE->Scale( hUnfold_interE->Integral() / hUnfold_initE->Integral() );

   TCanvas* c_initE= new TCanvas("canvas_initE","canvas_initE");
   hUnfold_initE->Draw();
   hMeas_initE->SetFillColorAlpha(kBlack,0.1);
   hMeas_initE->Draw("SAME");
   hTrue_initE->SetLineColor(8);
   hTrue_initE->Draw("SAME");
   c_initE->BuildLegend();
   c_initE->Write();
   //c_initE->SaveAs("bla.pdf");

   TCanvas* c_interE= new TCanvas("canvas_interE","canvas_interE");
   hUnfold_interE->Draw();
   hMeas_interE->SetFillColorAlpha(kBlack,0.1);
   hMeas_interE->Draw("SAME");
   hTrue_interE->SetLineColor(8);
   hTrue_interE->Draw("SAME");
   c_interE->BuildLegend();
   c_interE->Write();

   TCanvas* c_totInel= new TCanvas("canvas_totInel","canvas_totInel");
   hUnfold_totInel->Draw();
   hMeas_totInel->SetFillColorAlpha(kBlack,0.1);
   hMeas_totInel->Draw("SAME");
   hTrue_totInel->SetLineColor(8);
   hTrue_totInel->Draw("SAME");
   c_totInel->BuildLegend();
   c_totInel->Write();

   TCanvas* c_abs= new TCanvas("canvas_abs","canvas_abs");
   hUnfold_abs->Draw();
   hMeas_abs->SetFillColorAlpha(kBlack,0.1);
   hMeas_abs->Draw("SAME");
   hTrue_abs->SetLineColor(8);
   hTrue_abs->Draw("SAME");
   c_abs->BuildLegend();
   c_abs->Write();

   //Build the Incident histos in order to compare

   TH1D* hTrue_incident= new TH1D ("trueMC_incident", "MC Truth incident",    nBin_int, eEnd, eStart);
   TH1D* hMeas_incident= new TH1D ("meas_incident", "Measured incident", nBin_int, eEnd, eStart);
   TH1D* hUnfold_incident= new TH1D ("hUnfold_incident", "Unfolded incident", nBin_int, eEnd, eStart);

   build_incidentHist(hTrue_initE, hTrue_interE, hTrue_incident);
   build_incidentHist(hMeas_initE, hMeas_interE, hMeas_incident);
   build_incidentHist(hUnfold_initE, hUnfold_interE, hUnfold_incident);

   TCanvas* c_incident= new TCanvas("canvas_incident","canvas_incident");
   hUnfold_incident->Draw();
   hMeas_incident->SetFillColorAlpha(kBlack,0.1);
   hMeas_incident->Draw("SAME HIST");
   hTrue_incident->SetLineColor(8);
   hTrue_incident->Draw("SAME HIST");
   //h_mcReco_incident->SetLineColor(38);
   //h_mcReco_incident->Draw("SAME HIST");
   c_incident->BuildLegend();
   c_incident->Write();

   if(doXS){
   //access Jakes GeantFile in folder
   TFile f1("exclusive_xsec.root");
   TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
   TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
   TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
   f1.Close();

   output->cd();

   TH1D* h_betheMean_muon = new TH1D("h_betheMean_muon", "Mean Energy Loss", nBin_int, eEnd, eStart);

   //fill histo with Mean dEdX of bin center
   for(int i = 1; i <= nBin_int; i++){
      h_betheMean_muon->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size_int  , mass_muon) );
   };


   TH1D* h_unfoldXS_totInel = new TH1D("h_unfoldXS_totInel_data" ,"Data Unfold XS total Inelastic; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
      
   do_XS_log(  h_unfoldXS_totInel, hUnfold_totInel, hUnfold_incident, h_betheMean_muon );   
   do_XS_log_binomial_error( h_unfoldXS_totInel, hUnfold_totInel, hUnfold_incident, h_betheMean_muon );   

   TCanvas *c_unfold_totInel = new TCanvas("c_unfold_totInel", "c_unfold_totInel");
   gPad->SetGrid(1,1);
   h_unfoldXS_totInel->SetTitle( "Unfolded Total Inelastic XS; Energy [MeV]; #sigma [mb]");
   h_unfoldXS_totInel->GetXaxis()->SetRangeUser(200,1000);
   h_unfoldXS_totInel->GetXaxis()->SetNdivisions(1020);
   h_unfoldXS_totInel->GetYaxis()->SetNdivisions(1020);

   totInel_KE->SetTitle( "Unfolded Total Inelasitc XS");
   totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   totInel_KE->SetLineColor(kRed);
   totInel_KE->SetLineWidth(3);
   totInel_KE->Draw("AC");
   h_unfoldXS_totInel->SetMarkerSize(0.7);
   h_unfoldXS_totInel->Draw("PE0 SAME");

   c_unfold_totInel->Write();

   TH1D* h_unfoldXS_abs = new TH1D("h_unfoldXS_abs_data" ,"Data Unfold XS total Inelastic; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
      
   do_XS_log(  h_unfoldXS_abs, hUnfold_abs, hUnfold_incident, h_betheMean_muon );   
   do_XS_log_binomial_error( h_unfoldXS_abs, hUnfold_abs, hUnfold_incident, h_betheMean_muon );   

   TCanvas *c_unfold_abs = new TCanvas("c_unfold_abs", "c_unfold_abs");
   gPad->SetGrid(1,1);
   h_unfoldXS_abs->SetTitle( "Unfolded Total Inelastic XS; Energy [MeV]; #sigma [mb]");
   h_unfoldXS_abs->GetXaxis()->SetRangeUser(200,1000);
   h_unfoldXS_abs->GetXaxis()->SetNdivisions(1020);
   h_unfoldXS_abs->GetYaxis()->SetNdivisions(1020);

   abs_KE->SetTitle( "Unfolded Total Inelasitc XS");
   abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
   abs_KE->SetLineColor(kBlue);
   abs_KE->SetLineWidth(3);
   abs_KE->Draw("AC");
   h_unfoldXS_abs->SetMarkerSize(0.7);
   h_unfoldXS_abs->Draw("PE0 SAME");

   c_unfold_abs->Write();

   }



}

#ifndef __CINT__
int main () { unfold_bayes_wBG( path); return 0; }  // Main program when run stand-alone
#endif
