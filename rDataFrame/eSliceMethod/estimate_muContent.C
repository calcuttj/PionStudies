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
#include "RooUnfold.h"
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

const string path = "eSliceMethod_Prod4a_mc_1GeV_all_06_11_21.root";
const string dataPath = "eSliceMethod_Prod4a_58XX_1GeV_all_07_14_21.root";


//==============================================================================
// Estimate Mu-Excess content in Data wrt to MC
// //==============================================================================

void estimate_muContent(const string mcFilepath = path, bool doXS = true)
{
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   ROOT::RDataFrame data_inputFrame(pionTree, dataPath);

   gStyle->SetNdivisions(1020);
   string output_name;
   output_name = "out_estimate_muContent" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile( output_name.c_str() , "RECREATE"); //maybe save with binning?

   //=======================  Frame Definitions ==============================================================      
   auto frame = inputFrame
      .Define("true_equalBin", equalBin, {"true_initKE", "true_interKE"})
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"})
      //Definitions for Response Matrix generation

      .Define("truePion_response", 
            " true_beam_PDG == 211 && selected_incidentPion == 1"
            " && !true_equalBin && !reco_equalBin"
            " && reco_initKE != -999. && reco_interKE != -999.")

      .Define("truePion_miss", //either Missed Signal, or mis Reconstructed Signal
            " true_beam_PDG == 211 && !true_equalBin "
            " && (selected_incidentPion == 0 || reco_equalBin || reco_initKE == -999. || reco_interKE == -999.) ")

      .Define("truePion_fake", 
            " selected_incidentPion == 1 && (true_beam_PDG != 211 || true_equalBin)")

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
            " true_absSignal && selected_abs && true_pion_daughter == 0"
            " && !true_equalBin && !reco_equalBin"
            " && reco_initKE != -999. && reco_interKE != -999.")

      .Define("truePion_abs_miss", 
            " true_absSignal && !true_equalBin && true_pion_daughter == 0"
            " && ( reco_initKE == -999. || reco_interKE == -999. || reco_equalBin || !selected_abs )")

      .Define("truePion_abs_fake", 
            " selected_abs "
            " && ( !( true_absSignal && true_pion_daughter == 0) || true_equalBin )")
      .Define("true_initKE_fix", "return 881.;")
      .Define("true_interKE_fix", "return true_initKE_fix - (true_initKE - true_interKE);")
      .Define("reco_initKE_fix", "return 881.;")
      .Define("reco_interKE_fix", "return reco_initKE_fix - (reco_initKE - reco_interKE);");

   auto dataFrame = data_inputFrame
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"})
      .Filter("reco_initKE < 1200")
      .Define("reco_initKE_fix", "return 881.;")
      .Define("reco_interKE_fix", "return reco_initKE_fix - (reco_initKE - reco_interKE);");

   //auto count_data = dataFrame.Count();
   //std::cout << "count data = " << *count_data << std::endl;

   //=====================================================================================      

   //Frame for Measurement 
   auto frame_meas = frame 
      .Filter("true_beam_endZ > 0");
   

   //=======ADDRESS Muon Content in MC=====================================================
   //1) get muon content in selected_incidentPion
   //2) Tag a certain amount of muons --> can add them to the hMeas_mc and see how that matches hMeas_data
   //
   //LATER: they will then be also fed into Response matrix through "Fake function"
   //
   auto countAll_incidentPionSelection = frame_meas
                                        .Filter("selected_incidentPion")
                                         .Count();

   auto countMu_incidentPionSelection = frame_meas
                                        .Filter("selected_incidentPion && primaryMuon")
                                        .Count();

   double muContent = (double) *countMu_incidentPionSelection / *countAll_incidentPionSelection;

   cout << "==============================Muon Content  =====================================" << endl;
   cout << "N-ev selected incident Pion = " << *countAll_incidentPionSelection << endl;
   cout << "N-Mu selected incident Pion = " << *countMu_incidentPionSelection << endl;
   cout << "Mu Content in event Selection = " << muContent << endl;

   double muContentTarget = 0.17;
   unsigned int muExtra = (int) (*countMu_incidentPionSelection - muContentTarget * (*countAll_incidentPionSelection) ) / (muContentTarget - 1);
   cout << "Need to add " << muExtra << " to the eventSelection" << endl;

   //newMus that I need to add --> newMuContent = (alreadySelMu + newMus) / (totalSelectedPi + newMu)
   //solve for newMu
   //newMu = (alreadySelMu - muContent*totSelPi) / (muContent - 1)


   cout << "================================================================================" << endl;

   //auto count_mc = frame_meas.Count();
   //std::cout << "count MC = " << *count_mc << std::endl;
   cout << "==================================== TEST =====================================" << endl;

   //Binning
   eEnd = 0; nBin_int = 60;
   
   TH1D* hMeas_mc_initE= new TH1D ("meas_mc_initE", "Measured MC initE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_initE= new TH1D ("meas_data_initE", "Measured data initE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   TH1D* hMeas_mc_interE= new TH1D ("meas_mc_interE", "Measured MC interE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_interE= new TH1D ("meas_data_interE", "Measured data interE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   
   TH1D* hExtraMu_mc_interE= new TH1D ("hExtraMu_mc_interE", "Enhanced Muon Content to 20%; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   
   TH1D* hTrue_mc_muon= new TH1D ("hTrue_mc_muon", "MC muon Distribution; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   cout << "==================================== Building the Measured Distribution =====================================" << endl;

   frame_meas
      .Filter("selected_incidentPion")
      .Foreach( [ hMeas_mc_interE,  hExtraMu_mc_interE ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){

            hMeas_mc_interE->Fill( reco_interE );
            hExtraMu_mc_interE->Fill( reco_interE );
            }

            }
            ,{"reco_initKE_rwData", "reco_interKE_rwData", "reco_equalBin"}
            );

   //Enhancing Muon Content to Measured MC interE
   frame_meas
      .Filter("primaryMuon")
      .Range(muExtra)
      .Foreach( [ hExtraMu_mc_interE]( double reco_initE, double reco_interE, bool recoEqualBin) {

            //Extra Muons in the measured int!
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){

            hExtraMu_mc_interE->Fill( reco_interE );
            }

            }
            ,{"reco_initKE_rwData", "reco_interKE_rwData", "reco_equalBin"}
            );


   dataFrame
      .Filter("selected_incidentPion")
      .Foreach( [ hMeas_data_initE, hMeas_data_interE ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){

            hMeas_data_initE->Fill( reco_initE );
            hMeas_data_interE->Fill( reco_interE );
            }

            }
            ,{"reco_initKE", "reco_interKE", "reco_equalBin"}
            );
   
   hMeas_mc_interE->Write();
   hMeas_data_interE->Write();
   hMeas_data_initE->Write();

   cout << "==================================== Difference Data MC =====================================" << endl;

   //============= Scale MC to Data with Integral of range between 600 - eStart (bin 30 - 60)
   //
   TH1D* hMeas_mc_interE_normalise = (TH1D*)hMeas_mc_interE->Clone("hMeas_mc_interE_normalise");
   hMeas_mc_interE_normalise->SetTitle("Normalised MC to Data in Range 600 - 1200 MeV");
   hMeas_mc_interE_normalise->Scale( hMeas_data_interE->Integral(30,60) / hMeas_mc_interE->Integral(30,60) );
   
   TH1D* hDiff_data_mc_interE = (TH1D*)hMeas_data_interE->Clone("hDiff_data_mc_interE");
   hDiff_data_mc_interE->SetTitle("Difference Data - MC (Muon Excess in Data?)");
   hDiff_data_mc_interE->Add( hMeas_mc_interE_normalise, -1);

   //============= DIFFERENCE FOR ENHANCED MUON CONTENT
   //
   TH1D* hExtraMu_mc_interE_normalise = (TH1D*)hExtraMu_mc_interE->Clone("hExtraMu_mc_interE_normalise");
   hExtraMu_mc_interE_normalise->SetTitle("Normalised muon Enhanced MC to Data in Range 600 - 1200 MeV");
   hExtraMu_mc_interE_normalise->Scale( hMeas_data_interE->Integral(30,60) / hExtraMu_mc_interE->Integral(30,60) );
   
   TH1D* hDiff_enhancedMuMC = (TH1D*)hMeas_data_interE->Clone("hDiff_data_mc_interE");
   hDiff_enhancedMuMC->SetTitle("Difference Data - MC (Muon Excess in Data?)");
   hDiff_enhancedMuMC->Add( hExtraMu_mc_interE_normalise, -1);

   cout << "============= Integral of initE histo and Difference Histo " << endl;

   auto interE_integral = hMeas_data_interE->Integral();
   auto initE_integral = hMeas_data_initE->Integral();
   auto diff_integral = hDiff_data_mc_interE->Integral();

   cout << "Integral of initE histo Data = " << initE_integral << endl;
   cout << "Integral of interE histo Data = " << interE_integral << endl;
   cout << "Integral of diff MC/mu histo = " << diff_integral << endl;

   
   TCanvas* c_sel_data_mc= new TCanvas("c_sel_data_mc","c_sel_data_mc");
   gPad->SetGrid(1,1);
   hMeas_data_interE->SetFillColorAlpha(kBlack, 0.2);
   hMeas_data_interE->Draw("HIST");
   hMeas_mc_interE->SetFillColorAlpha(kBlue, 0.2);
   hMeas_mc_interE->SetLineColor(kBlue);
   hMeas_mc_interE->Draw("HIST SAME");
   c_sel_data_mc->BuildLegend();
   c_sel_data_mc->Write();
   //c_initE->SaveAs("bla.pdf");

   TCanvas* c_norm_data_mc= new TCanvas("c_norm_data_mc","c_norm_data_mc");
   c_norm_data_mc->Divide(1,2);
   c_norm_data_mc->cd(1);
   gPad->SetGrid(1,1);
   hMeas_data_interE->SetFillColorAlpha(kBlack, 0.2);
   hMeas_data_interE->Draw("HIST");
   hMeas_mc_interE_normalise->SetFillColorAlpha(kBlue, 0.2);
   hMeas_mc_interE_normalise->SetLineColor(kBlue);
   hMeas_mc_interE_normalise->Draw("HIST SAME");
   c_norm_data_mc->BuildLegend();
   c_norm_data_mc->cd(2);
   gPad->SetGrid(1,1);
   hDiff_data_mc_interE->SetLineColor(kRed);
   hDiff_data_mc_interE->SetFillColorAlpha(kRed, 0.3);
   hDiff_data_mc_interE->Draw("HIST");
   c_norm_data_mc->BuildLegend();

   c_norm_data_mc->Write();

   TCanvas* c_enhancedMu= new TCanvas("c_enhancedMu","c_enhancedMu");
   c_enhancedMu->Divide(1,2);
   c_enhancedMu->cd(1);
   gPad->SetGrid(1,1);
   hMeas_data_interE->SetFillColorAlpha(kBlack, 0.2);
   hMeas_data_interE->Draw("HIST");
   hExtraMu_mc_interE_normalise->SetFillColorAlpha(kYellow + 1, 0.2);
   hExtraMu_mc_interE_normalise->SetLineColor(kYellow + 1);
   hExtraMu_mc_interE_normalise->Draw("HIST SAME");
   c_enhancedMu->BuildLegend();
   c_enhancedMu->cd(2);
   gPad->SetGrid(1,1);
   hDiff_enhancedMuMC->SetLineColor(kRed);
   hDiff_enhancedMuMC->SetFillColorAlpha(kViolet, 0.3);
   hDiff_enhancedMuMC->Draw("HIST");
   c_enhancedMu->BuildLegend();

   c_enhancedMu->Write();


}

#ifndef __CINT__
int main () { estimate_muContent( path); return 0; }  // Main program when run stand-alone
#endif



/*
// Checking wether categorisation is not giving one event two contributions
// events not in event selection can not have a value and are just not contributing to unfolding
auto h_test = frame
.Define("sum", "truePion_response + truePion_miss + truePion_fake")
.Histo1D("sum");
h_test->Write();

frame
.Define("sum", "truePion_response + truePion_miss + truePion_fake")
.Foreach([](int sum, int pdg, bool selPi, bool recoBin, bool trueBin, double reco_initKE, double reco_interKE, bool res, bool miss, bool fake){

if(sum == 0){
std::cout << "PDG     = " << pdg << std::endl;
std::cout << "SelPi   = " << selPi << std::endl;
std::cout << "recoBin = " << recoBin << std::endl;
std::cout << "trueBin = " << trueBin << std::endl;
std::cout << "rInitKE = " << reco_initKE << std::endl;
std::cout << "rInterE = " << reco_interKE << std::endl;
std::cout << "response = " << res << std::endl;
std::cout << "miss = " << miss << std::endl;
std::cout << "fake = " << fake << std::endl;


}

},{"sum", "true_beam_PDG", "selected_incidentPion", "reco_equalBin", "true_equalBin", "reco_initKE", "reco_interKE", "truePion_response", "truePion_miss", "truePion_fake"});

exit;*/

   /*
   //==============================================================================================      
   //                         Do Ratios between MC True and Reco and apply to Data Reco
   //==============================================================================================      
  
   TH1D* hScale_data_wMCratio_totInel= new TH1D ("hScale_data_wMCratio_totInel", "Data totInel with MC true / reco; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hScale_data_wMCratio_abs= new TH1D ("hScale_data_wMCratio_abs", "Data abs with MC true / reco; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   for(int i=1; i<= nBin_int; i++){

      if(hMeas_mc_totInel->GetBinContent(i) != 0){
         double ratio_totInel = hTrue_totInel->GetBinContent(i) / hMeas_mc_totInel->GetBinContent(i);
         hScale_data_wMCratio_totInel->SetBinContent( i, hMeas_data_totInel->GetBinContent(i) * ratio_totInel );
      }
      else hScale_data_wMCratio_totInel->SetBinContent(i, hMeas_data_totInel->GetBinContent(i) );

      if(hMeas_mc_abs->GetBinContent(i) != 0){
         double ratio_abs = hTrue_abs->GetBinContent(i) / hMeas_mc_abs->GetBinContent(i);
         hScale_data_wMCratio_abs->SetBinContent(i, hMeas_data_abs->GetBinContent(i) * ratio_abs );
      }
      else hScale_data_wMCratio_abs->SetBinContent(i, hMeas_data_abs->GetBinContent(i) );


   }

   //==============================================================================================      
   //                         Do Ratios between MC True and Reco and apply to Data Reco
   //==============================================================================================      
   
   TH1D* hScale_data_wMCratio_incident= new TH1D ("hScale_data_wMCratio_incident", "Data incident with MC true / reco; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   for(int i=1; i<= nBin_int; i++){

      if(hMeas_mc_incident->GetBinContent(i) != 0){
         double ratio_incident = hTrue_incident->GetBinContent(i) / hMeas_mc_incident->GetBinContent(i);
         hScale_data_wMCratio_incident->SetBinContent( i, hMeas_data_incident->GetBinContent(i) * ratio_incident );
      }
      else hScale_data_wMCratio_incident->SetBinContent( i, hMeas_data_incident->GetBinContent(i) );
   }

   //Scale MC to Data
   hMeas_mc_incident->Scale( hMeas_data_incident->Integral() / hMeas_mc_incident->Integral() );
   hTrue_incident->Scale( hMeas_data_incident->Integral() / hTrue_incident->Integral() );

      //==============================================================================================      
      //RATIO XS
      TH1D* h_wRatioXS_totInel = new TH1D("h_wRatioXS_totInel" ,"DATA wRatio totInel XS; Energy [MeV]; #sigma [mb]", nBin_int, eEnd, eStart);

      do_XS_log(  h_wRatioXS_totInel, hScale_data_wMCratio_totInel, hScale_data_wMCratio_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_wRatioXS_totInel, hScale_data_wMCratio_totInel, hScale_data_wMCratio_incident, h_betheMean_muon );  
      //==============================================================================================      
      //RATIO XS
      TH1D* h_wRatioXS_abs = new TH1D("h_wRatioXS_abs" ,"DATA wRatio abs XS; Energy [MeV]; #sigma [mb]", nBin_int, eEnd, eStart);

      do_XS_log(  h_wRatioXS_abs, hScale_data_wMCratio_abs, hScale_data_wMCratio_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_wRatioXS_abs, hScale_data_wMCratio_abs, hScale_data_wMCratio_incident, h_betheMean_muon );  
*/
