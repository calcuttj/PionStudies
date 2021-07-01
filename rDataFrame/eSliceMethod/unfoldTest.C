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

void unfoldTest(const string mcFilepath = path)
{
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   ROOT::RDataFrame data(pionTree, dataPath);

   //no need to filter for only reconstructed events as with the Miss function one can take into account the not-reconstructed events
   auto frame = inputFrame
      .Filter("true_beam_PDG == 211");
  

   cout << "==================================== TRAIN ====================================" << endl;
  RooUnfoldResponse response (nBin_int, eEnd, eStart);

  //here I need to fill the response with the events that have true AND reco energy and note in the missed ones, the ones that are not reco'ed
  //replace the breit wigner filling

  frame
     //.Range(50000) //first 50k events
     .Foreach( [&response](double true_interE, double reco_interE){

           Double_t recoE = reco_interE;
           Double_t trueE = true_interE;

           if(reco_interE != -999.) {
            response.Fill (recoE, trueE);
            }
           else {
            response.Miss (trueE);
           };

           }
           ,{"true_interactingKE_fromLength", "reco_interactingKE"}
           );

  cout << "==================================== TEST =====================================" << endl;
  //Test on the remaining 20k events their true and reco distribution
  TH1D* hTrue= new TH1D ("true", "Test Truth",    nBin_int, eEnd, eStart);
  TH1D* hMeas= new TH1D ("meas", "Test Measured", nBin_int, eEnd, eStart);
  TH1D* hMCreco= new TH1D ("mcReco", "Test MC Reco", nBin_int, eEnd, eStart);
  
  // Test with a Gaussian, mean 500 and width 300.

  
  data
     .Filter("selected_incidentPion")
     .Foreach( [hMeas]( double reco_interE) {

           hMeas->Fill(reco_interE);

           }
           ,{"reco_interactingKE"}
           );

  frame
      .Foreach( [hTrue, hMCreco](double true_interE, double reco_interE) {

           hTrue->Fill(true_interE);
           hMCreco->Fill(reco_interE);
           
           }
           ,{"true_interactingKE_fromLength", "reco_interactingKE"}
           );

  hTrue->Scale( hMeas->Integral() / hTrue->Integral() );
  hTrue->Sumw2(0);
  hMCreco->Scale( hMeas->Integral() / hMCreco->Integral() );
  hMCreco->Sumw2(0);

  /*frame
     .Range(50001,0) // from 50k+1 until end
     .Foreach( [hTrue, hMeas](double true_interE, double reco_interE) {
           
           hMeas->Fill(reco_interE);
           hTrue->Fill(true_interE);
           
           }
           ,{"true_interactingKE_fromLength", "reco_interactingKE"}
           );
*/
  cout << "==================================== UNFOLD ===================================" << endl;
  RooUnfoldBayes   unfold (&response, hMeas, 4);    // OR
//RooUnfoldSvd     unfold (&response, hMeas, 20);   // OR
//RooUnfoldTUnfold unfold (&response, hMeas);       // OR
//RooUnfoldIds     unfold (&response, hMeas, 1);

  TH1D* hReco= (TH1D*) unfold.Hreco();

  TCanvas* c1= new TCanvas("canvas","canvas");

  //print infos with PrintTable command
  unfold.PrintTable (cout, hTrue);

  hReco->Draw();
  hMeas->Draw("SAME");
  hTrue->SetLineColor(8);
  hTrue->Draw("SAME");
  hMCreco->SetLineColor(kBlue);
  hMCreco->Draw("SAME");
  c1->BuildLegend();
  c1->SaveAs("test_unfoldMC.pdf");
   
  auto* R = response.HresponseNoOverflow();
  auto* c2 = new TCanvas();
  R->SetStats(0);
  R->Draw("colz");
  c2->Draw();


  auto cov = unfold.Ereco();
  auto* c3 = new TCanvas();
  cov.Draw("colz");
  c3->Draw();
  //c1->SaveAs("response.png");
}

#ifndef __CINT__
int main () { unfoldTest( path); return 0; }  // Main program when run stand-alone
#endif
