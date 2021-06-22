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
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"
#include "TMatrixDBase.h"
#include "TArray.h"
#include "../lambda.h"
#include "../betheBloch.h"
#include "eSlice.h"
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include <vector>

//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;

//Applying the eSlice Method

//***********************
//Main Function
// Test fitting gaussians to the smearing matrices

int smearing_gaus_test(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");

   TFile *output = new TFile ("output_smear_gaus.root", "RECREATE");

 
   auto frame = inputFrame      
      .Define("true_initKE", "true_firstEntryIncident")
      .Define("true_interKE", "true_interactingKE_fromLength")
      .Filter("true_beam_endZ > 0");

   //Build the True Process and TrueE Int and Inc Histograms that we need to compare unsmeared things to
   //
   //all available after beamCuts
   auto mc_preReco_allPions = frame.Filter("true_beam_PDG == 211");
   auto eventSel_post_pandoraReco = frame.Filter("primary_isBeamType");


   TH2D* h2_smearing_abs_int = new TH2D("h2_smearing_abs_int", "Smearing Absorption Interacting; reco InteractingKE; true InteractingKE", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_smearing_totInel_int = new TH2D("h2_smearing_totInel_int", "Smearing TotInel Interacting; reco InteractingKE; true InteractingKE", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_smearing_incident_initE = new TH2D("h2_smearing_incident_initE", "Smearing TruePion initialE; reco IncidentKE; true IncidentKE", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   TH2D* h2_smearing_incident_interE = new TH2D("h2_smearing_incident_interE", "Smearing TruePion interE; reco InteractingKE; true InteractingKE", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   //Fill smearing matrix with entry 1 on the diagonal to avoid them being "uninvertible"
   for(int i = 1; i <= nBin_int; i++){
      h2_smearing_incident_initE->SetBinContent(i,i,1);
      h2_smearing_incident_interE->SetBinContent(i,i,1);
      h2_smearing_abs_int->SetBinContent(i,i,1);
      h2_smearing_totInel_int->SetBinContent(i,i,1);
   };

    eventSel_post_pandoraReco 
      .Filter("true_beam_PDG == 211")
      .Foreach( [h2_smearing_incident_initE, h2_smearing_incident_interE] (double true_initE, double true_interE, double reco_initE, double reco_interE ){

            fill_smearing_matrix(h2_smearing_incident_initE, true_initE, true_interE, reco_initE, reco_interE, false);
            
            fill_smearing_matrix(h2_smearing_incident_interE, true_initE, true_interE, reco_initE, reco_interE, true);

            }
            ,{"true_initKE", "true_interKE", "reco_firstEntryIncident", "reco_interactingKE"}); 
  
   //normalise_smearing(h2_smearing_incident_initE);
   //normalise_smearing(h2_smearing_incident_interE);

   h2_smearing_incident_initE->Sumw2(0);
   h2_smearing_incident_initE->Write();

   h2_smearing_incident_interE->Sumw2(0);
   h2_smearing_incident_interE->Write();

   eventSel_post_pandoraReco
      .Filter("true_absSignal")
      //.Filter("true_absSignal && true_reco_initE_eq_interE")
      .Foreach( [h2_smearing_abs_int] (double true_initE, double true_interE, double reco_initE, double reco_interE ){

               fill_smearing_matrix( h2_smearing_abs_int, true_initE, true_interE, reco_initE, reco_interE, true );
            }
            ,{"true_initKE", "true_interKE", "reco_firstEntryIncident", "reco_interactingKE"}); 

   //normalise_smearing(h2_smearing_abs_int);

   h2_smearing_abs_int->Sumw2(0);
   h2_smearing_abs_int->Write();

   eventSel_post_pandoraReco
      .Filter("true_primPionInel")
      .Foreach( [h2_smearing_totInel_int] (double true_initE, double true_interE, double reco_initE, double reco_interE ){

               fill_smearing_matrix( h2_smearing_totInel_int, true_initE, true_interE, reco_initE, reco_interE, true );
            }
            ,{"true_initKE", "true_interKE", "reco_firstEntryIncident", "reco_interactingKE"}); 

   //normalise_smearing(h2_smearing_totInel_int);

   h2_smearing_totInel_int->Sumw2(0);
   h2_smearing_totInel_int->Write();


   //Test fitting with FitSlicesX(0,1,nBin_int)
   //0 --> gaus func
   //all bins of true_KE 1 to nBins, ignore overflow and underflow
   //
   h2_smearing_incident_initE->FitSlicesX(0,1,nBin_int, 0, "QNR");
   h2_smearing_incident_interE->FitSlicesX(0,1,nBin_int, 0, "QNR");
   h2_smearing_abs_int->FitSlicesX(0,1,nBin_int, 0, "QNR");
   h2_smearing_totInel_int->FitSlicesX(0,1,nBin_int, 0, "QNR");

//----------------------------------------------------------------
   TCanvas *c1_incident_initE = new TCanvas("c1_incident_initE", "");
   c1_incident_initE->Divide(1,2);
   c1_incident_initE->cd(1);
   c1_incident_initE->cd(1)->SetGrid();
   //Get Mean
   TH2D *h2_smearing_incident_initE_1 = (TH2D*)output->Get("h2_smearing_incident_initE_1");
   h2_smearing_incident_initE_1->SetTitle("True Pion initial Energy, fit Mean;reco incidentKE;");  
   h2_smearing_incident_initE_1->GetYaxis()->SetRangeUser(0,1000);
   h2_smearing_incident_initE_1->Draw();

   c1_incident_initE->cd(2);
   c1_incident_initE->cd(2)->SetGrid();
   //Get Mean
   TH2D *h2_smearing_incident_initE_2 = (TH2D*)output->Get("h2_smearing_incident_initE_2");
   h2_smearing_incident_initE_2->SetTitle("True Pion initial Energy, fit Sigma;reco incidentKE;");  
   h2_smearing_incident_initE_2->GetYaxis()->SetRangeUser(0,100);
   h2_smearing_incident_initE_2->Draw();

   TCanvas *c2_incident_initE = new TCanvas("c2_incident_initE", "");
   c2_incident_initE->SetGrid();
   h2_smearing_incident_initE->Draw("COLZ");

//----------------------------------------------------------------
   TCanvas *c1_incident_interE = new TCanvas("c1_incident_interE", "");
   c1_incident_interE->Divide(1,2);
   c1_incident_interE->cd(1);
   c1_incident_interE->cd(1)->SetGrid();
   //Get Mean
   TH2D *h2_smearing_incident_interE_1 = (TH2D*)output->Get("h2_smearing_incident_interE_1");
   h2_smearing_incident_interE_1->SetTitle("True Pion interacting Energy, fit Mean;reco interactingKE;");  
   h2_smearing_incident_interE_1->GetYaxis()->SetRangeUser(0,1000);
   h2_smearing_incident_interE_1->Draw();

   c1_incident_interE->cd(2);
   c1_incident_interE->cd(2)->SetGrid();
   //Get Mean
   TH2D *h2_smearing_incident_interE_2 = (TH2D*)output->Get("h2_smearing_incident_interE_2");
   h2_smearing_incident_interE_2->SetTitle("True Pion interacting Energy, fit Sigma;reco interactingKE;");  
   h2_smearing_incident_interE_2->GetYaxis()->SetRangeUser(0,100);
   h2_smearing_incident_interE_2->Draw();

   TCanvas *c2_incident_interE = new TCanvas("c2_incident_interE", "");
   c2_incident_interE->SetGrid();
   h2_smearing_incident_interE->Draw("COLZ");

//----------------------------------------------------------------

   TCanvas *c1_abs_int = new TCanvas("c1_abs_int", "");
   c1_abs_int->Divide(1,2);
   c1_abs_int->cd(1);
   c1_abs_int->cd(1)->SetGrid();
   
   //Get Mean
   TH2D *h2_smearing_abs_int_1 = (TH2D*)output->Get("h2_smearing_abs_int_1");
   h2_smearing_abs_int_1->SetTitle("True Absorption interacting Energy, fit Mean;reco interactingKE;");
   h2_smearing_abs_int_1->GetYaxis()->SetRangeUser(0,1000);
   h2_smearing_abs_int_1->Draw();

   c1_abs_int->cd(2);
   c1_abs_int->cd(2)->SetGrid();
   //Get Mean
   TH2D *h2_smearing_abs_int_2 = (TH2D*)output->Get("h2_smearing_abs_int_2");
   h2_smearing_abs_int_2->SetTitle("True Absorption interacting Energy, fit Sigma;reco interactingKE;");  
   h2_smearing_abs_int_2->GetYaxis()->SetRangeUser(0,100);
   h2_smearing_abs_int_2->Draw();

   TCanvas *c2_abs_int = new TCanvas("c2_abs_int", "");
   c2_abs_int->SetGrid();
   h2_smearing_abs_int->Draw("COLZ");
//----------------------------------------------------------------

   TCanvas *c1_totInel_int = new TCanvas("c1_totInel_int", "");
   c1_totInel_int->Divide(1,2);
   c1_totInel_int->cd(1);
   c1_totInel_int->cd(1)->SetGrid();
   
   //Get Mean
   TH2D *h2_smearing_totInel_int_1 = (TH2D*)output->Get("h2_smearing_totInel_int_1");
   h2_smearing_totInel_int_1->SetTitle("True PionInelastic interacting Energy, fit Mean;reco interactingKE;");  
   h2_smearing_totInel_int_1->GetYaxis()->SetRangeUser(0,1000);
   h2_smearing_totInel_int_1->Draw();

   c1_totInel_int->cd(2);
   c1_totInel_int->cd(2)->SetGrid();
   //Get Mean
   TH2D *h2_smearing_totInel_int_2 = (TH2D*)output->Get("h2_smearing_totInel_int_2");
   h2_smearing_totInel_int_2->SetTitle("True PionInelastic interacting Energy, fit Sigma;reco interactingKE;"); 
   h2_smearing_totInel_int_2->GetYaxis()->SetRangeUser(0,100),
   h2_smearing_totInel_int_2->Draw();

   TCanvas *c2_totInel_int = new TCanvas("c2_totInel_int", "");
   c2_totInel_int->SetGrid();
   h2_smearing_totInel_int->Draw("COLZ");

//----------------------------------------------------------------
return 0;
};


