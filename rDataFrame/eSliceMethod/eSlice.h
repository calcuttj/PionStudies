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
#include <ROOT/RDataFrame.hxx>


#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include <vector>


using namespace std;
using namespace ROOT::VecOps;

//definitions for eSlice Method
//
double mass_pion = 139.; // MeV
double mass_muon = 105.6;
double P_in_pion = 1000.; //MeV
double P_in_muon = 1000.; //MeV
double KE_in_pion = sqrt( pow(P_in_pion,2) + pow(mass_pion,2) ) - mass_pion;
double KE_in_muon = sqrt( pow(P_in_muon,2) + pow(mass_muon,2) ) - mass_muon;
//----------------------------------------
//
// E loss
//
//----------------------------------------
double eLoss_mc_trueE = 43; //MeV, kinetic energy
//----------------------------------------
//
//Histo Binning
//
//----------------------------------------
//binning, bc of beam spread we need to weight the bins, for incident Histo bin very fine and create it, rebin afterwards and divide by the rebinned number (otherwise it is multiple counting


double bin_size_int = 20.;
double bin_size_inc = 20.; //2

double eStart = 1200.; //1200
double eEnd = 0.;
int nBin_int = (eStart - eEnd) / bin_size_int;
int nBin_inc = (eStart - eEnd) / bin_size_inc;


//----------------------------------------
//
//ARGON PARAMETERS
//
//----------------------------------------

double atomic_mass = 39.948;
double density = 1.4;
double N_avogadro = 6.02*pow(10,23);
double factor_mbarn = pow(10,27);
double scale_factor = factor_mbarn * atomic_mass / ( density * N_avogadro * bin_size_int );

//----------------------------------------
//
// FIRST Incident Energy Entry
//
//----------------------------------------
//

auto firstIncident = [](const std::vector<double> &incidentEnergy){
   return incidentEnergy[0];
};

//----------------------------------------
//Function to BUILD Incident Histo
//from Initial E distribution and Inter E distribution
//
//----------------------------------------
void build_incidentHist(TH1D* initialE, TH1D* interE, TH1D* incident){

   int nBin = initialE->GetNbinsX();

   for(int i = nBin; i>=1; i--){

      double birth = initialE->GetBinContent( i );
      double birth_err = initialE->GetBinError( i );
      double death = interE->GetBinContent( i );
      double death_err = interE->GetBinError( i );

      for(int j = i -1; j >= 1; j--){
         incident->SetBinContent( j, incident->GetBinContent(j) + birth );
         incident->SetBinError( j, incident->GetBinError(j) + birth_err );
      };

      for(int j = i -1; j >= 1; j--){
         incident->SetBinContent( j, incident->GetBinContent(j) - death );
         incident->SetBinError( j, incident->GetBinError(j) - death_err );
      };

   };

   return;

};

auto equalBin(double init_KE, double inter_KE){
   int bin_initE = (int) init_KE / bin_size_inc + 1;
   int bin_interE = (int) inter_KE / bin_size_inc + 1;

   if(bin_initE == bin_interE) return true;
   else return false;

};

//----------------------------------------
// Function to check nitial E bin and Interacting Ebin
///----------------------------------------

auto checkBins(double init_KE, double inter_KE, int bin_init, int bin_inter){

   if(
         bin_init != bin_inter && 
         init_KE > inter_KE &&
         bin_init <= nBin_int &&
         bin_init > 0 &&
         bin_inter <= nBin_int &&
         bin_inter > 0 ){

      return true;
   }

   else return false;
};

//----------------------------------------
// Function to fill InitE, InterE histo and Interacting Histo
//----------------------------------------
void fill_initE_interE(TH1D* initE, TH1D* interE, double init_KE, double inter_KE){
   //make sure incident Pion does not interact in bin it was born
   int bin_initE = (int) init_KE / bin_size_inc + 1;
   int bin_interE = (int) inter_KE / bin_size_inc + 1;
   if( checkBins(init_KE, inter_KE, bin_initE, bin_interE) ){

      initE->SetBinContent( bin_initE, initE->GetBinContent( bin_initE ) + 1); 
      interE->SetBinContent( bin_interE, interE->GetBinContent( bin_interE ) + 1);    
   }
   return;
};

void fill_interacting(TH1D* interE, double init_KE, double inter_KE){
   //make sure incident Pion does not interact in bin it was born
   int bin_initE = (int) init_KE / bin_size_inc + 1;
   int bin_interE = (int) inter_KE / bin_size_inc + 1;
   if( checkBins(init_KE, inter_KE, bin_initE, bin_interE) ){
      interE->SetBinContent( bin_interE, interE->GetBinContent( bin_interE ) + 1);    
   }
   return;
};

//----------------------------------------
// Function to fill smearing histos initE and interE
//----------------------------------------

void fill_smearing_matrix(TH2D* h2_smear, double true_initKE, double true_interKE, double reco_initKE, double reco_interKE, bool doInter){
   //make sure incident Pion does not interact in bin it was born
   int bin_true_initE = (int) true_initKE / bin_size_inc + 1;
   int bin_true_interE = (int) true_interKE / bin_size_inc + 1;

   int bin_reco_initE = (int) reco_initKE / bin_size_inc + 1;
   int bin_reco_interE = (int) reco_interKE / bin_size_inc + 1;

   if( checkBins(true_initKE, true_interKE, bin_true_initE, bin_true_interE) && checkBins(reco_initKE, reco_interKE, bin_reco_initE, bin_reco_interE) ){

      //depends if initE smearing or interE smearing, maybe if we don't do unsmearing on initialE this can be taken away
      if(doInter) h2_smear->Fill(reco_interKE, true_interKE);
      else h2_smear->Fill(reco_initKE, true_initKE);

   }; 
   return;
};

void fill_smearing_matrix_inter(TH2D* h2_smear, double true_initKE, double true_interKE, double reco_initKE, double reco_interKE ){
   //make sure incident Pion does not interact in bin it was born
   int bin_true_initE = (int) true_initKE / bin_size_inc + 1;
   int bin_true_interE = (int) true_interKE / bin_size_inc + 1;

   int bin_reco_initE = (int) reco_initKE / bin_size_inc + 1;
   int bin_reco_interE = (int) reco_interKE / bin_size_inc + 1;

   if( checkBins(true_initKE, true_interKE, bin_true_initE, bin_true_interE) && checkBins(reco_initKE, reco_interKE, bin_reco_initE, bin_reco_interE) ){

      h2_smear->Fill(reco_initKE, true_initKE);

   }; 
   return;
};

//----------------------------------------
// Normalise Smearing histos, true is bin Y, reco is bin X
//----------------------------------------

void normalise_smearing(TH2D* smear){

   int nBin = smear->GetNbinsX();

   for(int i=1; i <= smear->GetNbinsX(); i++){

      int sum = smear->Integral( 1, nBin, i, i ); //sum up all recoE signals for one trueE signal

      if(sum !=0){
         for(int j = 1; j <= nBin; j++){

            smear->SetBinContent( j , i, smear->GetBinContent(j,i) / sum);

         };
      };
   }; 


   return;
};

//----------------------------------------
// CONSTRUCT XS from Histos 
//----------------------------------------

void do_XS_log(TH1D* xs, TH1D* interacting, TH1D* incident, TH1D* hist_bethe){

   TH1D* temp = (TH1D*) incident->Clone();
   temp->Add( interacting, -1);

   xs->Divide( incident, temp);
   for(int i = 1; i <= xs->GetNbinsX(); i++){

      if( xs->GetBinContent(i) > 0) xs->SetBinContent(i, log( xs->GetBinContent(i) ) );

      else xs->SetBinContent(i, -10);
   };

   xs->Multiply( hist_bethe);
   xs->Scale( scale_factor );

   return;
};

void do_XS_log_binomial_error(TH1D* xs, TH1D* interacting, TH1D* incident, TH1D* hist_bethe){

   for(int i = 1; i <= xs->GetNbinsX(); i++){

      double p = interacting->GetBinContent(i) / incident->GetBinContent(i);
      double nInc = incident->GetBinContent(i);
      double factor = scale_factor * hist_bethe->GetBinContent(i);

      xs->SetBinError( i, factor * sqrt( p * (1-p) / nInc ) );
   };


};

//----------------------------------------------------------
//Functions for Bethe Values MPV and Mean
//
void hist_bethe_mpv(double E_init, double mass_particle, TH1D* fit_mean, TH1D* h_bethe ){
   //initial beam energy
   int cnt =0;
   for(int i=1; i <= fit_mean->GetNbinsX(); i++){
      auto temp = fit_mean->GetBinContent(i);
      //should be from wire 68 on
      //first wire with hit
      //if(temp > 0) {
      h_bethe->SetBinContent( i, (1/0.51)*betheBloch_mpv(E_init, mass_particle)); 
      h_bethe->SetBinError(i, 0.001 );
      E_init = E_init - betheBloch(E_init, mass_particle)*0.51; 

      //0.51 is the pitch wrt beam angle for beam particles

      //energy at each passage is reduced by mean value of bethe bloch
      //}
      if(E_init <= 0) return;
   };
};

void hist_bethe_mean(double E_init, double mass_particle, TH1D* fit_mean, TH1D* h_bethe ){
   //initial beam energy
   int cnt =0;
   for(int i=1; i <= fit_mean->GetNbinsX(); i++){
      auto temp = fit_mean->GetBinContent(i);
      //should be from wire 68 on
      //first wire with hit
      //if(temp > 0) {
      h_bethe->SetBinContent( i, betheBloch(E_init, mass_particle)); 
      h_bethe->SetBinError(i, 0.001 );
      E_init = E_init - betheBloch(E_init, mass_particle)*0.51; //should put pitch apparent at that point
      //0.51 is the pitch wrt beam angle for beam particles


      //energy at each passage is reduced by mean value of bethe bloch
      //}

      if(E_init <= 0) return;

   };
};

void hist_bethe_mean_distance(double E_init, double mass_particle, TH1D* h_bethe ){

   for(int i=1; i <= h_bethe->GetNbinsX(); i++){
      h_bethe->SetBinContent( i, betheBloch(E_init, mass_particle)); 
      h_bethe->SetBinError(i, 0.001 );
      E_init = E_init - betheBloch(E_init, mass_particle);       
      if(E_init <= 0) return;

   };
};


//----------------------------------------
// Build Bethe Hist for XS calculation
//----------------------------------------

void fill_betheHisto( TH1D* bethe_hist, double mass){

   int nBin = bethe_hist->GetNbinsX();

   for(int i = 1; i <= nBin; i ++){
      bethe_hist->SetBinContent( i, betheBloch( eEnd + (i - 0.5)*bin_size_int, mass) ); 
   };

};


//----------------------------------------
//   MAKE INVERSE MATRIX
//----------------------------------------
// input should be smearing matrix, output is inverse matrix as TH2D
//
void invert_smearing( TH2D* smearing, TH2D* inverse){

   int nBin = smearing->GetNbinsX();
    
   TMatrixD matrix_help(nBin + 2, nBin + 2, smearing->GetArray(), "D");
   TMatrixD matrix = matrix_help.GetSub(1, nBin, 1, nBin); //avoid underflow and overflowbin

   //matrix_smearing_incident_initE.Print();

   Double_t det1;
   TMatrixD matrix_inverse = matrix;
   matrix_inverse.Invert(&det1);

   TMatrixD U1(matrix_inverse, TMatrixD::kMult, matrix);
   TMatrixDDiag diag1(U1); diag1 = 0.0;
   const Double_t U1_max_offdiag = (U1.Abs()).Max();
   std::cout << "  Smearing Function "  << std::endl;
   std::cout << "  Maximum off-diagonal = " << U1_max_offdiag << std::endl;
   std::cout << "  Determinant          = " << det1 << std::endl;   

   //put inverse into TH2
   for (int i = 1; i <= nBin; i++){
      for (int j= 1; j <= nBin; j++){
         inverse->SetBinContent(j, i, matrix_inverse(i-1,j-1)); //vector indices style for matrix 
      };
   };  
   return;
};

/*   //Check that matrix and inverse are unity
   TMatrixD unity_smearing_incident_initE = matrix_inverse_smearing_incident_initE * matrix_smearing_incident_initE ;
   TH2D* h2_prod_smearing_incident_initE = new TH2D("h2_prod_smearing_incident_initE", "", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);
   for (int i = 1; i <= nBin_int; i++){
      for (int j= 1; j <= nBin_int; j++){
         h2_prod_smearing_incident_initE->SetBinContent(j, i, unity_smearing_incident_initE(i-1,j-1)); 
      }
   }

   h2_prod_smearing_incident_initE->Write();

*/




auto deltaE = [](const std::vector<double> &dEdX, const std::vector<double> &pitch){

   std::vector<double> dEnergy;
   size_t i = 0;
   while(i < dEdX.size() && i < pitch.size()){
      dEnergy.push_back(dEdX[i]*pitch[i]);
      //std::cout << "dEnergy = " << dEdX[i]*pitch[i] << std::endl;
      i++;
   }
   return dEnergy;
};

auto relPos = [](const std::vector<double> &dEdX, const std::vector<double> &pitch){

   std::vector<double> relPos;
   size_t i = 0;
   double cnt = 0.;
   //std::cout << "dEdX size = " << dEdX.size() << std::endl;
   //std::cout << "pitch size = " << pitch.size() << std::endl;
   while(i < dEdX.size() && i < pitch.size()){

      relPos.push_back( cnt / (double) dEdX.size() );
      //std::cout << "cnt = " << cnt << std::endl;
      //std::cout << "dEdX size = " << (double) dEdX.size() << std::endl;
      //std::cout << "relPos = " << cnt / (double) dEdX.size() << std::endl;
      cnt++;
      i++;
   }
   return relPos;
};


