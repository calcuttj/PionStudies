#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TMath.h"
#include <ROOT/RDataFrame.hxx>
#include "TColor.h"

#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <numeric>


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


double bin_size_int = 25.;
double bin_size_inc = 25.; //2

double eStart = 1500.; //1200
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
      
      int birth = initialE->GetBinContent( i );
      int death = interE->GetBinContent( i );

      for(int j = i -1; j >= 1; j--){
         incident->SetBinContent( j, incident->GetBinContent(j) + birth );
      };

      for(int j = i -1; j >= 1; j--){
         incident->SetBinContent( j, incident->GetBinContent(j) - death );
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


//----------------------------------------
// Build Bethe Hist for XS calculation
//----------------------------------------

void fill_betheHisto( TH1D* bethe_hist, double mass){

   int nBin = bethe_hist->GetNbinsX();

   for(int i = 1; i <= nBin; i ++){
      bethe_hist->SetBinContent( i, betheBloch( eEnd + (i - 0.5)*bin_size_int, mass) ); 
   };

};


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



