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
//#include "../betheBloch.h"

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
//**************************
//
//Histo Binning
//
//**************************
//binning, bc of beam spread we need to weight the bins, for incident Histo bin very fine and create it, rebin afterwards and divide by the rebinned number (otherwise it is multiple counting


double bin_size_int = 20.;
double bin_size_inc = 1.;
double eStart = 1200.;
double eEnd = 0.;
int nBin_int = (eStart - eEnd) / bin_size_int;
int nBin_inc = (eStart - eEnd) / bin_size_inc;


//**************************
//
//ARGON PARAMETERS
//
//**************************

double atomic_mass = 39.948;
double density = 1.4;
double N_avogadro = 6.02*pow(10,23);
double factor_mbarn = pow(10,27);

//**************************
//
// FIRST Incident Energy Entry
//
//**************************
//

auto firstIncident = [](const std::vector<double> &incidentEnergy){
   return incidentEnergy[0];
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



