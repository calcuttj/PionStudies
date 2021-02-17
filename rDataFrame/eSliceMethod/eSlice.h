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


/*auto fillIncidentHisto = [h_true_pion_true_incidentE, &bin_size](double initialEnergy, double interactingEnergy){

      //bin that energy falls into is (int) energy/nbins + 1

      int binNumber_initEnergy = (int) initialEnergy / bin_size + 1;
      std::cout << "binNumber initEnergy " << binNumber_initEnergy << std::endl;

      int binNumber_interEnergy = (int) interactingEnergy / bin_size + 1;
      std::cout << "binNumber interEnergy " << binNumber_interEnergy << std::endl;  

      //how many bins need to be filled
      //start from initial energy value and subtract nBin at every time
      std::cout << "initial Energy = " << initialEnergy << std::endl;
      int cnt = 0;
      for(size_t i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){

         h_true_pion_true_incidentE->Fill( initialEnergy - cnt*bin_size);
         std::cout << initialEnergy - cnt*bin_size << std::endl;
         cnt++;
      };
      
      std::cout << "interacting Energy = " << interactingEnergy << std::endl;
};*/

/*auto bla = [&bin_size](double intialE){
      bin_size = bin_size + initialE;
};
*/

