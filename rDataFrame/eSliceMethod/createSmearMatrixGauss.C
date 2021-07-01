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
//After doing Fits of the Matrix reco-trueE interacting and initial one can proceed to create the smearing matrix with the values found for gaussian fits
//
//TODAY 22/06/21
//Used smearing_gaus_test.C and binSize of 20MeV for recoKE(x) - trueKE(y)
//THE FITS DEPEND ON BIN SIZE, the larger the bin size, the larger the sigmas of the fit!
//------------------------
//Interacting Energy
//Fitted Slices in Y (true)
//mean linear slope 1, sigma at around 30MeV
//------------------------
//initial Energy
//Fitted Slices in Y (true)
//mean linear slope 1, sigma at around 23MeV
//------------------------
//Decided to consider interval of +-3sigma --> -3*30 to + 3*30
//means to fill diagonal bin +- 4*bins
//all diagonal bins will be the same, so only need to produce one row and then just copy the values into the correct place for the other rows
//Strategy, create the filled bins for the normal distribution with mean = 0;
//Bin value is the integral of the standard normal function from bin-start to bin-end
//       --> for a gauss with mean = 0 and any sigma
//       probability of finding a value inside interval [-n*sigma, +n*sigma] = erf(n/sqrt(2))
//       can calculate probability for all bins from diag-4 to diag+4
//       --------------------------------------------
//!!!!!!!n depends on bin width of the smearing matrix!!! !!!!!!!!!!!!!
//--------------------------------------------------
//

int createSmearMatrixGauss(){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");

   string output_name = "smearMatrixGauss_" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile ( output_name.c_str() , "RECREATE");

   TH2D* smearMatrix_interacting = new TH2D("smearMatrix_gaussFit_interacting", "Smearing Matrix for Interacting E Pions, created from Gauss Fits; recoInteractingKE; trueInteractingKE", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   TH2D* smearMatrix_initial = new TH2D("smearMatrix_gaussFit_initial", "Smearing Matrix for Initial E Pions, created from Gauss Fits; recoInitialKE; trueInitialKE", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   TH1D* monoChromatic_trueKE_250 = new TH1D("monoChromatic_trueKE_250", "MonoChromatic Energy;trueIntearctingKE", nBin_int, eEnd, eStart);
   TH1D* monoChromatic_trueKE_460 = new TH1D("monoChromatic_trueKE_460", "MonoChromatic Energy;trueIntearctingKE", nBin_int, eEnd, eStart);
   TH1D* monoChromatic_trueKE_580 = new TH1D("monoChromatic_trueKE_580", "MonoChromatic Energy;trueIntearctingKE", nBin_int, eEnd, eStart);

   monoChromatic_trueKE_250->SetBinContent(14,3000);
   monoChromatic_trueKE_460->SetBinContent(24,3000);
   monoChromatic_trueKE_580->SetBinContent(30,3000);

   TH1D* monoChromatic_recoKE_250 = new TH1D("monoChromatic_recoKE_250", "MonoChromatic Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   TH1D* monoChromatic_recoKE_460 = new TH1D("monoChromatic_recoKE_460", "MonoChromatic Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   TH1D* monoChromatic_recoKE_580 = new TH1D("monoChromatic_recoKE_580", "MonoChromatic Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   
   TH1D* randNoise_recoKE_250 = new TH1D("randNoise_recoKE_250", "randNoise Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   TH1D* randNoise_recoKE_460 = new TH1D("randNoise_recoKE_460", "randNoise Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   TH1D* randNoise_recoKE_580 = new TH1D("randNoise_recoKE_580", "randNoise Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   
   TH1D* testSmearBack_recoKE_250 = new TH1D("testSmearBack_recoKE_250", "testSmearBack Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   TH1D* testSmearBack_recoKE_460 = new TH1D("testSmearBack_recoKE_460", "testSmearBack Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   TH1D* testSmearBack_recoKE_580 = new TH1D("testSmearBack_recoKE_580", "testSmearBack Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   
   TH1D* noiseSmearBack_recoKE_250 = new TH1D("noiseSmearBack_recoKE_250", "noiseSmearBack Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   TH1D* noiseSmearBack_recoKE_460 = new TH1D("noiseSmearBack_recoKE_460", "noiseSmearBack Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   TH1D* noiseSmearBack_recoKE_580 = new TH1D("noiseSmearBack_recoKE_580", "noiseSmearBack Energy;recoIntearctingKE", nBin_int, eEnd, eStart);
   //Pre-Work Set all values of TH2D's to 0

   for(int i=1; i <= nBin_int; i++){
      for(int j=1; j <= nBin_int; j++){

         smearMatrix_interacting->SetBinContent(i,j, 0.);
         smearMatrix_initial->SetBinContent(i,j, 0.);

      };
   };

   double sigma_interacting = 30., sigma_initial = 23.;

   int nBin_3sigma_interacting = (int) 3*sigma_interacting / bin_size_int;
   int nBin_3sigma_initial = (int) 3*sigma_initial / bin_size_int;

   //values that will go into matrix, diagonal is erf from -10MeV to + 10MeV, then from -30 to -10 etc
   //n = 10 / 30 = (half the bin-width / sigma)
   //Fill an array with 5 entries, 
   //0th entry is the diagonal probability value
   //1st entry is diag+1 and diag-1 bin entry
   //...
   //4th entry is diag+4 and diag-4 bin entry
   //use array to fill matrix
   //
   std::vector<double> prob_erf_interacting, prob_smearInteracting;
   double n_sigma_interacting;

   for(int i=0; i < nBin_3sigma_interacting + 1; i++){

      n_sigma_interacting = (i*bin_size_int + (bin_size_int / 2)) /sigma_interacting;
      prob_erf_interacting.push_back( TMath::Erf( n_sigma_interacting / sqrt(2) ) );
      //fill array that has binContents for smearing Matrix
      if(i == 0) prob_smearInteracting.push_back( prob_erf_interacting[0] );

      else{
         prob_smearInteracting.push_back( ( prob_erf_interacting[i] -  prob_erf_interacting[i-1] ) / 2 );
      }
      //std::cout << "Probability Erf = " << prob_erf_interacting[i] << std::endl;
      //std::cout << "Probability Bin Matrix = " << prob_smearInteracting[i] << std::endl;

   };

   std::vector<double> prob_erf_initial, prob_smearInitial;
   double n_sigma_initial;

   for(int i=0; i < nBin_3sigma_initial + 1; i++){

      n_sigma_initial = (i*bin_size_int + (bin_size_int / 2)) /sigma_initial;
      prob_erf_initial.push_back( TMath::Erf( n_sigma_initial / sqrt(2) ) );
      //fill array that has binContents for smearing Matrix
      if(i == 0) prob_smearInitial.push_back( prob_erf_initial[0]);

      else{
         prob_smearInitial.push_back( ( prob_erf_initial[i] -  prob_erf_initial[i-1] ) / 2 );
      }
      //std::cout << "Probability Erf = " << prob_erf_initial[i] << std::endl;
      //std::cout << "Probability Bin Matrix = " << prob_smearInitial[i] << std::endl;

   };


   //Fill Interacting Smearing Matrix
   //i goes through rows (true)
   //j goes through the columns (reco)
   for(int i = 1; i <= nBin_int; i++){

      //start at diagnoal bin
      for(int j = 0; j < prob_smearInteracting.size(); j++){

         //fill bins on i + j, up-going bins
         if( i + j <= nBin_int) smearMatrix_interacting->SetBinContent( i + j, i, prob_smearInteracting[j] );

         if( i != j && i - j >= 1) smearMatrix_interacting->SetBinContent( i - j, i, prob_smearInteracting[j] ); 

      }; 

   }
   smearMatrix_interacting->Write();

   //Fill initial Smearing Matrix
   //i goes through rows (true)
   //j goes through the columns (reco)
   for(int i = 1; i <= nBin_int; i++){

      //start at diagnoal bin
      for(int j = 0; j < prob_smearInitial.size(); j++){

         //fill bins on i + j, up-going bins
         if( i + j <= nBin_int) smearMatrix_initial->SetBinContent( i + j, i, prob_smearInitial[j] );

         if( i != j && i - j >= 1) smearMatrix_initial->SetBinContent( i - j, i, prob_smearInitial[j] ); 

      }; 

   }
   smearMatrix_initial->Write();


   //Test Inversion how it looks like!
   TH2D* inverseMatrix_interacting = new TH2D("inverseMatrix_gaussFit_interacting", "Inverse smearing  Matrix for Interacting E Pions, created from Gauss Fits; recoInteractingKE; trueInteractingKE", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   TH2D* inverseMatrix_initial = new TH2D("inverseMatrix_gaussFit_initial", "Inverse smearing Matrix for Initial E Pions, created from Gauss Fits; recoInitialKE; trueInitialKE", nBin_int, eEnd, eStart, nBin_int, eEnd, eStart);

   invert_smearing( smearMatrix_interacting, inverseMatrix_interacting );
   invert_smearing( smearMatrix_initial, inverseMatrix_initial );

   inverseMatrix_interacting->Write();
   inverseMatrix_initial->Write();



   //Test MonoChromatic Energy smear to reco
   //
   //y is fixed for true E
   for(int i=1; i <= nBin_int; i++){

      monoChromatic_recoKE_250->SetBinContent( i ,  monoChromatic_trueKE_250->GetBinContent(14)*smearMatrix_interacting->GetBinContent( i, 14) );
      monoChromatic_recoKE_460->SetBinContent( i ,  monoChromatic_trueKE_460->GetBinContent(24)*smearMatrix_interacting->GetBinContent( i, 24) );
      monoChromatic_recoKE_580->SetBinContent( i ,  monoChromatic_trueKE_580->GetBinContent(30)*smearMatrix_interacting->GetBinContent( i, 30) );

   };

   monoChromatic_recoKE_250->Write();
   monoChromatic_recoKE_460->Write();
   monoChromatic_recoKE_580->Write();

   //Create noise on recoKE histo
   //
   double ran_x, ran_y;
   double binC_250, binC_460, binC_580;
   for(int i=1; i <= nBin_int; i++){

      gRandom->Rannor(ran_x,ran_y);
      
      binC_250 = monoChromatic_recoKE_250->GetBinContent(i);

      if(binC_250 != 0) randNoise_recoKE_250->SetBinContent( i, binC_250 + sqrt(binC_250)*ran_x);
      //else randNoise_recoKE_250->SetBinContent( i, sqrt(1)*ran_x);
 
      binC_460 = monoChromatic_recoKE_460->GetBinContent(i);

      if(binC_460 != 0) randNoise_recoKE_460->SetBinContent( i, binC_460 + sqrt(binC_460)*ran_x);
      else randNoise_recoKE_460->SetBinContent( i, sqrt(1)*ran_x);
     
      binC_580 = monoChromatic_recoKE_580->GetBinContent(i);

      if(binC_580 != 0) randNoise_recoKE_580->SetBinContent( i, binC_580 + sqrt(binC_580)*ran_x);
      else randNoise_recoKE_580->SetBinContent( i, abs(sqrt(1)*ran_x));
      //std::cout << "Random x = " << ran_x << "  Random y = " << ran_y << std::endl;  
   };

   randNoise_recoKE_250->Write();
   randNoise_recoKE_460->Write();
   randNoise_recoKE_580->Write();

   // SMEAR Back Test
   //
   //y is fixed for true E
   for(int i=1; i <= nBin_int; i++){

      double sum_250 = 0, sum_noise_250 = 0;
      double sum_460 = 0, sum_noise_460 = 0;
      double sum_580 = 0, sum_noise_580 = 0;

      for(int j=1; j <= nBin_int; j++){

         sum_250 += monoChromatic_recoKE_250->GetBinContent(j)*inverseMatrix_interacting->GetBinContent( i , j);
         sum_460 += monoChromatic_recoKE_460->GetBinContent(j)*inverseMatrix_interacting->GetBinContent( i , j);
         sum_580 += monoChromatic_recoKE_580->GetBinContent(j)*inverseMatrix_interacting->GetBinContent( i , j);

         sum_noise_250 += randNoise_recoKE_250->GetBinContent(j)*inverseMatrix_interacting->GetBinContent( i , j);
         sum_noise_460 += randNoise_recoKE_460->GetBinContent(j)*inverseMatrix_interacting->GetBinContent( i , j);
         sum_noise_580 += randNoise_recoKE_580->GetBinContent(j)*inverseMatrix_interacting->GetBinContent( i , j);

      };

      testSmearBack_recoKE_250->SetBinContent( i , sum_250 );
      testSmearBack_recoKE_460->SetBinContent( i , sum_460 );
      testSmearBack_recoKE_580->SetBinContent( i , sum_580 );
   
      noiseSmearBack_recoKE_250->SetBinContent( i , sum_noise_250 );
      noiseSmearBack_recoKE_460->SetBinContent( i , sum_noise_460 );
      noiseSmearBack_recoKE_580->SetBinContent( i , sum_noise_580 );
   };

   testSmearBack_recoKE_250->Write();
   testSmearBack_recoKE_460->Write();
   testSmearBack_recoKE_580->Write();
   
   noiseSmearBack_recoKE_250->Write();
   noiseSmearBack_recoKE_460->Write();
   noiseSmearBack_recoKE_580->Write();


      return 0;
};


