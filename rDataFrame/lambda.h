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
//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;
//lambda.h includes the definitions used for the pionAnalysis 
//
//INPUTfiles
const std::string inputTree = "pionana/beamana";
//**********************************************************
//DEFINITIONS
//
//**********************************************************
//


int compString(std::string s1, std::string s2){

   return (s1 == s2);
   //if(s1 != s2) return 0;
   //if(s1 == s2) return 1;
};

Int_t palette_pid[] = {kRed, kOrange+7, kBlue+2, kRed+3, kGreen+2, kViolet-5, kCyan-7, kCyan+3, kPink+9, kOrange-2 };
Int_t palette_interactions[] = {kRed, kOrange, kBlue, kGreen, kCyan, kMagenta};
Int_t palette_pdg_reduced[] = {kViolet, kCyan, kGreen, kRed, kOrange, kBlue};

//good Reco primary Pions (true MC) with pi+inelastic interaction in the end
//
auto good_reco = [](double quality_backtrack_0, double quality_backtrack_1, double quality_backtrack_2, double quality_maxseg_0, double quality_maxseg_1, double quality_maxseg_2){
   if(quality_backtrack_0 > 15. || quality_backtrack_1 > 15. || quality_backtrack_2 > 15. ) return 0;
   else if (quality_maxseg_0 > 15. || quality_maxseg_1 > 15. || quality_maxseg_2 > 15.) return 0;

   return 1;
};


auto truePrimaryPionInel = [](int reco_beam_true_byHits_PDG, int reco_beam_true_byHits_origin, int good_reco, std::string reco_beam_true_byHits_process, std::string reco_beam_true_byHits_endProcess)
{

   std::string pionInel("pi+Inelastic");
   std::string prim ("primary");

   return good_reco == 1 && reco_beam_true_byHits_PDG == 211 && reco_beam_true_byHits_origin ==4 &&
      compString(reco_beam_true_byHits_process, prim) == 1 && compString(reco_beam_true_byHits_endProcess,pionInel) == 1;
};

//True Charge Exchange + Absorption Signal, has no piPlus or piMinus as daughters
auto trueChExAbsProcess = [](const int true_daughter_nPiPlus, const int true_daughter_nPiMinus,const int true_daughter_nPi0)
{
   return true_daughter_nPiPlus + true_daughter_nPiMinus == 0 && true_daughter_nPi0 < 2;
};

//True Charge Exchange Signal first filter for ChEx + Absoprtion signal, then ask for a pi0
auto trueChExProcess = [](int true_daughter_nPi0) {return true_daughter_nPi0 == 1;};

//True Absorption Signal, first filter for ChEx + Abs Signal, then ask for at least one proton Daughter
//
auto trueAbsProcess = [](int true_daughter_nPi0) { return true_daughter_nPi0 == 0;};

//True Background Signal, NOT Charge Exchange and NOT Absorption
//
auto trueBackGround = [](int true_daughter_nPiPlus, int true_daughter_nPiMinus, int true_daughter_nPi0){
   return !(true_daughter_nPiPlus + true_daughter_nPiMinus == 0 && true_daughter_nPi0 <2);
};

template <class F>
F primary_property(int pdg, int pdg_primary, const F primary_property){
   F property = -999;
   if(pdg == pdg_primary){
      property = primary_property;
   }
   return property;
};

//Define all particle types
//I haven't yet figured out how anything else than columns can be passed to lambdas. so for now this helps a define function to define a new column "proton" etc for the particles, filled with that value and then use it to pass to a function that will take that value and compare it to a PDG value
auto pdg_proton = [](){return 2212;};
auto pdg_piPlus = [](){return 211;};
auto pdg_piMinus = [](){return -211;};
auto pdg_muMinus = [](){return 13;};
auto pdg_muPlus = [](){return -13;};
auto pdg_kaon = [](){return 321;};
auto pdg_gamma = [](){return 22;};
auto pdg_electron = [](){return 11;};
auto pdg_positron = [](){return -11;};
auto pdg_nucleus = [](){return 9999;}; //watchout for nucleuss in comparing function
auto pdg_neutron = [](){return 2112;};

//COUNT particle Type in Daughters
auto count_type = [](int pdg, const std::vector<int> &pdg_vec){
   int cnt = 0;
   for(size_t pos = 0; pos < pdg_vec.size(); pos++){
      if( pdg != 9999 && pdg == pdg_vec[pos]) cnt++;
      else if (pdg == 9999 && pdg_vec[pos] > 3000) cnt++;
   };
   return cnt;
};

auto count_pi0_gamma = [](const std::vector<int> &pi0_gamma_ID, const std::vector<int> &reco_true_ID){
   int cnt_gamma = 0;
   if(pi0_gamma_ID.size() != 0){
      for(size_t pos = 0; pos < pi0_gamma_ID.size(); pos++){
         if(pos < reco_true_ID.size()){
            for(size_t cnt = 0; cnt < reco_true_ID.size(); cnt++){
               if(pi0_gamma_ID[pos] == reco_true_ID[cnt]) cnt_gamma++;
            };
         }
      };
   }
   return cnt_gamma;
};

//Find properties (stored in vector) of a specific daughter particle, special for nucleus daughters
template <class T>
T daughter_property(int pdg, const std::vector<int> &pdg_vec, const T &daughter_property){
   T return_vec; 
   for (size_t pos =0; pos < pdg_vec.size(); pos++){ 
      if(pdg!= 666 && pdg!= 9999 && pdg_vec[pos] == pdg && daughter_property.size()> pos){
         return_vec.push_back(daughter_property[pos]);}
      //for nuclei
      else if(pdg == 9999 && pdg_vec[pos] > 3000 && daughter_property.size() > pos) {
         return_vec.push_back(daughter_property[pos]);}
      //others than proton pion and photon
      else if(pdg == 666 && pdg_vec[pos] != 22 && pdg_vec[pos] != 211 && pdg_vec[pos] != 2212 && daughter_property.size() > pos){
         return_vec.push_back(daughter_property[pos]);};
   };
   return return_vec;
};

//shower property
template <class F>
F shower_property(const std::vector<double> &trackScore, const F &daughter_property){
   F return_vec; 
   for (size_t pos =0; pos < trackScore.size(); pos++){ 
      if(trackScore[pos] < 0.3 && trackScore[pos] != -999){
         return_vec.push_back(daughter_property[pos]);}
         };
   return return_vec;
};

template <class A>
A pion_daughter_property( const std::vector<int> &pdg_vec, const A &daughter_property){
   A return_vec; 
   for (size_t pos =0; pos < pdg_vec.size(); pos++){ 
      if(abs(pdg_vec[pos]) == 211 && daughter_property.size()> pos){
         return_vec.push_back(daughter_property[pos]);
      }
   };
   return return_vec;
};

template <class S>
S pi0_gamma_property (const std::vector<int> &pi0_gamma_ID, const std::vector<int> &reco_daugh_ID, const S &daughter_property){
   S return_vec;
   if(pi0_gamma_ID.size() != 0){

      for(size_t pos = 0; pos < pi0_gamma_ID.size(); pos++){
         if(pos < reco_daugh_ID.size()){
            for(size_t cnt = 0; cnt < reco_daugh_ID.size(); cnt++){
               if(pi0_gamma_ID[pos] == reco_daugh_ID[cnt] && cnt < daughter_property.size() ) return_vec.push_back(daughter_property[cnt]);
            };
         }
      };
   }
   return return_vec;
};

template <class B>
B nuclear_gamma_property (const std::vector<int> &pi0_gamma_ID, const std::vector<int> &reco_daugh_PDG, const std::vector<int> &reco_daugh_ID, const B &daughter_property){
   B return_vec;
   if(pi0_gamma_ID.size() != 0){

      for(size_t pos = 0; pos < pi0_gamma_ID.size(); pos++){
         if(pos < reco_daugh_ID.size()){
            for(size_t cnt = 0; cnt < reco_daugh_ID.size(); cnt++){
               if(pi0_gamma_ID[pos] != reco_daugh_ID[cnt] && reco_daugh_PDG[cnt] == 22 && cnt < daughter_property.size()) return_vec.push_back(daughter_property[cnt]);
            };
         }
      };
   }
   else {
      for(size_t cnt = 0; cnt < reco_daugh_PDG.size(); cnt++){
         if(reco_daugh_PDG[cnt] == 22 && cnt < daughter_property.size()) return_vec.push_back(daughter_property[cnt]);
      };
   }
   return return_vec;
};


//Find properties of an event if one of the daughters is a specific particle
auto event_property = [](int pdg, const std::vector<int> &pdg_vec, int ev_property){
   int return_value = 0; 
   for (size_t pos =0; pos < pdg_vec.size(); pos++){ 
      if(pdg!= 9999 && pdg_vec[pos] == pdg ){
         return_value = ev_property;
         break;}
      else if(pdg == 9999 && pdg_vec[pos] > 3000 ) {
         return_value = ev_property;
         break;};
   };
   return return_value;
};




//Means of dEdX Vector<Vector<double>> for the Daughter Particles
auto meanDaughterdEdX = [](const ROOT::RVec<std::vector<double>> &vecs){
   std::vector<double> means; 
   for (auto &&vec : vecs) means.push_back(accumulate(vec.begin(),vec.end(),0.0)/vec.size()); 
   return means;
};

auto medianDaughterdEdX = [](std::vector<std::vector<double>> &vecs){
   std::vector<double> return_vec;
   std::vector<double> help_vec;
   size_t size = 0;
   
   for(size_t i=0; i < vecs.size(); i++){

      help_vec.clear();

      for(size_t j=0; j < vecs[i].size(); j++){

         help_vec.push_back(vecs[i][j]);
      };

      sort( help_vec.begin(), help_vec.end() );
      size = help_vec.size();
      
      if(size == 0){
         return_vec.push_back(-9999);
         continue;
      }

      if (size % 2 == 0){
        return_vec.push_back( (help_vec[size / 2 - 1] + help_vec[size / 2]) / 2);
      }
      else{
         return_vec.push_back( help_vec[size / 2] );
      }

   };
   
   return return_vec;

};


//fill matrix Values into vector
auto matrix_to_vector = [](std::vector<std::vector<double>> &matrix){
   std::vector<double> return_vec;
   for(auto &&vec : matrix) {
      for(auto i : vec) {
         return_vec.push_back(i);
      };
   };
   return return_vec;
};


//truncated mean of SIGMA = cutting %
auto truncatedMean = [](double truncate_low, double truncate_high, std::vector<std::vector<double>> &vecs_dEdX){

   size_t size = 0;
   std::vector<double> trunc_mean;
   std::vector<double> help_vec;
   truncate_high = 1 - truncate_high; 
   int i_low = 0;
   int i_high = 0;

   //sort the dEdX vecotrs in matrix
   for(auto &&vec : vecs_dEdX){
      size = vec.size();
      help_vec.clear();

      //check dEdX vector isn't empty!
      if(vec.empty()){
         trunc_mean.push_back(-9999.);
         continue;
      }

      else{
         //Sort Vector
         sort(vec.begin(), vec.end());
       
         //Discard upper and lower part of signal
         //rint rounds to integer
         i_low = rint ( size*truncate_low);
         i_high = rint( size*truncate_high);
         
         
         for(int i = i_low; i <= i_high; i++){
               help_vec.push_back(vec[i]);
         };

         //Mean of help vector

         trunc_mean.push_back(accumulate(help_vec.begin(), help_vec.end(), 0.0) / help_vec.size());

      }


   };

   return trunc_mean;

};

auto truncatedMean_pos = [](double truncate_low, double truncate_high, std::vector<std::vector<double>> &vecs_dEdX){

   size_t size = 0;
   std::vector<double> trunc_mean;
   //std::vector<double> help_vec;
   truncate_high = 1 - truncate_high; 
   int i_low = 0;
   int i_high = 0;

   //sort the dEdX vecotrs in matrix
   for(auto &&vec : vecs_dEdX){
      //size = vec.size();
      std::vector<double> help_vec;


      //check dEdX vector isn't empty!
      if(vec.empty()){
         trunc_mean.push_back(-9999.);
         continue;
      }

      else{
         std::vector<double> temp_vec;
         for (double d : vec) {
           if (d < 0.) {
             continue;
           }
           temp_vec.push_back(d);
         }
         for (double d : temp_vec) {
           if (d < 0.) {
             std::cout << d << std::endl;
           }
         }
         if (temp_vec.empty()) {
           trunc_mean.push_back(-9999.);
           continue;
         }
         //Sort Vector
         sort(temp_vec.begin(), temp_vec.end());
       
         //Discard upper and lower part of signal
         //rint rounds to integer
         i_low = rint ( temp_vec.size()*truncate_low);
         i_high = rint( temp_vec.size()*truncate_high);
         
         
         //if (i_high >= temp_vec.size()) std::cout << "Warning: too high" << std::endl;
         for(int i = i_low; i </*=*/ i_high; i++){
           if (temp_vec[i] < 0) {
             std::cout << "added " << temp_vec[i] << " " << i << std::endl;
           }
           help_vec.push_back(temp_vec[i]);
         };

         //Mean of help vector

         trunc_mean.push_back(accumulate(help_vec.begin(), help_vec.end(), 0.0) / help_vec.size());
         if (trunc_mean.back() < 0 && trunc_mean.back() != -9999. && trunc_mean.back() != -999.) {
           std::cout << accumulate(help_vec.begin(), help_vec.end(), 0.0) << " " << help_vec.size() << std::endl;
           std::cout << temp_vec.size() << " " << i_low << " " << i_high << std::endl;
           for (size_t i = 0; i < help_vec.size(); ++i) {
             std::cout << "\t" << help_vec[i] << std::endl;
           }
         }
      }


   };

   return trunc_mean;

};

auto recoLength = [](const std::vector<double> &sX, const std::vector<double> &sY, const std::vector<double> &sZ, const std::vector<double> &eX, const std::vector<double> &eY, const std::vector<double> &eZ){
   std::vector<double> length;
   double dX = 0, dY = 0, dZ = 0;

   if( sX.size() == eZ.size()){

      for(size_t pos = 0; pos < sX.size(); pos++){
         dX = sX[pos] - eX[pos];
         dY = sY[pos] - eY[pos];
         dZ = sZ[pos] - eZ[pos];
         length.push_back(sqrt(dX*dX + dY*dY + dZ*dZ));
      };

   }
   else length.push_back(-999);
   return length;
};

auto select_daughter_matrix = [](int pdg, const std::vector<int> &pdg_vec, const std::vector<std::vector<double>> &matrix){

   std::vector<std::vector<double>> return_matrix;
   std::vector<double> help;
   for(size_t i=0; i < pdg_vec.size(); ++i){
      if(abs(pdg_vec[i]) == pdg) {

         for(size_t n=0; n < matrix[i].size(); ++n){
            help.push_back(matrix[i][n]);
         }

         return_matrix.push_back(help);
         help.clear();     
      };
   }
   return return_matrix;

};

auto matrix_vector = [](const std::vector<std::vector<double>> &matrix, const std::vector<std::vector<double>> &matrix_help){

   std::vector<double> return_vec;
   size_t temp;
   for(size_t i=0; i < matrix.size(); ++i){

      //to ensure that vectors have same length (resRange dEdX)
      if(matrix[i].size() <= matrix_help[i].size()) temp = matrix[i].size();
      else temp = matrix_help[i].size();

      for(size_t n=0; n < temp; ++n){
         return_vec.push_back(matrix[i][n]);
      }
   }

   return return_vec;
};

auto doChi2 = [](std::vector<double> &chi2, std::vector<int> &ndof){

   std::vector<double> return_vec;
   if(chi2.size() != 0){
      for(std::string::size_type pos = 0; pos < chi2.size(); pos++){
         return_vec.push_back(chi2.at(pos)/ndof.at(pos));
      }
      return return_vec;
   }

   return return_vec;
};


auto rad_deg = [](std::vector<double> &rad){
      std::vector<double> deg;
      for(size_t i=0; i<rad.size(); i++){
         deg.push_back(rad[i]*(180/3.14159));
      }
      return deg;
   };

auto energyDeposition = [](std::vector<double> &dEdX, std::vector<double> &len){

      std::vector<double> energy;
      for(size_t i=0; i< dEdX.size(); i++){
         energy.push_back( dEdX[i] * len[i]);
      };

      return energy;

   };



