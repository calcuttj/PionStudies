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
//using RDataFrame to cut and analyse PionTtrr

using namespace std;
using namespace ROOT::VecOps;

//**********************************************************
//DEFINITIONS
//
//**********************************************************
//
//Some Cut Values
double cutAPA3_Z = 220.; //updated 22aug2021
double cut_trackScore = 0.3;
double cut_michelScore = 0.55;
double cut_dEdX = 3.8;
int cut_nHits_shower_low = 40;
//daughter Distance cut
double cut_daughter_track_distance = 10.;
double cut_daughter_shower_distance_low = 2.;
double cut_daughter_shower_distance_high = 100.;
double cut_secondary_chi2 = 50.;
//Low Energy Values
double energy_limit = 15;
double dEdX_limit = 0.8;

//Daughter Pion Momentum (GeV)
double daughter_pion_momentum = 0.15;

//For MC from Owen Goodwins studies
double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

//For Data from Owen Goodwin
double data_xlow = 0., data_xhigh = 10., data_ylow= -5.;
double data_yhigh= 10., data_zlow=30., data_zhigh=35., data_coslow=.93;

//For Data from Owen Goodwin 
//Updated with Jake 19May2021
double mc_BI_xlow = -2., mc_BI_xhigh = 2., mc_BI_ylow= -1.5;
double mc_BI_yhigh= 2., mc_BI_zlow=28.5, mc_BI_zhigh=31., mc_BI_coslow=.97;

//Data and MC Beam Quality Cuts only from TPC info, Tingjun May 27,2021
//https://indico.fnal.gov/event/49253/contributions/216082/attachments/143665/181907/pioninel.pdf
//
//These are computed AFTER SCE corrections!!
//need to do plots too and cross-check values!
//Francesca values 58XX                                        //TJ values
double data_meanX = -28.6, data_sigmaX = 4.888;                //double data_meanX = -27.91, data_sigmaX = 4.7
double data_meanY = 424.6, data_sigmaY = 5.397;                //double data_meanY = 424.36, data_sigmaY = 5.1
double data_meanZ = 2.707, data_sigmaZ = 1.215;                  //double data_meanZ = 3.78, data_sigmaZ = 1.10;
double data_thetaX = 100.45 * TMath::Pi()/ 180;                //double data_thetaX = 100.45 * TMath::Pi()/ 18
double data_thetaY = 103.52 * TMath::Pi() / 180;               //double data_thetaY = 103.52 * TMath::Pi() / 1
double data_thetaZ = 17.83 * TMath::Pi() / 180;                //double data_thetaZ = 17.83 * TMath::Pi() / 18
                                                               //
                                                               //
double mc_meanX = -30.81, mc_sigmaX = 5.02;                    //double mc_meanX = -30.81, mc_sigmaX = 5.02;
double mc_meanY = 422.41, mc_sigmaY = 4.51;                    //double mc_meanY = 422.41, mc_sigmaY = 4.51;
double mc_meanZ = 0.11, mc_sigmaZ = 0.22;                      //double mc_meanZ = 0.11, mc_sigmaZ = 0.22;
double mc_thetaX = 101.58 * TMath::Pi()/ 180;                  //double mc_thetaX = 101.58 * TMath::Pi()/ 180;
double mc_thetaY = 101.19 * TMath::Pi()/ 180;                  //double mc_thetaY = 101.19 * TMath::Pi()/ 180;
double mc_thetaZ = 16.59 * TMath::Pi()/ 180;                   //double mc_thetaZ = 16.59 * TMath::Pi()/ 180;


double cut_beamQuality_TPC_xyz = 3.;
double cut_beamQuality_TPC_cosTheta = 0.95;

//Tag PrimaryPion without elastic Scattering
//
auto tagPrimPionInel= [](int true_beam_PDG, std::string true_beam_endProcess) {
  return int(true_beam_PDG == 211 && (true_beam_endProcess == "pi+Inelastic") );
};

auto tagPrimPionInel_withElastic = [](int true_beam_PDG,
    std::string true_beam_endProcess, int true_beam_nElasticScatters) {
  return int(true_beam_PDG == 211 && (true_beam_endProcess == "pi+Inelastic") 
      && true_beam_nElasticScatters > 0);
};

//Function to tag if in an event there is no Pion or Pion with momentum lower than daughter_pion_momentum
//for true Signal definition
auto tagDaughterPionMomentumHigh = [](std::vector<int> &true_daughter_PDG, std::vector<double> &true_daughter_startP,
      const int true_daughter_nPiPlus, const int true_daughter_nPiMinus){

   int daughter_pion = 0;
   if(true_daughter_nPiPlus + true_daughter_nPiMinus > 0) {
      for(size_t i=0; i < true_daughter_PDG.size(); i++){
         if( abs(true_daughter_PDG[i]) == 211 && true_daughter_startP[i] > daughter_pion_momentum) return daughter_pion = 1;}
   }
   return daughter_pion;
};

//True Charge Exchange + Absorption Signal, has no piPlus or piMinus daughters with Momentum bigger than daughter_pion_momentum
auto tagAbsChEx = [](int tagPrimPi, int piDaughterHighMomentum) {
  return int(tagPrimPi && piDaughterHighMomentum == 0);
};

//True Charge Exchange Signal first filter for ChEx + Absoprtion signal,
//then ask for a pi0
auto tagChEx = [](int tagAbsChEx, int true_daughter_nPi0) {
  return int( tagAbsChEx && true_daughter_nPi0 == 1 );
};

//True Absorption Signal, first filter for ChEx + Abs Signal,
//then ask for at least one proton Daughter
auto tagAbs = [](int tagAbsChEx, int true_daughter_nPi0) { 
  return int(tagAbsChEx && true_daughter_nPi0 == 0);
};

auto tagNpi0 = [](int tagAbsChEx, int true_daughter_nPi0) {
  return int(tagAbsChEx == 1 && true_daughter_nPi0 > 1);
};

auto tagBackGround = [](int tagPrimPi, int tagAbsChex ) {
  return int( !(tagPrimPi && tagAbsChex) );
};

//Beam Track ends in APA3
auto endAPA3 = [](double reco_beam_endZ) {
  return ( reco_beam_endZ < cutAPA3_Z );
};

auto primary_chi2 = [](double chi2_proton, int chi2_ndof){
  return( (chi2_proton/chi2_ndof > 140.) );
};

//Beam
auto isBeamType = [](int reco_beam_type){
  return (reco_beam_type == 13);
};

//TPC values based Beam Cut
//Tingjun: https://indico.fnal.gov/event/49434/contributions/217088/attachments/144364/183297/backgroundstudies.pdf
auto beamQuality_mc_TPCinfo = [](double calo_beam_startX, double calo_beam_startY,
                                    double calo_beam_startZ, double calo_beam_endX,
                                    double calo_beam_endY, double calo_beam_endZ)   {

   double diffX = calo_beam_endX - calo_beam_startX;
   double diffY = calo_beam_endY - calo_beam_startY;
   double diffZ = calo_beam_endZ - calo_beam_startZ;

   double cosTrk_thetaX = diffX / sqrt( diffX*diffX + diffZ*diffZ );
   double cosTrk_thetaY = diffY / sqrt( diffY*diffY + diffZ*diffZ );
   double cosTrk_thetaZ = diffZ / sqrt( diffX*diffX + diffZ*diffZ );

   double cosTheta = cos(mc_thetaX) * cosTrk_thetaX + cos(mc_thetaY) * cosTrk_thetaY + cos(mc_thetaZ) * cosTrk_thetaZ;

   if( sqrt( pow( (calo_beam_startX - mc_meanX ) / mc_sigmaX, 2) + pow( (calo_beam_startY - mc_meanY ) / mc_sigmaY , 2) ) > cut_beamQuality_TPC_xyz )
      return false;

   if( abs( (calo_beam_startZ - mc_meanZ ) / mc_sigmaZ ) > cut_beamQuality_TPC_xyz )
      return false;

   if( cosTheta < cut_beamQuality_TPC_cosTheta )
      return false;

   return true;
};

auto beamQuality_mc_TPCjustPosition = [](double calo_beam_startX, double calo_beam_startY,
                                    double calo_beam_startZ)   {

   if( sqrt( pow( (calo_beam_startX - mc_meanX ) / mc_sigmaX, 2) + pow( (calo_beam_startY - mc_meanY ) / mc_sigmaY , 2) ) > cut_beamQuality_TPC_xyz )
      return false;

   if( abs( (calo_beam_startZ - mc_meanZ ) / mc_sigmaZ ) > cut_beamQuality_TPC_xyz )
      return false;

   return true;
};

auto beamQuality_data_TPCinfo = [](double calo_beam_startX, double calo_beam_startY,
                                    double calo_beam_startZ, double calo_beam_endX,
                                    double calo_beam_endY, double calo_beam_endZ)   {

   double diffX = calo_beam_endX - calo_beam_startX;
   double diffY = calo_beam_endY - calo_beam_startY;
   double diffZ = calo_beam_endZ - calo_beam_startZ;

   double cosTrk_thetaX = diffX / sqrt( diffX*diffX + diffZ*diffZ );
   double cosTrk_thetaY = diffY / sqrt( diffY*diffY + diffZ*diffZ );
   double cosTrk_thetaZ = diffZ / sqrt( diffX*diffX + diffZ*diffZ );

   double cosTheta = cos(data_thetaX) * cosTrk_thetaX + cos(data_thetaY) * cosTrk_thetaY + cos(data_thetaZ) * cosTrk_thetaZ;

   if( sqrt( pow( (calo_beam_startX - data_meanX ) / data_sigmaX, 2) + pow( (calo_beam_startY - data_meanY ) / data_sigmaY , 2) ) > cut_beamQuality_TPC_xyz )
      return false;

   if( abs( (calo_beam_startZ - data_meanZ ) / data_sigmaZ ) > cut_beamQuality_TPC_xyz )
      return false;

   if( cosTheta < cut_beamQuality_TPC_cosTheta )
      return false;

   return true;
};

auto beamQuality_data_TPCjustPosition = [](double calo_beam_startX, double calo_beam_startY,
                                    double calo_beam_startZ)   {

   if( sqrt( pow( (calo_beam_startX - data_meanX ) / data_sigmaX, 2) + pow( (calo_beam_startY - data_meanY ) / data_sigmaY , 2) ) > cut_beamQuality_TPC_xyz )
      return false;

   if( abs( (calo_beam_startZ - data_meanZ ) / data_sigmaZ ) > cut_beamQuality_TPC_xyz )
      return false;

   return true;
};

auto manual_beamPos_mc = [](double beam_startX, double beam_startY,
                            double beam_startZ, double beam_dirX,
                            double beam_dirY,   double beam_dirZ, 
                            double true_dirX,   double true_dirY,
                            double true_dirZ,   double true_startX,
                            double true_startY, double true_startZ) {
  double projectX = (true_startX + -1*true_startZ*(true_dirX/true_dirZ) );
  double projectY = (true_startY + -1*true_startZ*(true_dirY/true_dirZ) );
  double cos = true_dirX*beam_dirX + true_dirY*beam_dirY + true_dirZ*beam_dirZ;

  if ( (beam_startX - projectX) < xlow )
    return false;
  
  if ( (beam_startX - projectX) > xhigh )
    return false;

  if ( (beam_startY - projectY) < ylow )
    return false;

  if ( (beam_startY - projectY) > yhigh )
    return false;
  
  if (beam_startZ < zlow || zhigh < beam_startZ)
    return false;
  
  if ( cos < coslow)
    return false;

  return true;

};

auto data_beam_PID = [](const std::vector<int>& pidCandidates){
  auto pid_search = std::find(pidCandidates.begin(), pidCandidates.end(), 211);
  return (pid_search != pidCandidates.end());
};

auto data_BI_quality = [](int data_BI_nMomenta, int data_BI_nTracks) {
  return (data_BI_nMomenta == 1 && data_BI_nTracks == 1);
};

auto manual_beamPos_data = [](double data_startX,
                              double data_startY,   double data_startZ,
                              double data_dirX,     double data_dirY,
                              double data_dirZ,     double data_BI_X,
                              double data_BI_Y,     double data_BI_dirX,
                              double data_BI_dirY,  double data_BI_dirZ,
                              int data_BI_nMomenta, int data_BI_nTracks) {

  double deltaX = data_startX - data_BI_X;
  double deltaY = data_startY - data_BI_Y;
  double cos = data_BI_dirX*data_dirX + data_BI_dirY*data_dirY +
               data_BI_dirZ*data_dirZ;

  if(data_BI_nMomenta != 1 || data_BI_nTracks != 1)
    return false;

  if( (deltaX < data_xlow) || (deltaX > data_xhigh) )
    return false;

  if ( (deltaY < data_ylow) || (deltaY > data_yhigh) )
    return false;

  if ( (data_startZ < data_zlow) || (data_startZ > data_zhigh) )
    return false;

  if (cos < data_coslow)
    return false;

  return true;

};

auto beam_cut_MC_BI = [](double startX,
                         double startY,   double startZ,
                         double dirX,     double dirY,
                         double dirZ,     double BI_X,
                         double BI_Y,     double BI_dirX,
                         double BI_dirY,  double BI_dirZ,
                         int BI_nMomenta, int BI_nTracks) {

  double deltaX = startX - BI_X;
  double deltaY = startY - BI_Y;
  double cos = BI_dirX*dirX + BI_dirY*dirY +
               BI_dirZ*dirZ;
  if(BI_nMomenta != 1 || BI_nTracks != 1)
    return false;

  if( (deltaX < mc_BI_xlow) || (deltaX > mc_BI_xhigh) )
    return false;

  if ( (deltaY < mc_BI_ylow) || (deltaY > mc_BI_yhigh) )
    return false;

  if ( (startZ < mc_BI_zlow) || (startZ > mc_BI_zhigh) )
    return false;

  if (cos < mc_BI_coslow)
    return false;

  return true;
};


//for marking cutflow in rows only needs condition before and tested
auto cutFlow = [](bool a, bool b){
  return (a && b);
};

//Distance to Vertex Daughter
auto compute_distanceVertex = [](double beam_endX,
                                 double beam_endY,
                                 double beam_endZ, 
                                 const std::vector<double> &d_startX,
                                 const std::vector<double> &d_startY,
                                 const std::vector<double> &d_startZ,
                                 const std::vector<double> &d_endX,
                                 const std::vector<double> &d_endY,
                                 const std::vector<double> &d_endZ) {
  std::vector<double> distance;
  double dummy = 0., dummy_1 = 0., dummy_2 = 0.;
  double diff_X_end = 0., diff_Y_end = 0., diff_Z_end = 0.;
  double diff_X_start = 0., diff_Y_start = 0., diff_Z_start = 0.;

  if(d_startX.empty()) return distance;

  for( size_t i = 0; i < d_startX.size(); ++i ) {
    diff_X_end = d_endX[i] - beam_endX;
    diff_Y_end = d_endY[i] - beam_endY;
    diff_Z_end = d_endZ[i] - beam_endZ;

    diff_X_start = d_startX[i] - beam_endX;
    diff_Y_start = d_startY[i] - beam_endY;
    diff_Z_start = d_startZ[i] - beam_endZ;

    dummy_1 = sqrt(diff_X_end*diff_X_end + diff_Y_end*diff_Y_end + 
                   diff_Z_end*diff_Z_end);

    dummy_2 = sqrt(diff_X_start*diff_X_start + diff_Y_start*diff_Y_start +
                   diff_Z_start*diff_Z_start);

    if(dummy_1 < dummy_2)
      distance.push_back(dummy_1);
    else 
      distance.push_back(dummy_2);
  }

  return distance;
};



//Removing primaryMuons in the Incident Pion sample
//by looking for MichelElectron with CNN score

auto candidate_primaryMuon = []( double michelScore, int nhits){

   if(michelScore == -999.) return false; //there is no michel score available keep event

   if( michelScore/nhits > cut_michelScore) return true; //there is a michel candidate


   return false;
};



//Using dEdX truncated mean and trackscore
auto secondary_noPion= [](
                           const std::vector<double> &track_score, 
                           const std::vector<int> &trackID,
                           const std::vector<double> &dEdX,
                           const std::vector<double> &chi2,
                           const std::vector<int> &ndof) {
  for( size_t i = 0; i < track_score.size(); ++i ) {
    if ((trackID[i] != -999) && (track_score[i] > cut_trackScore)) {
      //if (dEdX[i] <= cut_dEdX)) {
      if (dEdX[i] < 2.8 && dEdX[i] > 0.5) {
        return false;
      }
      //else if (dEdX[i] > 2.8 && dEdX[i] < 3.4) {
      else if (dEdX[i] < 3.4) {
        if (ndof[i] > 0 && chi2[i]/ndof[i] > 70.) {
          return false;
        }
      }
    }
  }

  return true;
};

auto is_secondary_pion= [](
                           const std::vector<double> &track_score, 
                           const std::vector<int> &trackID,
                           const std::vector<double> &dEdX,
                           const std::vector<double> &chi2,
                           const std::vector<int> &ndof) {
  std::vector<bool> results;
  for( size_t i = 0; i < track_score.size(); ++i ) {
    if ((trackID[i] != -999) && (track_score[i] > cut_trackScore)) {
      //if (dEdX[i] <= cut_dEdX)) {
      if (dEdX[i] < 2.8 && dEdX[i] > 0.5) {
        results.push_back(true);
      }
      //else if (dEdX[i] > 2.8 && dEdX[i] < 3.4) {
      else if (dEdX[i] < 3.4) {
        if (ndof[i] > 0 && chi2[i]/ndof[i] > 70.) {
          results.push_back(true);
        }
      }
      else {
        results.push_back(false);
      }
    }
    else {
      results.push_back(false);
    }
  }
  return results;
};


auto has_shower_nHits = [](const std::vector<double> &track_score,
                           const std::vector<int> &nHits) {
  if(track_score.empty() || nHits.empty())
    return false;

  for(size_t i = 0; i < track_score.size(); ++i){
     if ((track_score[i] < cut_trackScore) &&
         (nHits[i] > cut_nHits_shower_low) &&
         (track_score[i] != -999.)) {
       return true;
     }
  }

  return false;
};

auto has_shower_dist_energy = [](const std::vector<double> &track_score,
                                 const std::vector<double> &shower_x,
                                 const std::vector<double> &shower_y,
                                 const std::vector<double> &shower_z,
                                 const std::vector<double> &energy,
                                 double & x, double & y, double & z) {
  for(size_t i = 0; i < track_score.size(); ++i){
     double dist = sqrt(std::pow((shower_x[i] - x), 2) +
                        std::pow((shower_y[i] - y), 2) +
                        std::pow((shower_z[i] - z), 2));
     if ((track_score[i] < cut_trackScore) &&
         (track_score[i] > 0.) &&
         (dist > 5. && dist < 1000.) &&
         (energy[i] > 80. && energy[i] < 1000.)) {
       return true;
     }
  }

  return false;
};

auto is_pi0_shower = [](const std::vector<double> &track_score,
                                 const std::vector<double> &shower_x,
                                 const std::vector<double> &shower_y,
                                 const std::vector<double> &shower_z,
                                 const std::vector<double> &energy,
                                 double & x, double & y, double & z) {
  std::vector<bool> results;
  for(size_t i = 0; i < track_score.size(); ++i){
    double dist = sqrt(std::pow((shower_x[i] - x), 2) +
                       std::pow((shower_y[i] - y), 2) +
                       std::pow((shower_z[i] - z), 2));
    if ((track_score[i] < cut_trackScore) &&
        (track_score[i] > 0.) &&
        (dist > 5. && dist < 1000.) &&
        (energy[i] > 80. && energy[i] < 1000.)) {
      results.push_back(true);
    }
    else {
     results.push_back(false);
    }
  }

  return results;
};

auto shower_dists = [](const std::vector<double> &track_score,
                       const std::vector<double> &shower_x,
                       const std::vector<double> &shower_y,
                       const std::vector<double> &shower_z,
                       double & x, double & y, double & z) {
  std::vector<double> results;
  for(size_t i = 0; i < track_score.size(); ++i){
    if ((track_score[i] < cut_trackScore) &&
        (track_score[i] > 0.)) {
      double dist = sqrt(std::pow((shower_x[i] - x), 2) +
                         std::pow((shower_y[i] - y), 2) +
                         std::pow((shower_z[i] - z), 2));
      results.push_back(dist);
    }
    else {
      results.push_back(-999.);
    }
  }

  return results;
};

auto leading_proton_momentum = [](const std::vector<double> & daughter_p,
                                  const std::vector<int> & daughter_pdg/*,
                                  double proton_threshold = 0.*/) {
  double max_p = -1.;
  for (size_t i = 0; i < daughter_pdg.size(); ++i) {
    if (daughter_pdg[i] == 2212) {
      if (daughter_p[i] < .2)
        continue;
      if (daughter_p[i] > max_p)
        max_p = daughter_p[i];
    }
  }

  return max_p;
};

auto leading_proton_det_theta = [](const std::vector<double> & daughter_p,
                                   const std::vector<int> & daughter_pdg,
                                   const std::vector<double> & daughter_pz) {
  double max_p = -1.;
  double max_theta = -999.; 
  for (size_t i = 0; i < daughter_pdg.size(); ++i) {
    if (daughter_pdg[i] == 2212) {
      if (daughter_p[i] < .2)
        continue;
      if (daughter_p[i] > max_p) {
        max_p = daughter_p[i];
        max_theta = daughter_pz[i]/daughter_p[i];
      }
    }
  }

  return max_theta;
};

auto n_track_daughters = [](const std::vector<double> &track_score, 
                            const std::vector<int> &trackID) {
  int results = 0;
  for (size_t i = 0; i < track_score.size(); ++i) {
    if ((trackID[i] != -999) && (track_score[i] > cut_trackScore)) {
      ++results; 
    }
  }
  return results;
};

auto n_shower_daughters = [](const std::vector<double> &track_score, 
                             const std::vector<int> &showerID) {
  int results = 0;
  for (size_t i = 0; i < track_score.size(); ++i) {
    if ((showerID[i] != -999) && (track_score[i] < cut_trackScore)) {
      ++results; 
    }
  }
  return results;
};

auto leading_proton_det_phi = [](const std::vector<double> & daughter_p,
                                   const std::vector<int> & daughter_pdg,
                                   const std::vector<double> & daughter_px,
                                   const std::vector<double> & daughter_py) {
  double max_p = -1.;
  double max_phi = -999.; 
  for (size_t i = 0; i < daughter_pdg.size(); ++i) {
    if (daughter_pdg[i] == 2212) {
      if (daughter_p[i] < .2)
        continue;
      if (daughter_p[i] > max_p) {
        max_p = daughter_p[i];
        max_phi = atan(daughter_px[i]/daughter_py[i]) * 180. / TMath::Pi();
      }
    }
  }

  return max_phi;
};

auto backtrack_beam = [](const std::string process,
                         const bool matched,
                         const int origin, const int PDG) {
  if (process == "primary" && matched && origin == 4 && PDG == 211) {
    return 1;
  }
  else if (process == "primary" && matched && origin == 4 && PDG == -13) {
    return 2;
  }
  else if (origin == 2) {
    return 3;
  }
  else if (process == "primaryBackground") {
    return 4;
  }
  else if (process.find("Inelastic") != std::string::npos) {
    return 5;
  }
  else if (process == "Decay") {
    return 6;
  }
  else {
    return 7;
  }
};

auto categorize_daughters = [](
    const int beam_ID,
    const std::vector<int> bt_origins, const std::vector<int> bt_IDs,
    const std::vector<int> bt_PDGs, const std::vector<int> bt_ParIDs,
    const std::vector<int> bt_ParPDGs,
    const std::vector<int> true_daughters,
    const std::vector<int> true_grand_daughters) {
  std::vector<int> results;
  //std::cout << "asdf" << std::endl;
  for (size_t i = 0; i < bt_origins.size(); ++i) {
    if (bt_IDs[i] == beam_ID) {
      results.push_back(1);
    }
    else if (bt_origins[i] == 2) {
      results.push_back(2);
    }
    else if (abs(bt_PDGs[i]) == 11 && (abs(bt_ParPDGs[i]) == 13)) {
      results.push_back(10); 
    }
    else if (std::find(true_daughters.begin(), true_daughters.end(), bt_IDs[i]) !=
             true_daughters.end()) {
      if (abs(bt_PDGs[i]) == 211) {
        results.push_back(3);
      }
      else if (abs(bt_PDGs[i]) == 13) {
        results.push_back(4);
      }
      else if (bt_PDGs[i] == 2212) {
        results.push_back(5);
      }
      else if (bt_PDGs[i] == 22) {
        results.push_back(6);
      }
      else if (bt_PDGs[i] > 2212) {
        results.push_back(7);
      }
      else {
        //std::cout << bt_PDGs[i] << " " << bt_ParPDGs[i] << " " << (abs(bt_PDGs[i]) == 11 && (abs(bt_ParPDGs[i]) == 13)) << std::endl;
        results.push_back(12);
      }
    }
    else if (std::find(true_grand_daughters.begin(), 
                       true_grand_daughters.end(), bt_IDs[i]) !=
             true_grand_daughters.end()) {
      if ((bt_PDGs[i] == 22 || abs(bt_PDGs[i]) == 11) && bt_ParPDGs[i] == 111) {
        results.push_back(11);
      }
      else {
        results.push_back(8);
      }
    }
    else if (std::find(true_grand_daughters.begin(), 
                       true_grand_daughters.end(), bt_ParIDs[i]) !=
             true_grand_daughters.end()) {
        results.push_back(9);
    }
    else {
        //std::cout << bt_PDGs[i] << " " << bt_ParPDGs[i] << " " << (abs(bt_PDGs[i]) == 11 && (abs(bt_ParPDGs[i]) == 13)) << std::endl;
        results.push_back(12);
    }
  }
  return results;
};

auto daughter_PDG_types(const std::vector<int> bt_PDGs) {
  std::vector<int> results;
  for (const int PDG : bt_PDGs) {
    if (abs(PDG) == 211) {
      results.push_back(1);
    }
    else if (abs(PDG) == 13) {
      results.push_back(2);
    }
    else if (abs(PDG) == 2212) {
      results.push_back(3);
    }
    else if (abs(PDG) == 22) {
      results.push_back(4);
    }
    else if (abs(PDG) > 2212) {
      results.push_back(5);
    }
    else if (abs(PDG) == 11) {
      results.push_back(6);
    }
    else {
      results.push_back(7);
    }
  }
  return results;
};
