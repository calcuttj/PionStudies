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

const Double_t cutdummy= -99999.0;
const string path = "eSliceMethod_Prod4a_mc_1GeV_all_08_02_21.root";
const string dataPath = "eSliceMethod_Prod4a_58XX_1GeV_all_08_02_21.root";
const Double_t dataMC_events = 68500; //root file DATA has 70121 ev, root file MC has 74149 ev


//==============================================================================
// Example Unfolding
//==============================================================================

void unfold_bayes_wBG(const string mcFilepath = path, bool doMC = true, bool doXS = true)
{
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   ROOT::RDataFrame data_inputFrame(pionTree, dataPath);

   gStyle->SetNdivisions(1020);
   //ROOT::RDataFrame data_pre(pionTree, dataPath);

   string output_name;
   //if(doMC) output_name = "unfold_wBG_mc_" + std::to_string((int) bin_size_int) + "MeV.root";
   //else 
   output_name = "ratio_unfold_" + std::to_string((int) bin_size_int) + "MeV.root";

   TFile *output = new TFile( output_name.c_str() , "RECREATE"); //maybe save with binning?

   //=======================  Frame Definitions ==============================================================      

   //no need to filter for only reconstructed events as with the Miss function one can take into account the not-reconstructed events
   auto frame = inputFrame
      .Filter("pass_trueBeamLocation")// bool that marks in true the particles that passed the true beam location, same range in x and y as for the TPCposition cut for reco but on true var
      //.Filter("primary_isBeamType")
      .Range(dataMC_events)
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
            " && ( !( true_absSignal && true_pion_daughter == 0) || true_equalBin )");

   auto dataFrame = data_inputFrame
      .Range(dataMC_events)
      .Define("reco_equalBin", equalBin, {"reco_initKE", "reco_interKE"});
   //.Filter("reco_initKE < 1200");

   //auto count_data = dataFrame.Count();
   //std::cout << "count data = " << *count_data << std::endl;

   //=====================================================================================      

   //Frame to Train response
   auto frame_train = frame
      .Filter("true_beam_endZ > 0");
   //.Filter("true_beam_endZ > 0 && true_initKE < 1200 ");

   //Frame for Measurement 
   auto frame_meas = frame
      .Filter("true_beam_endZ > 0");
   //.Filter("true_beam_endZ > 0 && true_initKE < 1200");


   cout << "==================================== TRAIN ====================================" << endl;
   //events that have reco_initBin == reco_interBin go into the Miss function as they are treated like a reco-ineff
   //events with true_initBin == true_interBin havre been filtered out already
   RooUnfoldResponse response_interE = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_interE", "");
   RooUnfoldResponse response_initE = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_initE", "");
   RooUnfoldResponse response_totInel = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_totInel", "");
   RooUnfoldResponse response_abs = RooUnfoldResponse(nBin_int, eEnd, eStart, 
         "response_abs", "");

   //If there is need to use Overflow do response_bla.UseOverflow(true)

   frame_train
      //.Range(half_mc)
      //.Range(53000) 
      .Foreach( [ &response_initE, &response_interE](double true_initE, double true_interE, 
               double reco_initE, double reco_interE, 
               bool ev_response, bool ev_miss, bool ev_fake){

            //Signal
            if(ev_response){
            response_initE.Fill( reco_initE, true_initE);
            response_interE.Fill( reco_interE, true_interE);
            }
            //MisReconstruction
            else if( ev_miss ) {
            response_initE.Miss (true_initE);
            response_interE.Miss (true_interE);
            }
            //BackGround will be subtracted before unfolding
            else if( ev_fake ){

            response_initE.Fake(reco_initE);
            response_interE.Fake(reco_interE);
            }

      }
   ,{"true_initKE", "true_interKE", "reco_initKE_rwData", "reco_interKE_rwData", 
      "truePion_response", "truePion_miss", "truePion_fake"});

   //Total Inelastic
   frame_train
      //.Range(half_mc)
      //.Range(53000) 
      .Foreach( [ &response_totInel](double true_initE, double true_interE, 
               double reco_initE, double reco_interE, 
               bool ev_response, bool ev_miss, bool ev_fake){

            //Signal
            if(ev_response){
            response_totInel.Fill(reco_interE, true_interE);
            }
            //MisReconstruction
            else if( ev_miss ) {
            response_totInel.Miss (true_interE);
            }
            //BackGround will be subtracted before unfolding
            else if( ev_fake ){
            response_totInel.Fake(reco_interE);
            }

            }
            ,{"true_initKE", "true_interKE", "reco_initKE_rwData", "reco_interKE_rwData", 
            "truePion_totInel_response", "truePion_totInel_miss", "truePion_totInel_fake"});

   //Total Inelastic
   frame_train
      //.Range(half_mc)
      //.Range(53000)
      .Foreach( [ &response_abs](double true_initE, double true_interE, 
               double reco_initE, double reco_interE, 
               bool ev_response, bool ev_miss, bool ev_fake){

            //Signal
            if(ev_response){
            response_abs.Fill(reco_interE, true_interE);
            }
            //MisReconstruction
            else if( ev_miss ) {
            response_abs.Miss (true_interE);
            }
            //BackGround will be subtracted before unfolding
            else if( ev_fake ){
            response_abs.Fake(reco_interE);
            }

            }
            ,{"true_initKE", "true_interKE", "reco_initKE_rwData", "reco_interKE_rwData", 
            "truePion_abs_response", "truePion_abs_miss", "truePion_abs_fake"});

   cout << "==================================== TEST =====================================" << endl;


   TH1D* hTrue_initE= new TH1D ("trueMC_initE", "MC Truth initE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hMeas_mc_initE= new TH1D ("meas_mc_initE", "Measured MC initE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_initE= new TH1D ("meas_data_initE", "Measured data initE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   TH1D* hTrue_interE= new TH1D ("trueMC_interE", "MC Truth interE; Energy [MeV]; ev/bin",    nBin_int, eEnd, eStart);
   TH1D* hMeas_mc_interE= new TH1D ("meas_mc_interE", "Measured MC interE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_interE= new TH1D ("meas_data_interE", "Measured data interE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   TH1D* hTrue_totInel= new TH1D ("trueMC_totInel", "MC Truth totInel; Energy [MeV]; ev/bin",    nBin_int, eEnd, eStart);
   TH1D* hMeas_mc_totInel= new TH1D ("meas_mc_totInel", "Measured MC totInel; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_totInel= new TH1D ("meas_data_totInel", "Measured data totInel; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   TH1D* hTrue_abs= new TH1D ("trueMC_abs", "MC Truth abs; Energy [MeV]; ev/bin",    nBin_int, eEnd, eStart);
   TH1D* hMeas_mc_abs= new TH1D ("meas_mc_abs", "Measured MC abs; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_abs= new TH1D ("meas_data_abs", "Measured data abs; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   //Build hMeas from same prefiltered sample that was supplied to response and fake (the Miss are misReco so won't catch them in reco)
   //Build hTrue from what was part of response and Miss (the fake are not part of hTrue

   cout << "==================================== Building the Measured Distribution =====================================" << endl;

   frame_meas

      .Filter("selected_incidentPion")
      .Foreach( [ hMeas_mc_initE, hMeas_mc_interE, hMeas_mc_totInel ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){

            hMeas_mc_initE->Fill( reco_initE );
            hMeas_mc_interE->Fill( reco_interE );
            hMeas_mc_totInel->Fill( reco_interE );
            }

            /*for( int i = 1; i <= nBin_int; i++){
              hMeas_mc_initE->SetBinError( i , hMeas_mc_initE->GetBinContent(i)*0.05 );
              hMeas_mc_interE->SetBinError( i , hMeas_mc_interE->GetBinContent(i)*0.05 );
              hMeas_mc_totInel->SetBinError( i , hMeas_mc_totInel->GetBinContent(i)*0.05 );

              }*/

            }
            ,{"reco_initKE_rwData", "reco_interKE_rwData", "reco_equalBin"}
            //,{"reco_initKE", "reco_interKE", "reco_equalBin"}
            );

   frame_meas

      .Filter("selected_abs")
      .Foreach( [ hMeas_mc_abs ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){
            hMeas_mc_abs->Fill( reco_interE );
            }

            /*for(int i = 1; i <= nBin_int; i++){
              hMeas_mc_abs->SetBinError( i , hMeas_mc_abs->GetBinContent(i)*0.05 );

              }*/
            }
            ,{"reco_initKE_rwData", "reco_interKE_rwData", "reco_equalBin"}
            //,{"reco_initKE", "reco_interKE", "reco_equalBin"}
            );

   //====================  DATA ==========================
   dataFrame
      .Filter("selected_incidentPion")
      .Foreach( [ hMeas_data_initE, hMeas_data_interE, hMeas_data_totInel ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){

            hMeas_data_initE->Fill( reco_initE );
            hMeas_data_interE->Fill( reco_interE );
            hMeas_data_totInel->Fill( reco_interE );
            }

             /*
            for(int i = 1; i <= nBin_int; i++){
            hMeas_data_initE->SetBinError( i , hMeas_data_initE->GetBinContent(i)*0.2 );
            hMeas_data_interE->SetBinError( i , hMeas_data_interE->GetBinContent(i)*0.2 );
            hMeas_data_totInel->SetBinError( i , hMeas_data_totInel->GetBinContent(i)*0.2 );

            }
            // */
            }
            ,{"reco_initKE", "reco_interKE", "reco_equalBin"}
            );

   dataFrame
      .Filter("selected_abs")
      .Foreach( [ hMeas_data_abs ]( double reco_initE, double reco_interE, bool recoEqualBin){

            //Measured
            if(reco_initE != -999. && reco_interE != -999. && !recoEqualBin){
            hMeas_data_abs->Fill( reco_interE );
            }

             /*
            for(int i = 1; i <= nBin_int; i++){
            hMeas_data_abs->SetBinError( i , hMeas_data_abs->GetBinContent(i)*0.2 );

            }
            // */
            }
            ,{ "reco_initKE", "reco_interKE", "reco_equalBin"}
            );

   hMeas_mc_initE->Write();
   hMeas_mc_interE->Write();
   hMeas_mc_totInel->Write();
   hMeas_mc_abs->Write();

   hMeas_data_initE->Write();
   hMeas_data_interE->Write();
   hMeas_data_totInel->Write();
   hMeas_data_abs->Write();


   cout << "==================================== Building the True Distribution =====================================" << endl;


   frame_meas

      .Foreach( [ hTrue_initE, hTrue_interE ]( double true_initE, double true_interE, 
               bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){

            hTrue_initE->Fill( true_initE );
            hTrue_interE->Fill( true_interE );
            }

            }
            ,{"true_initKE", "true_interKE", "truePion_response", "truePion_miss"}
            );
   frame_meas

      .Foreach( [ hTrue_totInel ]( double true_initE, double true_interE, 
               bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){
            hTrue_totInel->Fill( true_interE );
            }

            }
            ,{ "true_initKE", "true_interKE", "truePion_totInel_response", "truePion_totInel_miss"}
            );

   frame_meas

      .Foreach( [ hTrue_abs ]( double true_initE, double true_interE, 
               bool ev_response, bool ev_miss){

            if( ev_response || ev_miss ){
            hTrue_abs->Fill( true_interE );
            }

            }
            ,{ "true_initKE", "true_interKE", "truePion_abs_response", "truePion_abs_miss"}
            );

   hTrue_initE->Write();
   hTrue_interE->Write();
   hTrue_totInel->Write();
   hTrue_abs->Write();


   RooUnfoldBayes unfold_initE;
   RooUnfoldBayes unfold_interE;
   RooUnfoldBayes unfold_totInel;
   RooUnfoldBayes unfold_abs;

   if(doMC){
      cout << "==================================== UNFOLD MC===================================" << endl;
      unfold_initE = RooUnfoldBayes(&response_initE, hMeas_mc_initE, 8);    // OR
      unfold_interE = RooUnfoldBayes(&response_interE, hMeas_mc_interE, 8);    // OR
      unfold_totInel = RooUnfoldBayes(&response_totInel, hMeas_mc_totInel, 8);    // OR
      unfold_abs = RooUnfoldBayes(&response_abs, hMeas_mc_abs, 8);    // OR
   }

   else{
      cout << "==================================== UNFOLD MC ===================================" << endl;
      unfold_initE = RooUnfoldBayes(&response_initE, hMeas_data_initE, 15);    // OR
      unfold_interE = RooUnfoldBayes(&response_interE, hMeas_data_interE, 15);    // OR
      unfold_totInel = RooUnfoldBayes(&response_totInel, hMeas_data_totInel, 15);    // OR
      unfold_abs = RooUnfoldBayes(&response_abs, hMeas_data_abs, 15);    // OR
   }


   cout << "=========================== Covariance Matrices ===================================" << endl;

   TMatrixD cov_initE = unfold_initE.Ereco();
   //unfold_initE.SetMeasuredCov( cov_initE );

   TMatrixD cov_interE = unfold_interE.Ereco();
   //unfold_interE.SetMeasuredCov( cov_interE );

   TMatrixD cov_totInel = unfold_totInel.Ereco();
   //unfold_totInel.SetMeasuredCov( cov_totInel );

   TMatrixD cov_abs = unfold_abs.Ereco();
   //unfold_abs.SetMeasuredCov( cov_abs ); 

   cout << "------------------------------------ Init E -----------------------------------" << endl;
   //TH1D* hUnfold_initE = (TH1D*) unfold_initE.Hreco((RooUnfold::ErrorTreatment)RooUnfold::kCovariance);
   TH1D* hUnfold_initE = (TH1D*) unfold_initE.Hreco();
   cout << "------------------------------------ Inter E -----------------------------------" << endl;
   //TH1D* hUnfold_interE = (TH1D*) unfold_interE.Hreco((RooUnfold::ErrorTreatment)RooUnfold::kCovariance);
   TH1D* hUnfold_interE = (TH1D*) unfold_interE.Hreco();
   cout << "------------------------------------ TotInel -----------------------------------" << endl;
   //TH1D* hUnfold_totInel = (TH1D*) unfold_totInel.Hreco((RooUnfold::ErrorTreatment)RooUnfold::kCovariance);
   TH1D* hUnfold_totInel = (TH1D*) unfold_totInel.Hreco();
   cout << "------------------------------------ abs -----------------------------------" << endl;
   //TH1D* hUnfold_abs = (TH1D*) unfold_abs.Hreco((RooUnfold::ErrorTreatment)RooUnfold::kCovariance);
   TH1D* hUnfold_abs = (TH1D*) unfold_abs.Hreco();

   

   if(doMC){
      hUnfold_initE->SetNameTitle("hUnfold_initE", "Unfold MC initE; Energy [MeV]; ev/bin");
      hUnfold_interE->SetNameTitle("hUnfold_interE", "Unfold MC interE; Energy [MeV]; ev/bin");
      hUnfold_totInel->SetNameTitle("hUnfold_totInel", "Unfold MC totInel; Energy [MeV]; ev/bin");
      hUnfold_abs->SetNameTitle("hUnfold_abs", "Unfold MC abs; Energy [MeV]; ev/bin");
   }
   else{
      hUnfold_initE->SetNameTitle("hUnfold_initE", "Unfold Data initE; Energy [MeV]; ev/bin");
      hUnfold_interE->SetNameTitle("hUnfold_interE", "Unfold Data interE; Energy [MeV]; ev/bin");
      hUnfold_totInel->SetNameTitle("hUnfold_totInel", "Unfold Data totInel; Energy [MeV]; ev/bin");
      hUnfold_abs->SetNameTitle("hUnfold_abs", "Unfold Data abs; Energy [MeV]; ev/bin");
   }
   hUnfold_initE->Write();
   hUnfold_interE->Write();
   hUnfold_totInel->Write();
   hUnfold_abs->Write();
   //print infos with PrintTable command
   cout << "==================================== UNFOLD INIT E===================================" << endl;
   unfold_initE.PrintTable (cout, hTrue_initE);
   cout << "==================================== UNFOLD INTER E===================================" << endl;
   unfold_interE.PrintTable (cout, hTrue_interE);
   cout << "==================================== UNFOLD TOTINEL E===================================" << endl;
   unfold_totInel.PrintTable (cout, hTrue_totInel);
   cout << "==================================== UNFOLD abs E===================================" << endl;
   unfold_abs.PrintTable (cout, hTrue_abs);

   //cout << "=========================== SCALE MC MEas TO DATA Meas by AREA AFTER having Added MUON CONTENT===================================" << endl;
   //
   //SCALING IS STILL AN ISSUE .... ! this needs to be adressed and properly done



   double scale_initE = hMeas_data_initE->Integral() / hMeas_mc_initE->Integral();
   //double scale_interE = hMeas_data_interE->Integral() / hMeas_mc_interE->Integral();
   //double scale_totInel = hMeas_data_totInel->Integral() / hMeas_mc_totInel->Integral();
   //double scale_abs = hMeas_data_abs->Integral() / hMeas_mc_abs->Integral();

   if(!doMC){
      hMeas_mc_initE->Scale( scale_initE );
      hMeas_mc_interE->Scale( scale_initE );
      hMeas_mc_totInel->Scale( scale_initE );
      hMeas_mc_abs->Scale( scale_initE );

      hTrue_initE->Scale( scale_initE );
      hTrue_interE->Scale( scale_initE );
      hTrue_totInel->Scale( scale_initE );
      hTrue_abs->Scale( scale_initE );
   }
   // hMeas_mc_interE->Scale( scale_interE );
   // hMeas_mc_totInel->Scale( scale_totInel );
   // hMeas_mc_abs->Scale( scale_abs );

   // hTrue_initE->Scale( scale_initE );
   // hTrue_interE->Scale( scale_interE );
   // hTrue_totInel->Scale( scale_totInel );
   // hTrue_abs->Scale( scale_abs );


   cout << "=========================== CHECK INTEGRALS INITE / INTER E===================================" << endl;
   cout << "Integral of Measured MC initE = " << hMeas_mc_initE->Integral() << endl;
   cout << "Integral of Measured MC interE = " << hMeas_mc_interE->Integral() << endl;
   cout << "Integral of Measured Data initE = " << hMeas_data_initE->Integral() << endl;
   cout << "Integral of Measured Data interE = " << hMeas_data_interE->Integral() << endl;
   cout << "Integral of True initE = " << hTrue_initE->Integral() << endl;
   cout << "Integral of True interE = " << hTrue_interE->Integral() << endl;
   cout << "Integral of Unfold initE = " << hUnfold_initE->Integral() << endl;
   cout << "Integral of Unfold interE = " << hUnfold_interE->Integral() << endl;

   cout << "=========================== CHECK INTEGRALS INITE / INTER E===================================" << endl;
   cout << "Integral of Measured totInel = " << hMeas_mc_totInel->Integral() << endl;
   cout << "Integral of Measured abs = " << hMeas_mc_abs->Integral() << endl;
   cout << "Integral of True totInel = " << hTrue_totInel->Integral() << endl;
   cout << "Integral of True abs = " << hTrue_abs->Integral() << endl;
   cout << "Integral of Unfold totInel = " << hUnfold_totInel->Integral() << endl;
   cout << "Integral of Unfold abs = " << hUnfold_abs->Integral() << endl;

   cout << "=========================== Correlation Matrices ===================================" << endl;
   cout << "------------------------------------ init E -----------------------------------" << endl;

   Int_t nb_initE = cov_initE.GetNrows();
   TH2D* corr_initE = new TH2D("corr_initE", "Correlations Init E", nb_initE, eEnd, eStart, nb_initE, eEnd, eStart);
   correlationMatrix( cov_initE, corr_initE, nb_initE);
   corr_initE->Write();

   cout << "------------------------------------ inter E -----------------------------------" << endl;

   Int_t nb_interE = cov_interE.GetNrows();
   TH2D* corr_interE = new TH2D("corr_interE", "Correlations Inter E", nb_interE, eEnd, eStart, nb_interE, eEnd, eStart);
   correlationMatrix( cov_interE, corr_interE, nb_interE);
   corr_interE->Write();

   cout << "------------------------------------ TotInel -----------------------------------" << endl;

   Int_t nb_totInel = cov_totInel.GetNrows();
   TH2D* corr_totInel = new TH2D("corr_totInel", "Correlations TotInel E", nb_totInel, eEnd, eStart, nb_totInel, eEnd, eStart);
   correlationMatrix( cov_totInel, corr_totInel, nb_totInel);
   corr_totInel->Write();

   cout << "------------------------------------ Absorption -----------------------------------" << endl;

   Int_t nb_abs = cov_abs.GetNrows();
   TH2D* corr_abs = new TH2D("corr_abs", "Correlations Abs E", nb_abs, eEnd, eStart, nb_abs, eEnd, eStart);
   correlationMatrix( cov_abs, corr_abs, nb_abs);
   corr_abs->Write();

   cout << "=========================== Pulls  ===================================" << endl;
   TH1D* hPull_initE= new TH1D ("hPull_initE", "Pulls initE; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hPull_interE= new TH1D ("hPull_interE", "Pulls interE; Energy [MeV]; ev/bin",    nBin_int, eEnd, eStart);
   TH1D* hPull_totInel= new TH1D ("hPull_totInel", "Pulls totInel; Energy [MeV]; ev/bin",    nBin_int, eEnd, eStart);
   TH1D* hPull_abs= new TH1D ("hPull_abs", "Pulls abs; Energy [MeV]; ev/bin",    nBin_int, eEnd, eStart);

   pullHisto( hUnfold_initE, hTrue_initE, hPull_initE);
   pullHisto( hUnfold_interE, hTrue_interE, hPull_interE);
   pullHisto( hUnfold_totInel, hTrue_totInel, hPull_totInel);
   pullHisto( hUnfold_abs, hTrue_abs, hPull_abs);

   hPull_initE->Write(); hPull_interE->Write(); hPull_totInel->Write(); hPull_abs->Write();


   //==============================================================================================      
   TCanvas* c_initE= new TCanvas("canvas_initE","canvas_initE");
   gPad->SetGrid(1,1);
   hUnfold_initE->Draw();
   hMeas_mc_initE->SetFillColorAlpha(kBlue,0.1);
   hMeas_mc_initE->SetLineColor(kBlue);
   hMeas_mc_initE->Draw("SAME HIST");
   hMeas_data_initE->SetFillColorAlpha(kBlack,0.3);
   hMeas_data_initE->Draw("SAME HIST");
   hTrue_initE->SetLineColor(8);
   hTrue_initE->Draw("HIST SAME");
   c_initE->BuildLegend();
   c_initE->Write();
   //c_initE->SaveAs("bla.pdf");

   TCanvas* c_interE= new TCanvas("canvas_interE","canvas_interE");
   gPad->SetGrid(1,1);
   hUnfold_interE->Draw();
   hMeas_mc_interE->SetFillColorAlpha(kBlue,0.1);
   hMeas_mc_interE->SetLineColor(kBlue);
   hMeas_mc_interE->Draw("SAME HIST");
   hMeas_data_interE->SetFillColorAlpha(kBlack,0.3);
   hMeas_data_interE->Draw("SAME HIST");
   hTrue_interE->SetLineColor(8);
   hTrue_interE->Draw("SAME HIST");
   c_interE->BuildLegend();
   c_interE->Write();

   TCanvas* c_totInel= new TCanvas("canvas_totInel","canvas_totInel");
   gPad->SetGrid(1,1);
   hUnfold_totInel->Draw();
   hMeas_mc_totInel->SetFillColorAlpha(kBlue,0.1);
   hMeas_mc_totInel->SetLineColor(kBlue);
   hMeas_mc_totInel->Draw("SAME HIST");
   hMeas_data_totInel->SetFillColorAlpha(kBlack,0.3);
   hMeas_data_totInel->Draw("SAME HIST");
   hTrue_totInel->SetLineColor(8);
   hTrue_totInel->Draw("SAME HIST");
   c_totInel->BuildLegend();
   c_totInel->Write();

   TCanvas* c_abs= new TCanvas("canvas_abs","canvas_abs");
   gPad->SetGrid(1,1);
   hUnfold_abs->Draw();
   hMeas_mc_abs->SetFillColorAlpha(kBlue,0.1);
   hMeas_mc_abs->SetLineColor(kBlue);
   hMeas_mc_abs->Draw("SAME HIST");
   hMeas_data_abs->SetFillColorAlpha(kBlack,0.3);
   hMeas_data_abs->Draw("SAME HIST");
   hTrue_abs->SetLineColor(8);
   hTrue_abs->Draw("SAME HIST");
   c_abs->BuildLegend();
   c_abs->Write();

   //Build the Incident histos in order to compare

   TH1D* hTrue_incident= new TH1D ("trueMC_incident", "MC Truth incident; Energy [MeV]; ev/bin",    nBin_int, eEnd, eStart);
   TH1D* hMeas_mc_incident= new TH1D ("meas_mc_incident", "Measured MC incident; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hMeas_data_incident= new TH1D ("meas_data_incident", "Measured Data incident; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);
   TH1D* hUnfold_incident= new TH1D ("hUnfold_incident", "Unfolded incident; Energy [MeV]; ev/bin", nBin_int, eEnd, eStart);

   build_incidentHist(hTrue_initE, hTrue_interE, hTrue_incident);
   build_incidentHist(hMeas_mc_initE, hMeas_mc_interE, hMeas_mc_incident);
   build_incidentHist(hUnfold_initE, hUnfold_interE, hUnfold_incident);
   build_incidentHist(hMeas_data_initE, hMeas_data_interE, hMeas_data_incident);


   //==============================================================================================      
   //build_incidentHist(hTrue_initE, hTrue_interE, hTrue_incident);

   //USE TRUE-INITE to build INCIDENT E histo for MEAS and UNFOLD
   //hTrue_initE->Scale( hMeas_mc_interE->Integral() / hTrue_initE->Integral() );
   //build_incidentHist(hTrue_initE, hMeas_mc_interE, hMeas_mc_incident);


   //hTrue_initE->Scale( hUnfold_interE->Integral() / hTrue_initE->Integral() );
   //build_incidentHist(hTrue_initE, hUnfold_interE, hUnfold_incident);

   TCanvas* c_incident= new TCanvas("canvas_incident","canvas_incident");
   gPad->SetGrid(1,1);
   hUnfold_incident->Draw();
   hMeas_mc_incident->SetFillColorAlpha(kBlue,0.1);
   hMeas_mc_incident->SetLineColor(kBlue);
   hMeas_mc_incident->Draw("SAME HIST");
   hMeas_data_incident->SetFillColorAlpha(kBlack,0.3);
   hMeas_data_incident->Draw("SAME HIST");
   hTrue_incident->SetLineColor(8);
   hTrue_incident->Draw("SAME HIST");
   //h_mcReco_incident->SetLineColor(38);
   //h_mcReco_incident->Draw("SAME HIST");
   c_incident->BuildLegend();
   c_incident->Write();

   auto* R_initE = response_initE.HresponseNoOverflow();
   auto* c_response_initE = new TCanvas();
   R_initE->SetStats(0);
   R_initE->Draw("colz");
   R_initE->GetXaxis()->SetNdivisions(1020);
   R_initE->GetYaxis()->SetNdivisions(1020);
   R_initE->SetTitle("Response Matrix Pion Initial Energy Distribution; reco Energy [MeV]; true Energy [MeV]");
   gPad->SetGrid(1,1);
   c_response_initE->Draw();

   auto* R_interE = response_interE.HresponseNoOverflow();
   auto* c_response_interE = new TCanvas();
   R_interE->SetStats(0);
   R_interE->Draw("colz");
   R_interE->GetXaxis()->SetNdivisions(1020);
   R_interE->GetYaxis()->SetNdivisions(1020);
   R_interE->SetTitle("Response Matrix Pion Interacting Energy Distribution; reco Energy [MeV]; true Energy [MeV]");
   gPad->SetGrid(1,1);
   c_response_interE->Draw();

   auto* R_totInel = response_totInel.HresponseNoOverflow();
   auto* c_response_totInel = new TCanvas();
   R_totInel->SetStats(0);
   R_totInel->Draw("colz");
   R_totInel->GetXaxis()->SetNdivisions(1020);
   R_totInel->GetYaxis()->SetNdivisions(1020);
   R_totInel->SetTitle("Response Matrix Pion Total Inelastic Energy Distribution; reco Energy [MeV]; true Energy [MeV]");
   gPad->SetGrid(1,1);
   c_response_totInel->Draw();

   auto* R_abs = response_abs.HresponseNoOverflow();
   auto* c_response_abs = new TCanvas();
   R_abs->SetStats(0);
   R_abs->Draw("colz");
   R_abs->GetXaxis()->SetNdivisions(1020);
   R_abs->GetYaxis()->SetNdivisions(1020);
   R_abs->SetTitle("Response Matrix Pion Absorption Energy Distribution; reco Energy [MeV]; true Energy [MeV]");
   gPad->SetGrid(1,1);
   c_response_abs->Draw();

   if(doXS){
      //access Jakes GeantFile in folder
      TFile f1("exclusive_xsec.root");
      TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
      TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
      TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
      f1.Close();

      output->cd();

      TH1D* h_betheMean_muon = new TH1D("h_betheMean_muon", "Mean Energy Loss", nBin_int, eEnd, eStart);

      //fill histo with Mean dEdX of bin center
      for(int i = 1; i <= nBin_int; i++){
         h_betheMean_muon->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size_int  , mass_muon) );
      };


      TH1D* h_unfoldXS_totInel = new TH1D("h_unfoldXS_totInel" ,"DATA Unfold Total Inelastic XS; Energy [MeV]; #sigma [mb]", nBin_int, eEnd, eStart);

      do_XS_log(  h_unfoldXS_totInel, hUnfold_totInel, hUnfold_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_unfoldXS_totInel, hUnfold_totInel, hUnfold_incident, h_betheMean_muon );   

      TH1D* h_TrueXS_totInel = new TH1D("h_TrueXS_totInel" ,"True Total Inelastic XS; Energy [MeV]; #sigma [mb]", nBin_int, eEnd, eStart);

      do_XS_log(  h_TrueXS_totInel, hTrue_totInel, hTrue_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_TrueXS_totInel, hTrue_totInel, hTrue_incident, h_betheMean_muon );   

      TH1D* h_rawXS_totInel = new TH1D("h_rawXS_totInel" ,"DATA Raw totInel XS; Energy [MeV]; #sigma [mb]", nBin_int, eEnd, eStart);

      do_XS_log(  h_rawXS_totInel, hMeas_data_totInel, hMeas_data_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_rawXS_totInel, hMeas_data_totInel, hMeas_data_incident, h_betheMean_muon );   


      TCanvas *c_unfold_totInel = new TCanvas("c_unfold_totInel", "c_unfold_totInel");
      gPad->SetGrid(1,1);
      h_unfoldXS_totInel->GetXaxis()->SetRangeUser(400,900);
      h_unfoldXS_totInel->GetXaxis()->SetNdivisions(1020);
      h_unfoldXS_totInel->GetYaxis()->SetNdivisions(1020);

      h_TrueXS_totInel->GetXaxis()->SetRangeUser(300,1000);

      totInel_KE->SetTitle( "Geant Total Inelasitc XS; Energy [MeV]; #sigma [mb]");
      totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
      totInel_KE->SetLineColor(kRed);
      totInel_KE->SetLineWidth(3);
      totInel_KE->Draw("AC");
      h_unfoldXS_totInel->SetMarkerSize(0.7);
      h_unfoldXS_totInel->Draw("PE0 SAME");
      h_TrueXS_totInel->SetMarkerColor(8);
      h_TrueXS_totInel->SetMarkerStyle(22);
      h_TrueXS_totInel->Draw("PE0 SAME");
      h_rawXS_totInel->SetMarkerColorAlpha(kOrange, 0.3);
      h_rawXS_totInel->Draw("P SAME");
      c_unfold_totInel->BuildLegend();

      c_unfold_totInel->Write();

      TH1D* h_unfoldXS_abs = new TH1D("h_unfoldXS_abs" ,"DATA Unfold Absorption XS; Energy [MeV]; #sigma [mb]", nBin_int, eEnd, eStart);

      do_XS_log(  h_unfoldXS_abs, hUnfold_abs, hUnfold_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_unfoldXS_abs, hUnfold_abs, hUnfold_incident, h_betheMean_muon );   

      TH1D* h_TrueXS_abs = new TH1D("h_TrueXS_abs" ,"True Absorption XS; Energy [MeV]; #sigma [mb]", nBin_int, eEnd, eStart);

      do_XS_log(  h_TrueXS_abs, hTrue_abs, hTrue_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_TrueXS_abs, hTrue_abs, hTrue_incident, h_betheMean_muon );   

      TH1D* h_rawXS_abs = new TH1D("h_rawXS_abs" ,"DATA Raw abs XS; Energy [MeV]; #sigma [mb]", nBin_int, eEnd, eStart);

      do_XS_log(  h_rawXS_abs, hMeas_data_abs, hMeas_data_incident, h_betheMean_muon );   
      do_XS_log_binomial_error( h_rawXS_abs, hMeas_data_abs, hMeas_data_incident, h_betheMean_muon );   

      TCanvas *c_unfold_abs = new TCanvas("c_unfold_abs", "c_unfold_abs");
      gPad->SetGrid(1,1);
      h_unfoldXS_abs->GetXaxis()->SetRangeUser(400,900);
      h_unfoldXS_abs->GetXaxis()->SetNdivisions(1020);
      h_unfoldXS_abs->GetYaxis()->SetNdivisions(1020);

      h_TrueXS_abs->GetXaxis()->SetRangeUser(300,1000);

      abs_KE->SetTitle( "Geant Absorption XS; Energy [MeV]; #sigma [mb]");
      abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
      abs_KE->SetLineColor(kBlue);
      abs_KE->SetLineWidth(3);
      abs_KE->Draw("AC");
      h_unfoldXS_abs->SetMarkerSize(0.7);
      h_unfoldXS_abs->Draw("PE0 SAME");
      h_TrueXS_abs->SetMarkerColor(8);
      h_TrueXS_abs->SetMarkerStyle(22);
      h_TrueXS_abs->Draw("PE0 SAME");
      h_rawXS_abs->GetXaxis()->SetRangeUser(400,900);
      h_rawXS_abs->SetMarkerColorAlpha(kOrange, 0.3);
      h_rawXS_abs->Draw("P SAME");
      c_unfold_abs->BuildLegend();

      c_unfold_abs->Write();

   }



}

#ifndef __CINT__
int main () { unfold_bayes_wBG( path); return 0; }  // Main program when run stand-alone
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

/*   //======= ADD Muon Content to initE and interE of pions and totInel --> FAKE!! (Extra BG) ===================================================
//Add Muons that were not selected by eventSelection in order to not double-count events
//ist that correct? probably should have many similar muons to those 8% already in the BG...(?)
frame_train
.Filter("primaryMuon && !selected_incidentPion")
.Filter("reco_initKE_rwData != -999. && reco_interKE_rwData != -999.")
.Range(muExtra_incidentPion)
.Foreach( [ &response_initE, &response_interE, &response_totInel] 
( double reco_initE, double reco_interE){

response_initE.Fake( reco_initE );
response_interE.Fake( reco_interE );
response_totInel.Fake( reco_interE );

},
{"reco_initKE_rwData", "reco_interKE_rwData"});

frame_train
.Filter("primaryMuon && !selected_abs")
.Filter("reco_initKE_rwData != -999. && reco_interKE_rwData != -999.")
.Range(muExtra_abs)
.Foreach( [ &response_abs] 
(double reco_interE){
response_abs.Fake( reco_interE );

},
{"reco_interKE_rwData"});
*/

