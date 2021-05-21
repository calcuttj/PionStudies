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

int eSliceMethod_trueProcess_trueE(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame frame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);

   //output file
   string outputNameMC = "output_eSliceMethod_trueProcess_trueE.root";
   string outputName;
   outputName = outputNameMC;


   TFile *output = new TFile ( outputName.c_str() , "RECREATE");

   //access Jakes GeantFile in folder
   TFile f1("exclusive_xsec.root");
   TGraph *totInel_KE = (TGraph*)f1.Get("total_inel_KE");
   TGraph *abs_KE = (TGraph*)f1.Get("abs_KE");
   TGraph *cex_KE = (TGraph*)f1.Get("cex_KE");
   f1.Close();

   string mg_title;
   mg_title = "Cross-Section MC; Reco kinetic Energy (MeV); #sigma (mbarn)";

   TMultiGraph *mg = new TMultiGraph();
   mg->Add(totInel_KE);
   mg->Add(abs_KE);
   //mg->Add(cex_KE);
   mg->SetTitle(mg_title.c_str());

   //switch to output-file
   output->cd();

   //--------------------------------------------------------
   // Initialise incident and interacting Histograms
   // may 21
   // build incident histo now in a different way after having unsmeared start and end distributions
   //
   //--------------------------------------------------------

   //Incident Histogram and Interacting Histogram, interacting for different Processes
   //Energies are kinetic Energy of beam particle



   TH1D* h_trueE_truePion_inc_initE = new TH1D("h_trueE_truePion_inc_initE", "Incident Selected Pion true_initE", nBin_int, eEnd, eStart);
   h_trueE_truePion_inc_initE->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_trueE_truePion_inc_interE = new TH1D("h_trueE_truePion_inc_interE", "Incident Selected Pion true_interE", nBin_int, eEnd, eStart);
   h_trueE_truePion_inc_interE->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_trueE_truePion_incident = new TH1D("h_trueE_truePion_incident", "Incident Selected Pion trueE", nBin_inc, eEnd, eStart);
   h_trueE_truePion_incident->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_trueE_trueAbs_interacting = new TH1D("h_trueE_trueAbs_interacting", "Interacting Selected ABS True Energy", nBin_int, eEnd, eStart);
   h_trueE_trueAbs_interacting->GetXaxis()->SetTitle("True KE (MeV)");

   TH1D* h_index_i= new TH1D("h_index_i", "Distribution Initial KE in LAr", 30, 0, 30);

   //Initial Filters for all events
   auto mcIncident_true_primaryPi = frame
    .Filter("true_beam_endZ > 0")
    //.Filter("primary_isBeamType && passBeamCut")      
    .Filter("true_beam_PDG == 211")
    //.Range(70,100)
    .Define("true_initialKE_LAr", [h_index_i](std::vector<double> &trajKE, std::vector<double> &trajZ){
          double energy;
          int i=0;
          while( trajZ[i] < 0){ //loop through trajZ until it goes above 0
            i++;
          }
         // std::cout << "..................................................................." << std::endl;
         // std::cout << "true_beam_traj_Z[" << i << "] = " << trajZ[i] << "; true_beam_traj_KE[" <<i<< "] = " << trajKE[i] << std::endl; 
         // std::cout << "Energy at that point            = " << trajKE[i] << std::endl;
         // //std::cout << "First Entry Incident            = " << trajKE[0] << std::endl;
         // std::cout << "Z+1 coordinate                  = " << trajZ[i + 1] << std::endl;
         // std::cout << "Energy at that point            = " << trajKE[i + 1] << std::endl;
         // h_index_i->Fill(trajKE[i + 1]);
         h_index_i->Fill(i);
         //return energy = 820;
         return energy = trajKE[0] - 50;
          }
          , {"true_beam_traj_KE", "true_beam_traj_Z_SCE"});
      
   //.Filter("true_primPionInel");

   auto mcInteracting_true_abs = mcIncident_true_primaryPi
      .Filter("true_absSignal && true_pion_daughter == 0");
   //========================================================
   //Build the Incident Histogram
   //---------
   //take the interacting energy, that slice is filled too
   mcIncident_true_primaryPi
      .Foreach( [h_trueE_truePion_inc_initE, h_trueE_truePion_inc_interE, h_trueE_truePion_incident] (double true_firstEntryIncident, double true_beam_interactingEnergy) { 
            //make sure incident Pion does not interact in bin it was born
            //if( !( initE_sameBin_interE( true_firstEntryIncident, true_beam_interactingEnergy ) ) ){
            //if(true ){
            int binNumber_initEnergy = (int) true_firstEntryIncident / bin_size_inc + 1;
            int binNumber_interEnergy = (int) true_beam_interactingEnergy / bin_size_inc + 1;


         if(binNumber_initEnergy != binNumber_interEnergy && true_firstEntryIncident > true_beam_interactingEnergy){
         if(binNumber_initEnergy <= nBin_int && binNumber_initEnergy >0 && binNumber_interEnergy <= nBin_int && binNumber_interEnergy > 0){

            //std::cout << "binNumber initial E = " << binNumber_initEnergy << std::endl;
            //std::cout << "binNumber inter   E = " << binNumber_interEnergy << std::endl;

            h_trueE_truePion_inc_initE->SetBinContent( binNumber_initEnergy, h_trueE_truePion_inc_initE->GetBinContent( binNumber_initEnergy ) + 1); 
            h_trueE_truePion_inc_interE->SetBinContent( binNumber_interEnergy, h_trueE_truePion_inc_interE->GetBinContent( binNumber_interEnergy ) + 1); 

             
            /*for(int i = binNumber_interEnergy; i <= binNumber_initEnergy; i++){      
            if( i > 0 && i <= nBin_inc){ //make sure we don't go outside of bin range
            h_trueE_truePion_incident->SetBinContent( i, h_trueE_truePion_incident->GetBinContent(i) + 1 ); 
            //increment previous bin content by one. should handle number of entries in hist. unlike AddBinContent
            };
            };*/
            };   
            };
            }
            ,{"true_initialKE_LAr", "true_KEint_fromEndP"});

      //h_trueE_truePion_incident->Sumw2(0);
      //h_trueE_truePion_incident->Rebin( bin_size_int/bin_size_inc );
      //h_trueE_truePion_incident->Scale( 1 / (bin_size_int/bin_size_inc) );


            h_trueE_truePion_inc_initE->Write();
            h_trueE_truePion_inc_interE->Write();


            //=====================================================
            //------------------------------------------------------
            //Interacting selected samples
            //------------------------------------------------------
            mcInteracting_true_abs
               .Foreach( [h_trueE_trueAbs_interacting] ( double true_firstEntryIncident, double true_beam_interactingEnergy ){

               int binNumber_initEnergy = (int) true_firstEntryIncident / bin_size_inc + 1;
               int binNumber_interEnergy = (int) true_beam_interactingEnergy / bin_size_inc + 1;

               if(binNumber_initEnergy != binNumber_interEnergy && true_firstEntryIncident > true_beam_interactingEnergy){
                     //if( !( initE_sameBin_interE( true_firstEntryIncident, true_beam_interactingEnergy ) ) ){
                     h_trueE_trueAbs_interacting->SetBinContent( binNumber_interEnergy, h_trueE_trueAbs_interacting->GetBinContent( binNumber_interEnergy ) + 1); 
                     //}
                     }
                     }            
                     ,{"true_initialKE_LAr","true_KEint_fromEndP"});

            h_trueE_trueAbs_interacting->Sumw2(0);
            h_trueE_trueAbs_interacting->Write();
            
            //=====================================================
            //           Debug
            //=====================================================

            
            std::cout << "Nbins Histogram = " << nBin_int << "  Energy from 0, 1200 MeV " << std::endl;
            std::cout << "Bin             Birth     Death     Cumulative Shifted      Ratio " << std::endl;

            
            int birth = 0;
            int death = 0;
            int cumulative_shifted = 0;
            double ratio = 0;
            int cum_before = 0;
            for(int i = nBin_int; i >= 1; i--){

               cumulative_shifted = birth + cum_before - death;
               cum_before = cumulative_shifted;
               birth = h_trueE_truePion_inc_initE->GetBinContent(i);
               death = h_trueE_truePion_inc_interE->GetBinContent(i);
               ratio = (double) death / cumulative_shifted;

               std::cout << i << "  " << birth << "     " << death << "   " << cumulative_shifted << "    " << ratio << std::endl;

              }

            TH1D* cum_initE = new TH1D("cum_initE", "cum_initE", nBin_int, eEnd, eStart);
            TH1D* cum_interE = new TH1D("cum_interE", "cum_interE", nBin_int, eEnd, eStart);
            for(int i = nBin_int; i >= 1; i--){
               cum_initE->SetBinContent(i,0);
               cum_interE->SetBinContent(i,0);
            };

            for(int i = nBin_int; i >= 1; i --){

               int nascita = h_trueE_truePion_inc_initE->GetBinContent( i );
               int morte = h_trueE_truePion_inc_interE->GetBinContent( i );
               
               for(int j = i - 1; j >= 1; j --){

                  h_trueE_truePion_incident->SetBinContent( j , h_trueE_truePion_incident->GetBinContent(j) + nascita);
               };

               for(int j = i - 1; j >= 1; j --){

                  h_trueE_truePion_incident->SetBinContent( j , h_trueE_truePion_incident->GetBinContent(j) - morte);
               };


            };

            h_trueE_truePion_incident->Write();

            for(int i = nBin_int; i >=1; i--) std::cout << "Entry[" << i << "] = " << h_trueE_truePion_incident->GetBinContent(i)  << std::endl;

            //Incident[i] = Birth[i-1] + Incident[i-1] - Death[i-1]
            /*h_trueE_truePion_incident->SetBinContent( nBin_int , 0);

            for(int i = nBin_int - 1; i >= 1; i--){

               int entry_incident = h_trueE_truePion_inc_initE->GetBinContent( i+1 )
                                  + h_trueE_truePion_incident->GetBinContent( i+1 )
                                  - h_trueE_truePion_inc_interE->GetBinContent( i+1 );

               h_trueE_truePion_incident->SetBinContent( i, entry_incident );

               


            }*/


            

            //=====================================================
            //            Prepare BetheBloch Mean for each Bin 
            //            QUESTION: take betheBloch of Pion or Muon?? Comparison to data fits better muon Bethe... 
            //            at hihger momentum ~400-800 they anway are almost the same
            //=====================================================
            TH1D* h_betheMean_muon = new TH1D("h_betheMean_muon", "Mean Energy Loss", nBin_int, eEnd, eStart);

            //fill histo with Mean dEdX of bin center
            for(int i = 1; i <= nBin_int; i++){
               h_betheMean_muon->SetBinContent(i , betheBloch( eEnd + (i - 0.5)*bin_size_int  , mass_pion) );
            };

            h_betheMean_muon->Write();


            //TH1* h_cum_initE = h_trueE_truePion_inc_initE->GetCumulative(false);
            //TH1* h_cum_interE = h_trueE_truePion_inc_interE->GetCumulative(false);
            //FILLING OF INCIDENT
            //--> if Pion is "born" in bin 10 and interacts in bin 4 then incident is filled from bin 9 to 4
            //cumulative sum of initialE: the difference between higher bin and the lower next shows how many pions were born 
            //cumulative sum of interE: the difference between the higher bin and lower shows how many pions died
            //incident histo is built such that it takes the difference of what was born in previous bin compared to what died in previous bin
            //i.e. incident histo cannot fill the first energy bin of KE 1200MeV

            /*for(int i= nBin_int-1; i >= 1; i--){
              double diff = h_cum_initE->GetBinContent(i+1) - h_cum_interE->GetBinContent(i+1);
              h_trueE_truePion_incident->SetBinContent( i, diff);
              };

              h_trueE_truePion_incident->Sumw2(0);
              h_trueE_truePion_incident->Write();
            */

              /*for(int i= nBin_int; i > 1; i--){
              double diff = h_cum_initE->GetBinContent(i) - h_cum_interE->GetBinContent(i);
              h_trueE_truePion_incident->SetBinContent( i-1, diff);
              };

              h_trueE_truePion_incident->Sumw2(0);
              h_trueE_truePion_incident->Write();*/
              
            //=====================================================
            //             Computing the XS
            //=====================================================
            //
            // xs(Ebin) = (A / (Na*density*bin_size)) * dEdX(Ebin) * hInteracting / hIncident
            //
            // More Accurate use log( Ninc / (Ninc - Nint ))
            //
            double scale_factor = factor_mbarn * atomic_mass / ( density * N_avogadro * bin_size_int );

            //------------------------------------------------------
            //    Absorption, Selected Interactions Reconstrucetd Energy
            //------------------------------------------------------
            //
            string xs_abs_name, xs_abs_title;
            xs_abs_name = "h_xs_trueE_trueAbs_mc";
            xs_abs_title = "Absorption MC";


            TH1D* h_xs_trueE_trueAbs = (TH1D*) h_trueE_truePion_incident->Clone(xs_abs_name.c_str());
            h_xs_trueE_trueAbs->SetTitle(xs_abs_title.c_str());

            TH1D* dummy_recoE_abs = (TH1D*) h_trueE_truePion_incident->Clone(xs_abs_name.c_str());
            dummy_recoE_abs->Add( h_trueE_trueAbs_interacting, -1);
            h_xs_trueE_trueAbs->Divide( dummy_recoE_abs );

            for(int i = 1; i <= nBin_int; i++) {
               if(h_xs_trueE_trueAbs->GetBinContent(i) > 0) h_xs_trueE_trueAbs->SetBinContent(i, log( h_xs_trueE_trueAbs->GetBinContent(i) ));
               else h_xs_trueE_trueAbs->SetBinContent(i, -10);

            };
            h_xs_trueE_trueAbs->Multiply( h_betheMean_muon );
            h_xs_trueE_trueAbs->Scale( scale_factor );

            for(int i=1; i <= nBin_int; i++){

               double p = h_trueE_trueAbs_interacting->GetBinContent(i) / h_trueE_truePion_incident->GetBinContent(i);
               double nInc_i = h_trueE_truePion_incident->GetBinContent(i);
               double help_factor = scale_factor*h_betheMean_muon->GetBinContent(i);

               h_xs_trueE_trueAbs->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
            };
            h_trueE_trueAbs_interacting->Write();

            h_xs_trueE_trueAbs->Write();


            //------------------------------------------------------
            //    Total Inelastic, Selected Interactions Reconstrucetd Energy
            //------------------------------------------------------
            //
            
            TH1D* h_xs_trueE_truePion_totInel = (TH1D*) h_trueE_truePion_incident->Clone("h_xs_trueE_truePion_totInel_mc");
            h_xs_trueE_truePion_totInel->SetTitle("Total Inelastic MC");

            TH1D* dummy_trueE_totInel = (TH1D*) h_trueE_truePion_incident->Clone("");
            dummy_trueE_totInel->Add( h_trueE_truePion_inc_interE, -1);
            h_xs_trueE_truePion_totInel->Divide( dummy_trueE_totInel );
            
            for(int i = 1; i <= nBin_int; i++) {
               if(h_xs_trueE_truePion_totInel->GetBinContent(i) > 0) h_xs_trueE_truePion_totInel->SetBinContent(i, log( h_xs_trueE_truePion_totInel->GetBinContent(i) ));
               else h_xs_trueE_truePion_totInel->SetBinContent(i, -10);
            };
         
            /*TH1D* h_xs_trueE_truePion_totInel = (TH1D*) h_trueE_truePion_inc_interE->Clone("h_xs_trueE_truePion_totInel_mc");
            h_xs_trueE_truePion_totInel->Divide( h_trueE_truePion_incident ); */
            

            h_xs_trueE_truePion_totInel->Multiply( h_betheMean_muon );
            h_xs_trueE_truePion_totInel->Scale( scale_factor );

            for(int i=1; i <= nBin_int; i++){

               double p = h_trueE_truePion_inc_interE->GetBinContent(i) / h_trueE_truePion_incident->GetBinContent(i);
               double nInc_i = h_trueE_truePion_incident->GetBinContent(i);
               double help_factor = scale_factor*h_betheMean_muon->GetBinContent(i);

               h_xs_trueE_truePion_totInel->SetBinError( i , help_factor*sqrt( p*(1-p)/nInc_i ));
            };
            h_trueE_truePion_inc_interE->Write();
            //h_trueE_truePion_incident->Write();
            h_xs_trueE_truePion_totInel->Write();
            //------------------------------------------------------

            //=====================================================*
            //            Plotting and Style
            //=====================================================
            //
            TCanvas *c_trueE_abs = new TCanvas("c_trueE_abs", "c_trueE_abs");
            gPad->SetGrid(1,1);
            h_xs_trueE_trueAbs->SetTitle( "True Absorption;true Kinetic Energy (MeV); #sigma (mb)");
            h_xs_trueE_trueAbs->GetXaxis()->SetRangeUser(400,1000);
            h_xs_trueE_trueAbs->GetXaxis()->SetNdivisions(1020);
            h_xs_trueE_trueAbs->GetYaxis()->SetNdivisions(1020);

            abs_KE->SetTitle( "Absorption;true Kinetic Energy (MeV); #sigma (mb)");
            abs_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
            abs_KE->SetLineColor(kRed);
            abs_KE->SetLineWidth(3);
            abs_KE->Draw("AC");
            h_xs_trueE_trueAbs->SetMarkerSize(0.7);
            h_xs_trueE_trueAbs->Draw("PE0 SAME");

            c_trueE_abs->Write();


            TCanvas *c_trueE_totInel = new TCanvas("c_trueE_totInel", "c_trueE_totInel");
            gPad->SetGrid(1,1);
            h_xs_trueE_truePion_totInel->SetTitle( "True Total Inelastic; true Kinetic Energy (MeV); #sigma (mb)");
            h_xs_trueE_truePion_totInel->GetXaxis()->SetRangeUser(400,1000);
            h_xs_trueE_truePion_totInel->GetXaxis()->SetNdivisions(1020);
            h_xs_trueE_truePion_totInel->GetYaxis()->SetNdivisions(1020);

            totInel_KE->SetTitle("Total Inelastic MC; true kinetic Energy (MeV); #sigma (mbarn)");
            totInel_KE->GetXaxis()->SetRangeUser(eEnd, eStart);
            totInel_KE->SetLineColor(kRed);
            totInel_KE->SetLineWidth(3);
            totInel_KE->Draw("AC");
            h_xs_trueE_truePion_totInel->SetMarkerSize(0.7);
            h_xs_trueE_truePion_totInel->Draw("PE0 SAME");

            c_trueE_totInel->Write();
            //output->Write();
            //f1.Close();
   TCanvas *c_all = new TCanvas("c_all", "c_all");
   gPad->SetGrid(1,1);
   
   //cex_KE->SetLineColor(kGreen);
   abs_KE->SetLineColor(kBlue);
   mg->GetXaxis()->SetRangeUser(0,1100);
   mg->Draw("AC");
   h_xs_trueE_truePion_totInel->Draw("PE0 SAME");
   h_xs_trueE_trueAbs->Draw("PE0 SAME");

   c_all->Write();
            return 0;
            }

