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


int print_infoTxt(const string mcFilepath){

   gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
   ROOT::RDataFrame inputFrame(pionTree, mcFilepath);
   gStyle->SetNdivisions(1020);
   gStyle->SetPaintTextFormat("3.2f");

   auto frame = inputFrame
      .Filter("true_beam_PDG == 211 && true_beam_endZ > 0");
      //.Range(20,50);

   //ofstream pionFile;
   //pionFile.open("mc_pionInfo.txt");
   std::cout << "PDG true_startKE true_endKE reco_startKE reco_endKE true_endZ reco_calo_endZ end_wire true_abs true_chex true_BG isDecay pandoraReco selected_incidentPion selected_abs" << std::endl;

   frame
      .Foreach([](int pdg, double tr_startKE, double tr_endKE, double rec_startKE, double rec_endKE, double tr_endZ, double rec_endZ, const std::vector<double> &reco_beam_calo_wire, int abs, int chex, bool dec, bool bc, bool selInc, bool selAbs){

            std::cout << pdg  << " " << tr_startKE << " " << tr_endKE << " " << rec_startKE << " " << rec_endKE << " " 
            << tr_endZ << " " << rec_endZ << " " << reco_beam_calo_wire[reco_beam_calo_wire.size() - 1] << " " << abs << " " << chex << " " << dec << " " << bc << " " << selInc << " " << selAbs << " " << std::endl;
            }

            ,{"true_beam_PDG", "true_firstEntryIncident", "true_interactingKE_fromLength", 
            "reco_firstEntryIncident", "reco_interactingKE", "true_beam_endZ", "reco_beam_calo_endZ", "reco_beam_calo_wire", "true_absSignal", "true_chexSignal", "isDecay", "primary_isBeamType", "selected_incidentPion", "selected_abs"});

  
   //pionFile.close();
   return 0;
}

