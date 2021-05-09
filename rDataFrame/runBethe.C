#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "betheBloch.h"
#include <math.h>
#include <iostream>
using namespace std;


void runBethe(double mass_particle) {
   int n_points = 20000;
   double energy[n_points], dEdX_mean[n_points], dEdX_mpv[n_points], betGam[n_points], ratio[n_points];
   

   for (int i = 0; i < n_points; i = i+1){
      energy[i]= 0.1 + (i+1)/3;

      dEdX_mean[i] = betheBloch( energy[i], mass_particle );
      dEdX_mpv[i] = betheBloch_mpv( energy[i], mass_particle );
      ratio[i] = dEdX_mpv[i] / dEdX_mean[i];
      betGam[i] = betaGamma( energy[i] );

   }


   TCanvas *c1 = new TCanvas("c1","Mean Energy Loss and MPV Bethe Bloch"); //,200,10,700,500);
   //c1->SetLogx();
   //c1->SetFillColor(42);
   c1->Divide(2,1);
   c1->cd(1)->SetGrid();
   c1->cd(2)->SetGrid();

   TMultiGraph *mg = new TMultiGraph();

   TGraph *gr1 = new TGraph (n_points, energy, dEdX_mean);
   gr1->SetTitle("Mean Energy Loss");
   gr1->SetLineColor(kBlue);
   gr1->SetMarkerColor(kBlue);
   TGraph *gr2 = new TGraph (n_points, energy, dEdX_mpv);
   gr2->SetTitle("MPV Energy Loss");
   gr2->SetLineColor(kRed);
   gr2->SetMarkerColor(kRed);
   
   TGraph *gr3 = new TGraph (n_points, energy, ratio);
   gr3->SetTitle("Ratio MPV/Mean Energy Loss");
   
   mg->Add(gr1);
   mg->Add(gr2);
   c1->cd(1);
   mg->Draw("alp");
   mg->GetXaxis()->SetTitle("kinetic Energy (MeV)");
   c1->cd(1)->BuildLegend();

   c1->cd(2);
   gr3->Draw("alp");
   gr3->GetXaxis()->SetTitle("kinetic Energy (MeV)");

}
