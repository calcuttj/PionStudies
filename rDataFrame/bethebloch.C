#include "TCanvas.h"
#include "TGraph.h"
#include <math.h>
#include <iostream>
using namespace std;


double bethe(double betagamma, double mc)
{
    double K,rho,Z,A, charge, me, I, beta, gamma, p,c,wmax;
    K = 0.307;
    rho = 1.394;
    charge = 1;
    Z = 18;
    A = 40;
    c = 299792458; 
    I = pow(10,-6)*10.5*18; //MeV
    me = 0.51; //MeV me*c^2
    
    p = betagamma*mc;
    beta = p/sqrt(pow(mc,2) + pow(p,2));
    gamma =  1/sqrt(1 - pow(beta,2));
    
    //wmax = 2*me*pow(c,2)*pow(beta,2)*pow(gamma,2)/(1+2*gamma*me/mc + pow(me,2)/pow(mc,2));
    wmax = 2*me*pow(beta,2)*pow(gamma,2)/(1+2*gamma*me/mc + pow(me,2)/pow(mc,2));
    
    
    double enloss;
    enloss = (K*rho*Z*pow(charge,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2));

   return enloss;
}


void bethebloch() {
    double betagamma[10000], dEdx[10000];
    
    int i;
    double fake;
    double melec, mprot, mmuo, mpi, mka;
    melec = 0.51;
    mprot = 938;
    mmuo = 105.7;
    mpi = 139;
    mka = 493.6;
    
    
  for (i = 0; i < 10000; i = i+1){
      betagamma[i]= 0.1 + (i+1)/100;
}
  
  for (i = 0; i<10000; i = i+1){
  
  	fake = betagamma[i];	
  	dEdx[i] = bethe(fake,mpi);
      
}

 TCanvas *c1 = new TCanvas("c1","A Simple Graph Example"); //,200,10,700,500);
 //c1->SetLogx();
 //c1->SetFillColor(42);
 c1->SetGrid();

TGraph *gr1 = new TGraph (10000, betagamma, dEdx);
gr1->Draw("alp");
  
}
