void gaus_fit(int minEntries, TH2D* h2, TH1D *fit_mean, TH1D *fit_std, TH1D *fit_chi2, TH1D *fit_ndf){

   double temp_maxEntry = 0., temp_std = 0.;

   for(int i = 1; i <= h2->GetNbinsX(); i++){
      //for(int i = 383; i <= 390; i++){

      TH1D *temp_projection = h2->ProjectionY("py",i,i);
  
      //fit only if Entries > 100
      if( temp_projection->GetEntries() > minEntries){
         
         //bin with highest entry WATCHOUT hard coded Bins..
         int max_bin = temp_projection->GetMaximumBin();
      
         temp_maxEntry = temp_projection->GetBinCenter( max_bin );
         temp_std = temp_projection->GetStdDev(1);
         TF1* f1 = new TF1("f1", "gaus", temp_maxEntry - temp_std , temp_maxEntry + temp_std);        
         temp_projection->Fit("f1", "RS QN");
         //gPad->WaitPrimitive();
         //save Parameters of Gaus Fit each to a histo
         fit_mean->SetBinContent( i , f1->GetParameter(1) );
         fit_mean->SetBinError( i , f1->GetParError(1) );
         fit_std->SetBinContent( i , f1->GetParameter(2) );
         fit_std->SetBinError( i , f1->GetParError(2) );
         fit_chi2->SetBinContent( i , f1->GetChisquare() );
         fit_ndf->SetBinContent( i , f1->GetNDF() );

         //Delete fitted function and histo
         delete f1;
         delete temp_projection;
      }
   };
};

void landau_fit(int minEntries, TH2D* h2, TH1D *fit_mpv, TH1D *fit_std, TH1D *fit_chi2, TH1D *fit_ndf){

   double temp_maxEntry = 0.;
    for(int i = 1; i <= h2->GetNbinsX(); i++){
   //   for(int i = 100; i <= 120; i++){

      TH1D *temp_projection = h2->ProjectionY("py",i,i);
      //temp_projection->Draw();
      //gPad->WaitPrimitive();
      //temp_projection->Close();

      //fit only if Entries > 100
      if( temp_projection->GetEntries() > minEntries){

         //bin with highest entry WATCHOUT hard coded Bins..
         int max_bin = temp_projection->GetMaximumBin();
         int widthLeft = max_bin - 100;
         int widthRight = max_bin + 100;

         temp_maxEntry = temp_projection->GetBinContent( max_bin );
         widthLeft = temp_projection->FindFirstBinAbove( 0.2*temp_maxEntry, 1, widthLeft, max_bin );
         widthRight = temp_projection->FindLastBinAbove( 0.2*temp_maxEntry, 1, max_bin , widthRight );

         TF1* f2 = new TF1("f2", "landau", temp_projection->GetBinCenter(widthLeft) , temp_projection->GetBinCenter(widthRight));
         //temp_projection->Fit("f2", "RS ");
         temp_projection->Fit("f2", "RS QN");
         //gPad->WaitPrimitive();
         //save Parameters of Landau Fit each to a histo
         //Landa Par1 MPV is not real MPV, do GetMaximumX for landau most probable value
         fit_mpv->SetBinContent( i , f2->GetMaximumX() );
         fit_mpv->SetBinError( i , f2->GetParError(1) );
         fit_std->SetBinContent( i , f2->GetParameter(2) );
         fit_std->SetBinError( i , f2->GetParError(2) );
         fit_chi2->SetBinContent( i , f2->GetChisquare() );
         fit_ndf->SetBinContent( i , f2->GetNDF() );

         //Delete fitted function and histo
         delete f2;
         delete temp_projection;
      }
   };
};


