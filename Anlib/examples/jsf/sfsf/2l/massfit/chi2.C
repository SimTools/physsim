// func->SetParameters(143.1693,119.2333,199.539);
func->SetParameters(143.4,119.412,199.539);

func->FixParameter(0,143.1693);
func->FixParameter(1,119.2333);

Double_t sum=0;                     
{
   for (Int_t i=0; i<hElFinal->GetNbinsX(); i++) { 
      Double_t yy = hElFinal->GetBinContent(i);      
      Double_t ey = hElFinal->GetBinError(i);        
      if (yy > 0) {                                   
         Double_t xx   = hElFinal->GetBinCenter(i);       
         Double_t f    = func->Eval(xx);                   
         Double_t dsum = TMath::Power((yy-f)/ey,2);              
         sum += dsum;
         cerr << " x " << xx  << " y " << yy << " dy " << ey << " dsum " << dsum << " sum " << sum << endl;
      }
   }
}

