// example code to run Bulk Graviton->ZZ->ZlepZhad selections on electron-channel

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1D.h>
#include <TFile.h>
#include <TF1.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TObject.h>
#include "setNCUStyle.C"

using namespace std;

float FWHM(TH1F* hist)
{
  int bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2);
  int bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2);
  float fwhm = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);
  // return (fwhm/2.36);
  return (fwhm     );
  //  return fwhm;
}


void plotAllMassVariables(std::string inputFile){

  setNCUStyle();
  TString outputFile;
  outputFile=gSystem->GetFromPipe(Form("file=%s; test=${file%%.root}; echo \"${test}.pdf\"",inputFile.data()));
  cout << "output file name = " << outputFile.Data() << endl;

  TString HeaderName;
  HeaderName=gSystem->GetFromPipe(Form("file=%s; test=${file%%.root}; echo \"${test}\"",inputFile.data()));

  const int NTYPES=5;
  const int NHISTOS=3;
  TFile *inf = new TFile(inputFile.data());
  TH1F* hmass[NTYPES][NHISTOS];
  TH1F* hdiffmass[NTYPES][NHISTOS];
  int COLORS[] ={1,2,4,8,6,7,8};


  std::string prefix[]={"leading","subleading","both"};
  std::string name[]={"SD","SDCorr","SDCorrThea","PR","PRCorr"};

  float max[3]={-9999,-9999,-9999};
  float maxdiff[3]={-9999,-9999,-9999};
  
  for(int i=0; i < NTYPES;i++){
    for(int k=0; k < NHISTOS; k++){
      hmass[i][k] = (TH1F*)inf->FindObjectAny(Form("h_%s_%s",name[i].data(),prefix[k].data()));
      hmass[i][k]->Scale(1.0/hmass[i][k]->Integral());
      if( hmass[i][k]->GetMaximum()>max[k])
	max[k]=hmass[i][k]->GetMaximum();

      hdiffmass[i][k] = (TH1F*)inf->FindObjectAny(Form("h_diff_%s_%s",name[i].data(),prefix[k].data()));						  
      hdiffmass[i][k]->Scale(1.0/hdiffmass[i][k]->Integral());
      if( hdiffmass[i][k]->GetMaximum()>maxdiff[k])
	maxdiff[k]=hdiffmass[i][k]->GetMaximum();
    }
  }
  
  TCanvas* c1 = new TCanvas("c1","",500,500);

  TLegend* leg = new TLegend(0.179,0.592,0.426,0.882);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);



  int nPage=0;
  for(int k=0; k < NHISTOS; k++){
    for(int i=0; i < NTYPES;i++){

      hmass[i][k]->SetMaximum(max[k]*1.1);
      hmass[i][k]->SetLineWidth(3);
      hmass[i][k]->SetLineColor(COLORS[i]);
      hmass[i][k]->SetXTitle("Mass [GeV]");
      hmass[i][k]->SetTitle(k==2? Form("%s, %s jets",HeaderName.Data(),prefix[k].data()):
			    Form("%s, %s jet",HeaderName.Data(),prefix[k].data())
			    );
      if(i==0)
	hmass[i][k] ->Draw();
      else
	hmass[i][k] ->Draw("same");
      if(k==0)
	leg->AddEntry(hmass[i][k], name[i].data(),"l");
      if(i==NTYPES-1)
	leg->Draw("same");
      if(i==NTYPES-1 && nPage==0)
	{
	  c1->Print(Form("%s(",outputFile.Data()),"pdf");
	  nPage++;
	}
      else if(i==NTYPES-1)
	{
	  c1->Print(Form("%s",outputFile.Data()),"pdf");
	  nPage++;
	}
    }

    TLegend* leg2= new TLegend(0.657,0.263,0.904,0.889);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.035);
    leg2->SetBorderSize(0);

    for(int i=0; i < NTYPES;i++){

      hdiffmass[i][k]->SetMaximum(maxdiff[k]*1.5);
      hdiffmass[i][k]->SetLineWidth(3);
      hdiffmass[i][k]->SetLineColor(COLORS[i]);
      hdiffmass[i][k]->SetXTitle("(Mass-125)/125 [GeV]");
      hdiffmass[i][k]->SetTitle(k==2? Form("%s, %s jets",HeaderName.Data(),prefix[k].data()):
			    Form("%s, %s jet",HeaderName.Data(),prefix[k].data())
			    );
      if(i==0)
     	hdiffmass[i][k] ->Draw();
      else
     	hdiffmass[i][k] ->Draw("same");
      ofstream fout;
      fout.open(Form("%s_%s.dat",prefix[k].data(),name[i].data()),ios::out | ios::app);

      string tagname="RMS = ";
      tagname += Form("%.3f",hdiffmass[i][k]->GetRMS());
      string tagname2="FWHM = ";
      tagname2 += Form("%.3f",FWHM(hdiffmass[i][k])); 
      string tagname3="Mean = ";
      tagname3 += Form("%.3f",hdiffmass[i][k]->GetMean());

      fout << hdiffmass[i][k]->GetMean() << " " << hdiffmass[i][k]->GetMeanError() << " " << hdiffmass[i][k]->GetRMS() << " " << hdiffmass[i][k]->GetRMSError()  << endl;
      fout.close();

      ofstream fout2;
      fout2.open(Form("rel_%s_%s.dat",prefix[k].data(),name[i].data()),ios::out | ios::app);

      fout2 << hmass[i][k]->GetMean() << " " << hmass[i][k]->GetMeanError() << " " << hmass[i][k]->GetRMS() << " " << hmass[i][k]->GetRMSError()  << endl;
      fout2.close();
 
      leg2->AddEntry(hdiffmass[i][k], name[i].data(),"l");
      leg2->AddEntry((TObject*)0, tagname3.data(),"");
      leg2->AddEntry((TObject*)0, tagname.data(),"");
      // leg2->AddEntry((TObject*)0, tagname2.data(),"");
      leg2->Draw("same");
      if(i==NTYPES-1)
	{
	  leg2->Draw("same");
	  if(k==2)
	    c1->Print(Form("%s)",outputFile.Data()),"pdf");
	  else
	    c1->Print(Form("%s",outputFile.Data()),"pdf");
	  nPage++;
	}
     }
  } // end loop over  jet type
 
 


}
