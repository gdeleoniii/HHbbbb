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
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include "setNCUStyle.C"
#include <string>

using namespace std;

void plotMultiGraphs(){

  std::string prefix[]={"leading","subleading","both"};
  std::string name[]={"SD","SDCorr","SDCorrThea","PR","PRCorr"};

  const int NTYPES=5;
  float mass[9]={1200,1400,1600,1800,2000,2500,3000,4000,4500};
  int MARKERS[6]={20,21,22,23,34,29};
  int COLORS[NTYPES]={1,4,2,kOrange,kGreen+2};
  setNCUStyle(true);

  TCanvas* c1 = new TCanvas("c1","",700,500);

  for(int i=0; i<3; i++){

    TGraphErrors* graph_mean[NTYPES];
    TGraphErrors* graph_RMS[NTYPES];
    TGraphErrors* graph_RMSMean[NTYPES];

    TMultiGraph *mg = new TMultiGraph();
    TMultiGraph *mg_h = new TMultiGraph();
    TMultiGraph *mg_a = new TMultiGraph();

    for(int j=0; j<NTYPES; j++)
      {
	ifstream fin;
	fin.open(Form("%s_%s.dat",prefix[i].data(),name[j].data()));
	
	float mean[9],meanerr[9],RMS[9],RMSerr[9];
	float rel[9], relerr[9];
	float masserr[9];
	for(int il=0; il<9; il++)
	  {
	    fin >> mean[il] >> meanerr[il] >> RMS[il] >> RMSerr[il];
	    masserr[il]=100;
	  }
	fin.close();

	
	graph_mean[j] = new TGraphErrors(9,mass,mean,masserr,meanerr);
	graph_mean[j]->SetName(Form("gr_Mean_%d",j));
	graph_mean[j]->SetMarkerStyle(MARKERS[j]);
	graph_mean[j]->SetMarkerColor(COLORS[j]);
	graph_mean[j]->SetLineColor(COLORS[j]);
	graph_mean[j]->SetMarkerSize(1.1);
	mg->Add(graph_mean[j]);

	
	
	graph_RMS[j] = new TGraphErrors(9,mass,RMS,masserr,RMSerr);
	graph_RMS[j]->SetName(Form("gr_RMS_%d",j));
	graph_RMS[j]->SetMarkerStyle(MARKERS[j]);
	graph_RMS[j]->SetMarkerColor(COLORS[j]);
	graph_RMS[j]->SetLineColor(COLORS[j]);
	graph_RMS[j]->SetMarkerSize(1.2);
	mg_h->Add(graph_RMS[j]);

	ifstream fin2;
	fin2.open(Form("rel_%s_%s.dat",prefix[i].data(),name[j].data()));
	
	for(int il=0; il<9; il++)
	  {
	    fin2 >> mean[il] >> meanerr[il] >> RMS[il] >> RMSerr[il];
	    float yield = RMS[il]/mean[il];
	    float err = yield*sqrt(pow(RMSerr[il]/RMS[il],2)+
				   pow(meanerr[il]/mean[il],2));
	    rel[il] = yield;
	    relerr[il] = err;
	  }
	fin2.close();

	
	graph_RMSMean[j] = new TGraphErrors(9,mass,rel,masserr,relerr);
	graph_RMSMean[j]->SetName(Form("gr_RMSMean_%d",j));
	graph_RMSMean[j]->SetMarkerStyle(MARKERS[j]);
	graph_RMSMean[j]->SetMarkerColor(COLORS[j]);
	graph_RMSMean[j]->SetLineColor(COLORS[j]);
	graph_RMSMean[j]->SetMarkerSize(1.1);
	mg_a->Add(graph_RMSMean[j]);

	
      } // end loop of mass types

    mg->SetTitle(i==2? Form("%s jets",prefix[i].data()):
		 Form("%s jet",prefix[i].data())
		 );


    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("M_{bulkG} [GeV]");
    mg->GetYaxis()->SetTitleOffset(1.1);
    mg->GetYaxis()->SetTitle("Mean of (Mass-125)/125");
    mg->GetYaxis()->SetRangeUser(-0.3,0.2);
    
    TLegend* leg = new TLegend(0.148,0.634,0.397,0.877);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.05);
    leg->SetBorderSize(0);
    for(int itype=0; itype<NTYPES; itype++)
      {
	leg->AddEntry(Form("gr_Mean_%d",itype), Form("%s",name[itype].data()),"p");  
      }
    
    leg->Draw("same");

    string output = Form("MassMean_%s",prefix[i].data());
    string final;
    final = output + ".gif";
    c1->Print(final.data());
    final = output + ".pdf";
    c1->Print(final.data());

    mg_h->SetTitle(i==2? Form("%s jets",prefix[i].data()):
		 Form("%s jet",prefix[i].data())
		 );
    
    mg_h->Draw("AP");
    mg_h->GetXaxis()->SetTitle("M_{bulkG} [GeV]");
    mg_h->GetYaxis()->SetTitleOffset(1.1);
    mg_h->GetYaxis()->SetTitle("RMS of (Mass-125)/125");
    mg_h->GetYaxis()->SetRangeUser(0.08,0.25);
    
    leg->Draw("same");
    
    output = Form("MassRMS_%s",prefix[i].data());
    final = output + ".gif";
    c1->Print(final.data());
    final = output + ".pdf";
    c1->Print(final.data());



    mg_a->SetTitle(i==2? Form("%s jets",prefix[i].data()):
		 Form("%s jet",prefix[i].data())
		 );
    
    mg_a->Draw("AP");
    mg_a->GetXaxis()->SetTitle("M_{bulkG} [GeV]");
    mg_a->GetYaxis()->SetTitleOffset(1.1);
    mg_a->GetYaxis()->SetTitle("RMS/Mean of mass");
    mg_a->GetYaxis()->SetRangeUser(0.08,0.25);
    
    leg->Draw("same");
    
    output = Form("MassRMSMean_%s",prefix[i].data());
    final = output + ".gif";
    c1->Print(final.data());
    final = output + ".pdf";
    c1->Print(final.data());

    
  } // end loop of jet types



}
