#include <iostream>
#include <fstream>
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLegend.h>
#include <string>
#include "TCanvas.h"
#include "macro_v1.C"
#include "untuplizer.h"

void xAna_all() {
  
  TFile* file[10];
    
  file[0] = TFile::Open("fatmass_1000.root");
  file[1] = TFile::Open("fatmass_1200.root");
  file[2] = TFile::Open("fatmass_1600.root");
  file[3] = TFile::Open("fatmass_1800.root");
  file[4] = TFile::Open("fatmass_2000.root");
  file[5] = TFile::Open("fatmass_2500.root");
  file[6] = TFile::Open("fatmass_3000.root");
  file[7] = TFile::Open("fatmass_3500.root");
  file[8] = TFile::Open("fatmass_4000.root");
  file[9] = TFile::Open("fatmass_4500.root");
  
  TH1D* h0 = (TH1D*)(file[0]->Get("invmass"));	    
  TH1D* h1 = (TH1D*)(file[1]->Get("invmass"));
  TH1D* h2 = (TH1D*)(file[2]->Get("invmass"));
  TH1D* h3 = (TH1D*)(file[3]->Get("invmass"));
  TH1D* h4 = (TH1D*)(file[4]->Get("invmass"));
  TH1D* h5 = (TH1D*)(file[5]->Get("invmass"));						   
  TH1D* h6 = (TH1D*)(file[6]->Get("invmass"));
  TH1D* h7 = (TH1D*)(file[7]->Get("invmass"));
  TH1D* h8 = (TH1D*)(file[8]->Get("invmass"));
  TH1D* h9 = (TH1D*)(file[9]->Get("invmass"));
  
  TCanvas* c = new TCanvas("c","",0,0,1000,800);
  c->cd();
  h1->Draw();
  h1->SetTitle("Mass reconstruction of 2 Fat jets");
  h1->SetStats(0);
  h1->SetLineColor(kRed-2);
  h0->Draw("same");
  h1->SetLineColor(kPink-2);
  h2->Draw("same");
  h2->SetLineColor(kViolet-2);
  h3->Draw("same");
  h3->SetLineColor(kBlue-2);
  h4->Draw("same");
  h4->SetLineColor(kAzure+8);
  h5->Draw("same");
  h5->SetLineColor(kTeal-7);
  h6->Draw("same");
  h6->SetLineColor(kGreen+3);
  h7->Draw("same");
  h7->SetLineColor(kYellow-7);
  h8->Draw("same");
  h1->SetLineColor(kOrange+1);
  h9->Draw("same");
  h1->SetLineColor(kBlack);

  TLegend *legend = new TLegend(0.6,0.12,0.9,0.3);
  legend->AddEntry(h0,"X mass = 1000 GeV","l");
  legend->AddEntry(h1,"X mass = 1200 GeV", "l");
  legend->AddEntry(h2,"X mass = 1600 GeV", "l");
  legend->AddEntry(h3,"X mass = 1800 GeV", "l");
  legend->AddEntry(h4,"X mass = 2000 GeV", "l");
  legend->AddEntry(h5,"X mass = 2500 GeV", "l");
  legend->AddEntry(h6,"X mass = 3000 GeV", "l");
  legend->AddEntry(h7,"X mass = 3500 GeV", "l");
  legend->AddEntry(h8,"X mass = 4000 GeV", "l");
  legend->AddEntry(h9,"X mass = 4500 GeV", "l");
  legend->SetFillColor(0);
  legend->Draw();
}
		     
