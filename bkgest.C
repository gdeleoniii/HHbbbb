#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <THStack.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TPad.h>
#include <THnSparse.h>
#include <TStyle.h>
#include <TStyle.h>
#include "setNCUStyle.C"
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>

/*
double quadfunc(double*v, double* p) {
  double x= v[0];
  return p[0]+p[1]*x+p[2]*x*x;
}	
*/		    

void backest(std::string var,std::string name){
  TFile* file = TFile::Open("JetHT-Run2015_bkgest.3.2.root");


  TH1F* h1 = (TH1F*)(file->Get(Form("%sone",var.data())));
  TH1F* h2 = (TH1F*)(file->Get(Form("%stwo",var.data())));

  Double_t norm1 = 1/h1->Integral();
  Double_t norm2 = 1/h2->Integral();

  h1->Scale(norm1);
  h2->Scale(norm2);

  setNCUStyle();
  
  TCanvas* c4 = new TCanvas("c4","",0,0,900,900);                                                                                                  
  c4->cd();
  /*
  h1->SetLineColor(kBlack);
  h1->SetLineWidth(2);
  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(0.8);
  h2->SetLineColor(kOrange+4);
  h2->SetFillColor(kOrange+5);
  h2->SetLineWidth(2);
  h1->GetXaxis()->SetTitle(Form("%s",name.data()));
  h2->GetXaxis()->SetTitle(Form("%s",name.data()));
  //h1->Draw("el");
  h2->Draw("hist");
  //h2->Draw("histsame");
  h1->Draw("elsame");
  */
  h1->SetLineWidth(3);
  h1->SetLineColor(kGray+3);
  h2->SetLineWidth(3);
  h2->SetLineColor(kTeal+5);
  h1->GetXaxis()->SetTitle(Form("%s",name.data()));                                                                                               
  h2->GetXaxis()->SetTitle(Form("%s",name.data()));
  h2->Draw("hist");
  h1->Draw("histsame");
  c4->RedrawAxis();
  
  TLegend *leg = new TLegend(0.61, 0.80, 0.79, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  //leg->SetHeader("Pass-fail ratio (leading jet)");
  leg->AddEntry(h1, "one anti-tag", "l");
  leg->AddEntry(h2, "two anti-tag", "l");
  leg->Draw();

  c4->Print(Form("%s_bkgest.2.1.pdf",var.data()));
}
