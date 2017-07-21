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
#include <TMatrixD.h>
#include <TGraphAsymmErrors.h>
#include <TFitResult.h>
#include "setNCUStyle.C"
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>

void ile() {
  TFile* file  = TFile::Open("JetHT-Run2015_passfailratio.4.root");
 
  TProfile* h_data1     = (TProfile*)(file->Get("tau21vsPR_700to1500_before"));
  TProfile* h_data2     = (TProfile*)(file->Get("tau21vsPR_700to1500_after"));
 
  setNCUStyle();
  TCanvas* a = new TCanvas("a","a",0,0,750,550);
  
  h_data1->SetXTitle("h1 jet pruned mass");
  //h_data1->SetYTitle(Form("%s",title.data()));
  h_data1->GetYaxis()->SetTitleOffset(1.4);
  h_data1->SetMarkerStyle(20);
  h_data1->SetMarkerSize(0.8); 
  h_data1->SetLineColor(kGreen+2);
  h_data1->SetMarkerColor(kGreen+2);
  h_data2->SetMarkerStyle(21);
  h_data2->SetMarkerSize(0.8);
  h_data2->SetLineColor(kBlue+2);
  h_data2->SetMarkerColor(kBlue+2);
  h_data1->Draw();
  h_data2->Draw("same");

  TLegend *legendx = new TLegend(0.5,0.5,0.8,0.65);
  //legendx->SetHeader(Form("%s",header.data()));
  legendx->AddEntry(h_data2,"700 GeV < dijet mass < 1500 GeV","lp");
  legendx->AddEntry(h_data1,"dijet mass > 1500 GeV", "lp");
  legendx->SetFillStyle(0);
  legendx->SetTextSize(0.035);
  legendx->Draw();
}

//----------for fitting errors----------------------
TGraphAsymmErrors* fitErrorBand(const TF1* f, TMatrixD* corrMatrix, double (*fitFunc)(double*, double*)){

  double par[3] = {0};

  for( int i = 0; i < 3; i++ )
    par[i] = f->GetParameter(i);

  TF1* posFit[3];
  TF1* negFit[3];

  for( int i = 0; i < 3; i++ ){

    double partemp[3] = {par[0],par[1],par[2]};

    posFit[i] = new TF1(Form("posFit%d",i), fitFunc, -75, 75, 3);
    partemp[i]  = par[i] + f->GetParError(i);
    posFit[i]->SetParameters(partemp[0],partemp[1],partemp[2]);

  }

  for( int i = 0; i < 3; i++ ){

    double partemp[3] = {par[0],par[1],par[2]};

    negFit[i] = new TF1(Form("negFit%d",i), fitFunc, -75, 75, 3);
    partemp[i]  = par[i] - f->GetParError(i);
    negFit[i]->SetParameters(partemp[0],partemp[1],partemp[2]);

  }

  TMatrixD posColM(3,1);
  TMatrixD negColM(3,1);
  TMatrixD posRowM(1,3);
  TMatrixD negRowM(1,3);

  int    NBINS = 75;
  float x     = -75.0;
  double width = 2; //(75-x)/NBINS;

  //double x[] = {-75,-60,-32.5,-0.5,25,57.5,75};
  //const int NBINS=sizeof(x)/sizeof(x[0]);

  double funcX[NBINS];
  double funcY[NBINS];
  double posUnc[NBINS];
  double negUnc[NBINS];

  for( int n = 0; n < NBINS; n++){

    for(int i = 0; i < 3; i++){
    
      posColM(i,0) = fabs(f->Eval(x) - posFit[i]->Eval(x));
      negColM(i,0) = fabs(f->Eval(x) - negFit[i]->Eval(x));
      posRowM(0,i) = posColM(i,0);
      negRowM(0,i) = negColM(i,0);
    
    }

    TMatrixD posTemp = posRowM*(*corrMatrix*posColM);
    TMatrixD negTemp = negRowM*(*corrMatrix*negColM);
    
    posUnc[n] = TMath::Sqrt(posTemp(0,0));
    negUnc[n] = TMath::Sqrt(negTemp(0,0));

    funcX[n] = x;
    funcY[n] = f->Eval(x);

    x += width;

  }


  double cent = 0;
  double upup = 0;
  double dwdw = 0;
  float weight[NBINS]; 
  float var = 0;
  for(int i = 0; i <NBINS; i++) {
    weight[i] = 1/(posUnc[i]*posUnc[i]);
    var += weight[i];
    cent += funcY[i]*weight[i];
    cout<<funcX[i]+125<<" "<<funcY[i]<<" "<<funcY[i]+posUnc[i]<<" "<<funcY[i]-negUnc[i]<<endl;
  }
  upup = (cent/var) + sqrt(1/var);  
  dwdw = (cent/var) - sqrt(1/var);

  //cout<<"average = "<<cent/var<<" +/-  "<<sqrt(1/var)<<" ; "<<upup<<" "<<dwdw<<endl;
  TGraphAsymmErrors* g = new TGraphAsymmErrors(NBINS, funcX, funcY, 0, 0, negUnc, posUnc);

  return g;

}

//------------------------------------------------

double quadfunc(double* v, double* p) {
  //return p[0];
  double x= v[0];
  return p[0]+p[1]*x+p[2]*x*x;
}			    



void pfratio() {
  //TFile* file  = TFile::Open("/Users/gdeleoniii/Documents/GREG/ROOT/ssss/JetHT-Run2015_bkgest.2.1.root");
  TFile* file  = TFile::Open("JetHT-Run2015_passfailratio.7.1.root");
  //TFile* file  = TFile::Open("JetHT-Run2015_passfailratio.9.root");  

  TH1F* h1 = (TH1F*)(file->Get("pass"));
  TH1F* h2 = (TH1F*)(file->Get("fail"));
  TH1F* h4 = (TH1F*)(file->Get("pfratio"));

  
  TH1F* h3 = (TH1F*)h1->Clone("h3");
  h3->Reset();
  h3->Divide(h1,h2);
  
  TF1* f_fitprmass = new TF1("f_fitprmass", quadfunc,-75,75,3);

  //setNCUStyle();
  
  TCanvas* c4 = new TCanvas("c4","",0,0,600,600);                                                                                                  
  c4->cd();
  h3->Fit("f_fitprmass","","",-75,75);
  //h3->GetYaxis()->SetRangeUser(0.06,0.182);
  h3->GetYaxis()->SetRangeUser(0,0.5);
  h3->SetLineColor(kBlack);                                                                                                                  
  h3->GetYaxis()->SetTitle("R_{p/f}");                                                                                                                           
  h3->GetXaxis()->SetTitle("m_{J} - m_{H} (GeV)");                                                                                                      
  h3->GetYaxis()->SetTitleOffset(1.6);                                                                                                                          
  //h3->SetStats(0);                                                                                                                                              
  h3->Draw("e1");

  cout<<"integral = " << (f_fitprmass->Integral(-20,10))/30 << " "<< (f_fitprmass->IntegralError(-20,10))/30 <<endl;

  int nb = h3->GetNbinsX();
  float mean[nb]; //= h3->GetBinContent[nb];
  float error[nb]; // = h3->GetBinError[ 
  float mu = 0;
  float sig = 0;

  for(int i=1;i<=nb;i++) {
    mean[i] = h3->GetBinContent(i);
    error[i] = h3->GetBinError(i);
    if( mean[i] > 0) {
      mu += mean[i]*(1/(error[i]*error[i]));
      sig += 1/(error[i]*error[i]);
      cout<< mu <<" "<<sig<<endl;
    }
  }
  cout<<" weighted ave = "<<mu/sig<< " +/- "<<sqrt(1/sig)<<endl;


  //for(int i=0;i<3;i++) {
  //float relchange = abs(f_fitprmass->GetParameter(i) - f_fitprmass->GetParError(i))/f_fitprmass->GetParameter(i);
  //cout<<"relative change = "<< f_fitprmass->GetParameter(i) << " " << f_fitprmass->GetParError(i) << " " << relchange<< endl;
  //}
  // Calculate the fitting uncertainties by using error matrix

  //TGraphAsymmErrors* g_errorBands = fitErrorBand(f_fitprmass, &corrMatrix, quadfunc);
  
    
  TFitResultPtr fitptr = h3->Fit(f_fitprmass, "QS");
  TFitResult fitresult = (*fitptr);
  TMatrixD corrMatrix  = fitresult.GetCorrelationMatrix();  
  
  TGraphAsymmErrors* g_errorBands = fitErrorBand(f_fitprmass, &corrMatrix, quadfunc);
 
  // *g_errorBands = fitErrorBand(f_fitprmass, &corrMatrix, quadfunc);
  
  g_errorBands->SetLineColor(kWhite);
  g_errorBands->SetFillStyle(1001);
  g_errorBands->SetFillColorAlpha(kYellow-9,0.75);
  g_errorBands->Draw("e3same");
  
  h3->Draw("e1same");
  
  c4->RedrawAxis();
  //h1->SetLineColor(kOrange-4);
  //h1->Draw("e1same");
  
  TLegend *legend = new TLegend(0.2,0.7,0.5,0.85);
  legend->AddEntry(h3,"events used in fit","le");
  legend->AddEntry(f_fitprmass,"fit", "l");
  //legend->AddEntry(g_errorBands,"fit error", "f");
  legend->SetFillStyle(0);
  legend->SetTextSize(0.035);
  legend->Draw();
  
  /*cout<<h3->GetBinContent(1)<<" "<<h4->GetBinContent(1)<<endl;
  cout<<h3->GetBinContent(2)<<" "<<h4->GetBinContent(2)<<endl;
  cout<<h3->GetBinContent(3)<<" "<<h4->GetBinContent(3)<<endl;
  cout<<h3->GetBinContent(4)<<" "<<h4->GetBinContent(4)<<endl;
  cout<<h3->GetBinContent(5)<<" "<<h4->GetBinContent(5)<<endl;
  TCanvas* c2 = new TCanvas("c2","",0,0,600,600);
  c2->cd();
  h4->Fit("f_fitprmass","","",-75,75);
  h4->GetYaxis()->SetTitle("R_{p/f}");
  h4->GetXaxis()->SetTitle("m_{J} - m_{H} (GeV)");
  h4->GetYaxis()->SetTitleOffset(1.6);
  h4->Draw("e1");*/

}
