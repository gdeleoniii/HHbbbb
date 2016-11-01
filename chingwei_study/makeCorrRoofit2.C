#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include"setNCUStyle.C"
#include "TF1.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooTFnBinding.h" 
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "TVirtualFFT.h"
#include "RooBifurGauss.h"
#include "RooHistPdf.h"
#include "RooBreitWigner.h"
#include "RooNovosibirsk.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TText.h"
#include "RooChi2Var.h"
#include "RooCBShape.h"
#include "RooVoigtian.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooBifurGauss.h"

TCanvas* c1;

//using namespace RooFit ;

void makeCorrRoofit2(){
	setNCUStyle(true);
	c1 = new TCanvas("c1","",1360,768);
	
	TFile *f;
	f=TFile::Open("corr2/corr.root");
	
	TFile* tf1[10];
	int masspoint[10]={600,800,1000,1200,1400,1600,1800,2000,2500,3000};
	for(int i=0;i<10;i++){
		tf1[i]=TFile::Open(Form("corr2/BulkGravToWW_%d.root",masspoint[i]));
		//f[i+11]=TFile::Open(Form("R%s.root",masspoint[i].data()));
	}
	//double xsec[9]={1.90,0.763,0.33,0.155,7.65e-2,1.58e-2,3.73e-3,2.08e-4,4.91e-5};
		
	double ptBins[14]={200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2500};
	double ptBinsCenter[14]={350,450,550,650,750,850,950,1125,1375,1625,1875};
	double ptBinsCenterE[14]={350,450,550,650,750,850,950,1125,1375,1625,1875};
	double ptBinsError[14]={0};
	double ptBinsErrorE[14]={0};
	
	double mean[6][15];
	double sigma[6][15];
	
	for(int i=0;i<10;i++){
		TH1D* th1=(TH1D	*)tf1[i]->Get(Form("ptBarel%.0f",ptBins[0]));
		for(int j=1;j<14;j++){
			TH1D* th2=(TH1D	*)tf1[i]->Get(Form("ptBarel%.0f",ptBins[j]));
			th1->Add(th2);
		}
		
		ptBinsCenter[i]=th1->GetMean();
		ptBinsError[i]=th1->GetRMS();
		ptBinsError[i]=20;
		cout<<i<<","<<ptBinsCenter[i]<<","<<ptBinsError[i]<<endl;
	}
	
	
	
	for(int i=0;i<10;i++){
		TH1D* th1=(TH1D	*)tf1[i]->Get(Form("ptEndcap%.0f",ptBins[0]));
		for(int j=1;j<14;j++){
			TH1D* th2=(TH1D	*)tf1[i]->Get(Form("ptEndcap%.0f",ptBins[j]));
			th1->Add(th2);
		}
		
		ptBinsCenterE[i]=th1->GetMean();
		ptBinsErrorE[i]=th1->GetRMS();
		ptBinsErrorE[i]=20;
		cout<<i<<","<<ptBinsCenterE[i]<<","<<ptBinsErrorE[i]<<endl;
	}
	
	
	//RooRealVar x("x","x",0,150) ;
	
	for(int i=0;i<10;i++){
		TH1D* th1=(TH1D*)tf1[i]->Get(Form("recoBarelMass%.0f",ptBins[i]));
		for(int j=1;j<14;j++){
			TH1D* th2=(TH1D	*)tf1[i]->Get(Form("recoBarelMass%.0f",ptBins[j]));
			th1->Add(th2);
		}
		th1->SetTitle(Form("%d",masspoint[i]));
		th1->SetMaximum(th1->GetMaximum()*1.3);
		
		/*
		RooDataHist dh("dh","dh",x,Import(*th1)) ;
		RooRealVar mean("mean","mean",125,50,150) ;
		RooRealVar sigmaR("sigmaR","sigmaR",4,0.1,9) ;
		RooRealVar sigmaL("sigmaL","sigmaL",10,0.1,18) ;
		RooBifurGauss Bi("bi","bi",x,mean,sigmaL,sigmaR);
		
		
		RooRealVar m0("m0","m0",80,0,150);
		RooRealVar sig("sigma","sigma",12,0,100);
		RooRealVar alpha("alpha","alpha",1,0,10);
		RooRealVar n("n","n",1,0,50);
		RooCBShape CB("CB","CB",x,m0,sig,alpha,n);
		RooPlot* frame=x.frame(Title(Form("%d",masspoint[i])));
		dh.plotOn(frame);
		CB.fitTo(dh) ; 
		CB.plotOn(frame,LineColor(kRed)) ;
		CB.paramOn(frame,Layout(0.55)) ;
		mean[4][i]=80/m0.getValV();
		sigma[4][i]=m0.getError()*(80/m0.getValV())/m0.getValV();
		//mean[4][i]=125/th1->GetMean();
		//sigma[4][i]=th1->GetMeanError()/th1->GetMean();
		//th1->Draw();
		frame->SetMaximum(th1->GetMaximum()*1.3);
		frame->Draw() ;
		//cout<<i<<"="<<mean[4][i]<<endl;
		*/
		mean[4][i]=80/th1->GetMean();
		sigma[4][i]=th1->GetMeanError()/th1->GetMean();

		TF1* f2 = new TF1("f2","[4]*ROOT::Math::crystalball_function(x,[0],[1],[2],[3])");
		f2->SetParameters(0.15,100,0.2,80,10);
		th1->Fit("f2","","",40,100);
		
	if(i==0)c1->Print("plotsR/recoBarel.pdf(");
		else if(i==9)c1->Print("plotsR/recoBarel.pdf)");
		else  c1->Print("plotsR/recoBarel.pdf");
	}
	
	for(int i=0;i<10;i++){
		
		TH1D* th1=(TH1D*)tf1[i]->Get(Form("recoEndcapMass%.0f",ptBins[i]));
		for(int j=1;j<14;j++){
			TH1D* th2=(TH1D	*)tf1[i]->Get(Form("recoEndcapMass%.0f",ptBins[j]));
			th1->Add(th2);
		}
		th1->SetTitle(Form("%d",masspoint[i]));
		th1->SetMaximum(th1->GetMaximum()*1.3);
		/*
		RooDataHist dh("dh","dh",x,Import(*th1)) ;
		RooRealVar mean("mean","mean",125,50,150) ;
		RooRealVar sigmaR("sigmaR","sigmaR",4,0.1,9) ;
		RooRealVar sigmaL("sigmaL","sigmaL",10,0.1,18) ;
		RooBifurGauss Bi("bi","bi",x,mean,sigmaL,sigmaR);
		
		
		RooRealVar m0("m0","m0",80,0,150);
		RooRealVar sig("sigma","sigma",12,0,100);
		RooRealVar alpha("alpha","alpha",1,0,10);
		RooRealVar n("n","n",1,0,50);
		RooCBShape CB("CB","CB",x,m0,sig,alpha,n);
		RooPlot* frame=x.frame(Title(Form("%d",masspoint[i])));
		dh.plotOn(frame);
		CB.fitTo(dh) ; 
		CB.plotOn(frame,LineColor(kRed)) ;
		CB.paramOn(frame,Layout(0.55)) ;
		mean[5][i]=80/m0.getValV();
		sigma[5][i]=m0.getError()*(80/m0.getValV())/m0.getValV();
		//mean[4][i]=125/th1->GetMean();
		//sigma[4][i]=th1->GetMeanError()/th1->GetMean();
		//th1->Draw();
		frame->SetMaximum(th1->GetMaximum()*1.3);
		frame->Draw() ;
		//cout<<i<<"="<<mean[5][i]<<endl;
		*/
		  mean[5][i]=80/th1->GetMean();
		sigma[5][i]=th1->GetMeanError()/th1->GetMean();

		TF1* f2 = new TF1("f2","[4]*ROOT::Math::crystalball_function(x,[0],[1],[2],[3])");
		f2->SetParameters(0.15,100,0.2,80,10);
                th1->Fit("f2","","",40,100);
		
	if(i==0)c1->Print("plotsR/recoEndcap.pdf(");
		else if(i==9)c1->Print("plotsR/recoEndcap.pdf)");
		else  c1->Print("plotsR/recoEndcap.pdf");
	
	}
	
	
	
	
	TGraphErrors* tg1[6];
	
	tg1[4]=new TGraphErrors(10,ptBinsCenter,mean[4],ptBinsError,sigma[4]);
	tg1[5]=new TGraphErrors(10,ptBinsCenterE,mean[5],ptBinsErrorE,sigma[5]);
	
	for(int i=0;i<14;i++)cout<<i<<"="<<mean[4][i]<<endl;
	for(int i=0;i<14;i++)cout<<i<<"="<<mean[5][i]<<endl;
	
	TLegend *leg = new TLegend(0.68, 0.65, 0.94, 0.90);
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
	
	
	tg1[4]->GetXaxis()->SetTitle("jet Pt");
	tg1[4]->GetYaxis()->SetTitle("M_{PDG}/M_{Reco}");
	tg1[4]->SetTitle("Gen Correction");
	tg1[4]->SetMinimum(1);
	tg1[4]->SetMaximum(1.4);
	tg1[4]->Draw("APL");
	tg1[4]->SetFillColor(0);
	tg1[5]->SetFillColor(0);
	tg1[5]->SetLineColor(2);
	tg1[5]->SetMarkerColor(2);
	tg1[5]->Draw("PLsame");
	
	TFile* outFile = new TFile("pdgToReco2.root","recreate");
	tg1[4]->SetName("barel");
	tg1[5]->SetName("endcap");
	tg1[4]->Write();
	tg1[5]->Write();
	outFile->Close();
	/*
	TF1* recoOneBarel = new TF1("genBarel","[0]+[1]*pow(x*[2],-[3])");
	  recoOneBarel->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
	TF1* recoOneEndcap = new TF1("genEndcap","[0]+[1]*pow(x*[2],-[3])");
	  recoOneEndcap->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
	*/
	
	  TF1* f2 = new TF1("f2","[0]/pow(x,0.5)+[1]/pow(x,[7])+[2]*x+[3]+[4]/x+[5]/pow(x,[6])");

f2->SetParameters(
-1356.74,
12127.2,
-0.00129311,
15.6148,
-18669.2,
5.14592e+4,
6.50081,
0.7995);

TF1* f3 =  new TF1("fendcap","[0]+[1]/sqrt(x)+[2]/x+[3]/pow(x,0.75)");
/*f3->SetParameters(
-26.9944,
434.986,
1.62699);*/
	
	f2->SetLineColor(1);
	f3->SetLineColor(2);
	//tg1[4]->Fit(f2,"","",200,2000);
	//tg1[4]->Fit(f2,"","",200,2000);
	//tg1[4]->Fit(f2,"","",200,2000);
	//tg1[5]->Fit(f3,"","",200,2000);
	//tg1[5]->Fit(f3,"","",200,2000);
leg->Clear();


 TF1* puppisd_corrGEN      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGEN->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
  TF1* puppisd_corrRECO_cen = new TF1("puppisd_corrRECO_cen","([0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5))*([6]+[7]*pow(x*[8],-[9]))",200,2000);
  puppisd_corrRECO_cen->SetParameters(
   				      1.05807,
   				      -5.91971e-05,
   				      2.296e-07,
   				      -1.98795e-10,
   				      6.67382e-14,
   				      -7.80604e-18,
					 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				      );

  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","([0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5))*([6]+[7]*pow(x*[8],-[9]))",200,2000);
  puppisd_corrRECO_for->SetParameters(
   				      1.26638,
   				      -0.000658496,
   				      9.73779e-07,
   				      -5.93843e-10,
   				      1.61619e-13,
   				      -1.6272e-17,
						 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454);
					
					
				TF1* n=new TF1("n","puppisd_corrRECO_cen*puppisd_corrGEN");	
				TF1* n2=new TF1("n","puppisd_corrRECO_for*puppisd_corrGEN");	
  
  leg->AddEntry(tg1[4],"reco barel");
  leg->AddEntry(tg1[5],"reco endcap");
  leg->AddEntry(puppisd_corrRECO_cen,"Thea barel");
  leg->AddEntry(puppisd_corrRECO_for,"Thea endcap");

	leg->Draw("same");
	//genBarel->Draw("same");
//	genEndcap->Draw("same");
	puppisd_corrRECO_cen->SetLineColor(3);
	puppisd_corrRECO_for->SetLineColor(4);
	//tg1[2]->Draw("APL");
	puppisd_corrRECO_cen->Draw("same");
	puppisd_corrRECO_for->Draw("same");
	
	c1->Print("plotsR/recoOne.pdf");
	
}
