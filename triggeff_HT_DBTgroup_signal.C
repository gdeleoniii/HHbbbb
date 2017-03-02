#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TPad.h>
#include <THnSparse.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include "untuplizer_dbtgroup.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "readSample.h"
#include "setNCUStyle.C"
 
using namespace std;
void triggeff_HT_DBTgroup_signal(std::string inputFile) {

  int total = 0;

  std::vector<string> infiles;

  readSample(inputFile, infiles);
    
  TreeReader data(infiles);

  TGraphAsymmErrors* hAE_Mjj = new TGraphAsymmErrors();
  TGraphAsymmErrors* hAE_Mjjred = new TGraphAsymmErrors();

  TH1F* h_Mjj1 =  new TH1F("","",44,800,3000);
  TH1F* h_Mjj2 =  new TH1F("","",44,800,3000);
  TH1F* h_Mjj3 =  new TH1F("","",44,800,3000);
  TH1F* h_Mjjred1 =  new TH1F("","",44,800,3000);
  TH1F* h_Mjjred2 =  new TH1F("","",44,800,3000);
  TH1F* h_Mjjred3 =  new TH1F("","",44,800,3000);
  h_Mjj1->Sumw2();
  h_Mjj2->Sumw2();
  h_Mjj3->Sumw2();
  h_Mjjred1->Sumw2();
  h_Mjjred2->Sumw2();
  h_Mjjred3->Sumw2();


  //----------------------------------
  Int_t nbins = 44;
  double mjj_central[nbins];
  double mjj_up[nbins];
  double mjj_down[nbins];
  double mjjred_central[nbins];
  double mjjred_up[nbins];
  double mjjred_down[nbins];

  float mjj_mass[nbins];
  float mjjred_mass[nbins];

  ifstream fin;
  fin.open("Mjj_TriggerEff_Weights.dat");
  ifstream fin1;
  fin1.open("Mjjred_TriggerEff_Weights.dat");
  for(int i = 0;i <nbins;i++) {
    //ifstream fin;
    //fin.open("Mjj_TriggerEff_Weights.dat");
    fin >> mjj_mass[i] >> mjj_central[i] >> mjj_up[i]>> mjj_down[i];
    //fin.close();

    //ifstream fin1;
    //fin1.open("Mjjred_TriggerEff_Weights.dat");
    fin1 >> mjjred_mass[i] >> mjjred_central[i] >> mjjred_up[i]>> mjjred_down[i];
    //fin1.close();
  }
  fin.close();
  fin1.close();
  //----------------------------------

  Int_t x = 0;
  total += data.GetEntriesFast();

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);

    int nFATJet         = data.GetInt("nFatjetAK08ungroomed");
    Float_t nPVs = data.GetFloat("nPVs");
    const int nFJets=nFATJet;
    Float_t*  FatjetMass = data.GetPtrFloat("FatjetAK08ungroomed_mass");
    Float_t*  FatjetPRmass = data.GetPtrFloat("FatjetAK08ungroomed_mpruned");
    Float_t*  FatjetPt = data.GetPtrFloat("FatjetAK08ungroomed_pt");
    Float_t*  FatjetEta = data.GetPtrFloat("FatjetAK08ungroomed_eta");
    Float_t*  FatjetPhi = data.GetPtrFloat("FatjetAK08ungroomed_phi");
    Float_t*  FatjetIDTight = data.GetPtrFloat("FatjetAK08ungroomed_id_Tight");
    Float_t*  FatjetPRmassCorr = data.GetPtrFloat("FatjetAK08ungroomed_mprunedcorr");
    Float_t*  FatjetTau1 = data.GetPtrFloat("FatjetAK08ungroomed_tau1");
    Float_t*  FatjetTau2 = data.GetPtrFloat("FatjetAK08ungroomed_tau2");
    Float_t*  FatjetDBT = data.GetPtrFloat("FatjetAK08ungroomed_bbtag");
    int passedHLT800 = data.GetInt("HLT_BIT_HLT_PFHT800_v");
    int passedHLT350 = data.GetInt("HLT_BIT_HLT_PFHT350_v");
    //if(passedHLT == 0)continue;

    //if(nPVs<1)continue;

    vector<int> fatjet;
    for(int ij=0;ij<nFJets; ij++) {
      if(FatjetPt[ij]<200)continue;
      if(fabs(FatjetEta[ij])>2.4)continue;
      if(FatjetIDTight[ij] == 0)continue;
      if(FatjetPRmassCorr[ij]<50)continue;
       
      fatjet.push_back(ij);
    }
  
    if(fatjet.size()<2)continue;

    int lead = fatjet[0];
    int subl = fatjet[1];

    TLorentzVector Jet1(0,0,0,0);
    TLorentzVector Jet2(0,0,0,0);
    Jet1.SetPtEtaPhiM(FatjetPt[lead],FatjetEta[lead],FatjetPhi[lead],FatjetMass[lead]);
    Jet2.SetPtEtaPhiM(FatjetPt[subl],FatjetEta[subl],FatjetPhi[subl],FatjetMass[subl]);

    Double_t DelEta = fabs(Jet1.Eta() - Jet2.Eta());
    if(DelEta>1.3)continue;

    //TLorentzVector Mjj(0,0,0,0);
    Double_t Mjj = (Jet1 + Jet2).M();
    Double_t msubt = Mjj-(FatjetPRmassCorr[lead]-125)-(FatjetPRmassCorr[subl]-125);

    /*
    cout << " " <<passedHLT800 << " " << passedHLT350 << endl;
    if(passedHLT800 == 0) continue;
    h_Mjj1->Fill(Mjj);
    h_Mjjred1->Fill(msubt);

    if(passedHLT350 == 0) continue;
    h_Mjj2->Fill(Mjj);
    h_Mjjred2->Fill(msubt);
    //x+=1;
    */
    

    Float_t width = 50;
    cout<<Mjj<<endl;
    for(int i = 0;i <nbins;i++) {
      
      if(Mjj > mjj_mass[i] && Mjj < (mjj_mass[i]+width)) {
	h_Mjj1->Fill(Mjj,mjj_central[i]);
	h_Mjj2->Fill(Mjj,mjj_up[i]);
	h_Mjj3->Fill(Mjj,mjj_down[i]);
	cout<<i<<" "<<Mjj<<" "<<mjj_mass[i]<<" "<<mjj_mass[i]+width<<endl;
      }
      
      if(msubt > mjjred_mass[i] && msubt < (mjjred_mass[i]+width)) {
	h_Mjjred1->Fill(msubt,mjjred_central[i]);
	h_Mjjred2->Fill(msubt,mjjred_up[i]);
	h_Mjjred3->Fill(msubt,mjjred_down[i]);
      }

      
    }
    
    
    
  } //end of the event loop
  
  setNCUStyle(true);
  TCanvas* c3 =  new TCanvas("c3","c3",0,0,800,600);
  c3->cd();
  h_Mjj1->SetLineWidth(2);
  h_Mjj2->SetLineWidth(2);
  h_Mjj3->SetLineWidth(2);
  h_Mjj2->Draw("hist");
  h_Mjj2->GetYaxis()->SetTitle("");                                                                                                   
  h_Mjj2->GetXaxis()->SetTitle("dijet mass");
  h_Mjj2->GetXaxis()->SetTitleSize(0.04);
  h_Mjj2->GetYaxis()->SetTitleSize(0.04);
  h_Mjj2->SetLineColor(kBlue-1);
  h_Mjj1->SetLineColor(kRed-1);
  h_Mjj1->Draw("histsame");
  h_Mjj3->SetLineColor(kGreen-1);
  h_Mjj3->Draw("histsame");
  TLegend *leg = new TLegend(0.67, 0.72, 0.79, 0.88);                                                                                                        
  leg->SetBorderSize(0);                                                                                                                                     
  leg->SetFillColor(0);                                                                                                                                    
  leg->SetFillStyle(0);                                                                                                                                   
  leg->SetTextSize(0.03);                                                                                                                                 
  leg->AddEntry(h_Mjj1, "central", "l");                                                                                                                     
  leg->AddEntry(h_Mjj2, "up", "l");                                                                                                                       
  leg->AddEntry(h_Mjj3, "down", "l");                                                                                                              
  leg->Draw();

 

  TCanvas* c2 =  new TCanvas("c2","c2",0,0,800,600);
  c2->cd();
  h_Mjjred1->SetLineWidth(2);
  h_Mjjred2->SetLineWidth(2);
  h_Mjjred3->SetLineWidth(2);
  h_Mjjred2->Draw("hist");
  h_Mjjred2->GetYaxis()->SetTitle("");
  h_Mjjred2->GetXaxis()->SetTitle("dijet mass");
  h_Mjjred2->GetXaxis()->SetTitleSize(0.04);
  h_Mjjred2->GetYaxis()->SetTitleSize(0.04);
  h_Mjjred2->SetLineColor(kBlue-1);
  h_Mjjred1->SetLineColor(kRed-1);
  h_Mjjred1->Draw("histsame");
  h_Mjjred3->SetLineColor(kGreen-1);
  h_Mjjred3->Draw("histsame");
  TLegend *leg1 = new TLegend(0.67, 0.72, 0.79, 0.88);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(h_Mjjred1, "central", "l");
  leg1->AddEntry(h_Mjjred2, "up", "l");
  leg1->AddEntry(h_Mjjred3, "down", "l");
  leg1->Draw();

  Double_t plus = h_Mjj2->GetBinContent(3);
  Double_t minus = h_Mjj3->GetBinContent(3);
  Double_t center = h_Mjj1->GetBinContent(3);

  Double_t plus1 = h_Mjjred2->GetBinContent(3);
  Double_t minus1 = h_Mjjred3->GetBinContent(3);
  Double_t center1 = h_Mjjred1->GetBinContent(3);

  Float_t up = abs(center - plus)/center;
  Float_t dw = abs(center - minus)/center;

  Float_t up1 = abs(center1 - plus1)/center1;
  Float_t dw1 = abs(center1 - minus1)/center1;

  cout<<plus<<" "<<minus<<" "<<center<<" "<<up<< " "<<dw<<endl; 
  cout<<plus1<<" "<<minus1<<" "<<center1<<" "<<up1<< " "<<dw1<<endl;

  std::string bulkg_name[]={"900","1000","1400","1600","2000","2500","3000","4000"};

  for(int i=0;i<8;i++){
    bool bulkmass=(inputFile.find(Form("%s",bulkg_name[i].data()))!= std::string::npos); 
    if(bulkmass) {
      c3->Print(Form("TriggerEff_SysUnc_Mjj_BulkGrav_%s.pdf",bulkg_name[i].data()));
      c2->Print(Form("TriggerEff_SysUnc_Mjjred_BulkGrav_%s.pdf",bulkg_name[i].data()));
    }
  }

}
