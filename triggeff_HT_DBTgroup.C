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
#include <TStyle.h>
#include "untuplizer_dbtgroup.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "readSample.h"
#include "setNCUStyle.C"
 
using namespace std;
void triggeff_HT_DBTgroup(std::string inputFile) {

  int total = 0;

  std::vector<string> infiles;

  readSample(inputFile, infiles);
    
  TreeReader data(infiles);

  TH1F* h_Mjj1 =  new TH1F("","",80,700,3000);
  TH1F* h_Mjj2 =  new TH1F("","",80,700,3000);
  TH1F* h_Mjj3 =  new TH1F("","",80,700,3000);
  TH1F* h_Mjjred1 =  new TH1F("","",80,700,3000);
  TH1F* h_Mjjred2 =  new TH1F("","",80,700,3000);
  TH1F* h_Mjjred3 =  new TH1F("","",80,700,3000);
  h_Mjj1->Sumw2();
  h_Mjj2->Sumw2();
  h_Mjj3->Sumw2();
  h_Mjjred1->Sumw2();
  h_Mjjred2->Sumw2();
  h_Mjjred3->Sumw2();

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

    if(nPVs<1)continue;

    vector<int> fatjet;
    for(int ij=0;ij<nFJets; ij++) {
      //if(FatjetPt[ij]<200)continue;
      //if(fabs(FatjetEta[ij])>2.4)continue;
      //if(FatjetIDTight[ij] == 0)continue;
      //if(FatjetPRmassCorr[ij]<50 || FatjetPRmassCorr[ij]>200)continue;
       
      fatjet.push_back(ij);
    }
  
    if(fatjet.size()<2)continue;

    int lead = fatjet[0];
    int subl = fatjet[1];

    TLorentzVector Jet1(0,0,0,0);
    TLorentzVector Jet2(0,0,0,0);
    Jet1.SetPtEtaPhiM(FatjetPt[lead],FatjetEta[lead],FatjetPhi[lead],FatjetMass[lead]);
    Jet2.SetPtEtaPhiM(FatjetPt[subl],FatjetEta[subl],FatjetPhi[subl],FatjetMass[subl]);

    TLorentzVector Mjj(0,0,0,0);
    Mjj = Jet1 + Jet2;
    Float_t msubt = Mjj.M()-(FatjetPRmassCorr[lead]-125)-(FatjetPRmassCorr[subl]-125);
    //if(msubt<800)continue;

    if(passedHLT350==1 && passedHLT800==1) {
      h_Mjj1->Fill(Mjj.M());
      h_Mjjred1->Fill(msubt);
      }

    if(passedHLT350 == 1) {
      h_Mjj2->Fill(Mjj.M());
      h_Mjjred2->Fill(msubt);
    }

  } //end of the event loop

  setNCUStyle(true);
  TCanvas* c3 =  new TCanvas("c3","c3",0,0,800,600);
  c3->cd();
  //h_Mjj3->GetXaxis()->SetRangeUser(400,5000);
  h_Mjj3->GetYaxis()->SetTitle("HT350&&HT800/HT350");                                                                                                   
  h_Mjj3->GetXaxis()->SetTitle("dijet mass");
  h_Mjj3->GetXaxis()->SetTitleSize(0.05);
  h_Mjj3->GetYaxis()->SetTitleSize(0.05);

  h_Mjj3->SetLineWidth(2);
  h_Mjj3->SetLineColor(kRed+3);
  h_Mjj3->Divide(h_Mjj1,h_Mjj2);
  h_Mjj3->Draw();

  c3->Print("TriggerEff_Mjj_JetHT-Run2015_Prescaledx.pdf");

  TCanvas* c =  new TCanvas("c","c",0,0,800,600);
  c->cd();
  //h_Mjjred3->GetXaxis()->SetRangeUser(400,5000);
  h_Mjjred3->GetYaxis()->SetTitle("HT350&&HT800/HT350");
  h_Mjjred3->GetXaxis()->SetTitle("reduced dijet mass");
  h_Mjjred3->GetXaxis()->SetTitleSize(0.05);
  h_Mjjred3->GetYaxis()->SetTitleSize(0.05);

  h_Mjjred3->SetLineWidth(2);
  h_Mjjred3->SetLineColor(kGreen+3);
  h_Mjjred3->Divide(h_Mjjred1,h_Mjjred2);
  h_Mjjred3->Draw();

  c->Print("TriggerEff_Mjjred_JetHT-Run2015_Prescaledx.pdf");

  float stddev_mjj = h_Mjj3->GetStdDev();
  float stddev_mjjred = h_Mjjred3->GetStdDev();
  cout<<stddev_mjj<<" "<<stddev_mjjred<<endl;
}
