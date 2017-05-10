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
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "setNCUStyle.C"


using namespace std;
void dbtsf(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());

  
  TH1F* h_Mjjred1 =  new TH1F("","",84,800,5000);
  TH1F* h_Mjjred2 =  new TH1F("","",84,800,5000);
  TH1F* h_Mjjred3 =  new TH1F("","",84,800,5000);

  TH1F* h_leadPt=new TH1F("","",40,200,1000);
  TH1F* h_sublPt=new TH1F("","",40,200,1000);

  Float_t nPass[20]={0};
  
  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    
    nPass[0]++;

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);
    
    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCSV    = data.GetPtrFloat("FATjetCSV");
    Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Int_t*   FATnSubSDJet   = data.GetPtrInt("FATnSubSDJet");
    vector<float>* FATsubjetSDCSV       = data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);
    Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
    Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
    vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));    

    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doubleSV = data.GetPtrFloat("ADDjet_DoubleSV");

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
	bool results = trigResult[it];

	// std::cout << thisTrig << " : " << results << std::endl;
	
	if( (thisTrig.find("HLT_PFHT800")!= std::string::npos && results==1)
	    )
	  {
	    passTrigger=true;
	    break;
	  }


      }


    if(!passTrigger)continue;

    nPass[1]++;   

    //3. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
       
    vector<int> fatjet;
    vector<pair<int,int>> Mjj;
    for(int ij=0; ij<nFJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      if(thisJet->Pt()<200)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!FATjetPassIDTight[ij])continue;
    
      
      fatjet.push_back(ij);	
    }
    
    if(fatjet.size()<2)continue;

    nPass[2]++;

    int aa = fatjet[0]; //Mjj[0].second;
    int ee = fatjet[1]; //Mjj[0].first;
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa); 
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);

    Double_t dEta = fabs(Jet1->Eta() - Jet2->Eta());
    if(dEta>1.3)continue;

    nPass[3]++;

    Float_t mff=(*Jet1+*Jet2).M();
    Float_t msubt = mff-(fatjetPRmassL2L3Corr[aa]-125)-(fatjetPRmassL2L3Corr[ee]-125);
    if(msubt<800)continue;

    nPass[4]++;

    if(fatjetPRmassL2L3Corr[aa]<105 || fatjetPRmassL2L3Corr[aa]>135)continue;
    if(fatjetPRmassL2L3Corr[ee]<105 || fatjetPRmassL2L3Corr[ee]>135)continue;

    nPass[5]++;

    Double_t leadtau21 = (fatjetTau2[aa]/fatjetTau1[aa]);
    Double_t subltau21 = (fatjetTau2[ee]/fatjetTau1[ee]);
    if(leadtau21>0.6)continue;
    if(subltau21>0.6)continue;

    nPass[6]++;

    int addJetIndex[2]={-1,-1}; 
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;

    if(addjet_doubleSV[addJetIndex[0]]<0.6)continue;
    if(addjet_doubleSV[addJetIndex[1]]<0.6)continue;

    nPass[7]++;

    h_leadPt->Fill(Jet1->Pt());
    h_sublPt->Fill(Jet2->Pt());

    if(Jet1->Pt()<400 || Jet2->Pt()<400) {
      //0.92 +/- 0.09
      h_Mjjred1->Fill(msubt,0.92);
      h_Mjjred2->Fill(msubt,1.01);
      h_Mjjred3->Fill(msubt,0.83);
    }
   
    if((Jet1->Pt()>400 && Jet1->Pt()<500) || (Jet2->Pt()>400 && Jet2->Pt()<500)) {
      //0.99 +/- 0.12
      h_Mjjred1->Fill(msubt,0.99);
      h_Mjjred2->Fill(msubt,1.11);
      h_Mjjred3->Fill(msubt,0.87);
    }
    
    if((Jet1->Pt()>500 && Jet1->Pt()<600) || (Jet2->Pt()>500 && Jet2->Pt()<600)) {
      //0.94 +/- 0.19
      h_Mjjred1->Fill(msubt,0.94);
      h_Mjjred2->Fill(msubt,1.13);
      h_Mjjred3->Fill(msubt,0.75);
    }

    if((Jet1->Pt()>600) || (Jet2->Pt()>600)) {
      //1.05 +/- 0.21
      h_Mjjred1->Fill(msubt,1.05);
      h_Mjjred2->Fill(msubt,1.26);
      h_Mjjred3->Fill(msubt,0.84);
    }

  } //end of the event loop

  setNCUStyle();


  std::string bulkg_name[]={"1000","1200","1400","1600","1800","2000","2500","3000","4000","4500"};
  std::string radion_name[]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4500"};
  float bmass[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
  float rmass[10]={1000,1200,1400,1600,1800,2000,2500,3000,3500,4500};

  int nBM = sizeof(bmass)/sizeof(bmass[0]);
  int nRM = sizeof(rmass)/sizeof(rmass[0]);

  h_Mjjred1->Scale((2.6837)/nPass[0]);
  h_Mjjred2->Scale((2.6837)/nPass[0]);
  h_Mjjred3->Scale((2.6837)/nPass[0]);

  TCanvas* c2 =  new TCanvas("c2","c2",0,0,800,600);
  c2->cd();
  h_Mjjred1->SetLineWidth(2);
  h_Mjjred2->SetLineWidth(2);
  h_Mjjred3->SetLineWidth(2);
  h_Mjjred1->SetLineStyle(10);
  h_Mjjred2->SetLineStyle(1);
  h_Mjjred3->SetLineStyle(5);
  h_Mjjred2->Draw("hist");
  h_Mjjred2->GetYaxis()->SetTitle("");
  h_Mjjred2->GetXaxis()->SetTitle("reduced dijet mass");
  h_Mjjred2->GetXaxis()->SetTitleSize(0.04);
  h_Mjjred2->GetYaxis()->SetTitleSize(0.04);
  h_Mjjred2->SetLineColor(kOrange+2);
  h_Mjjred1->SetLineColor(kViolet+2);
  h_Mjjred1->Draw("histsame");
  h_Mjjred3->SetLineColor(kGreen+2);
  h_Mjjred3->Draw("histsame");
  TLegend *leg1 = new TLegend(0.52, 0.72, 0.90, 0.88);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(h_Mjjred1, "central", "pl");
  leg1->AddEntry(h_Mjjred2, "up", "pl");
  leg1->AddEntry(h_Mjjred3, "down", "pl");
  leg1->Draw();

  Double_t central1 = h_Mjjred1->Integral();
  Double_t up1      = h_Mjjred2->Integral();
  Double_t down1    = h_Mjjred3->Integral();

  float sysmjjred = abs(up1-down1)/(2*central1);
  float sys1up = abs(up1-central1)/central1;
  float sys1down = abs(down1-central1)/central1;

  cout<<sysmjjred<<" "<<sys1up<<" "<<sys1down<<" "<<endl;

  TCanvas* c =  new TCanvas("c","c",0,0,800,600);
  h_leadPt->SetLineWidth(2);
  h_sublPt->SetLineWidth(2);
  h_leadPt->SetLineColor(kOrange+2);
  h_sublPt->SetLineColor(kViolet+2);
  h_leadPt->Draw("hist");
  h_sublPt->Draw("histsame");

}
