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
#include "readSample.h"
 
using namespace std;
void backgroundestimate(std::string inputFile) {

  int total = 0;
  //read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
    
  TreeReader data(infiles);
  TH1F* h_dsv[2][4];
  TH1F* h_prmass[2][4];
  TH1F* h_pt[2][4];
  TH1F* h_eta[2][4];
  TH1F* h_tau21[2][4];

  TH1F* h_mjj[4];
  TH1F* h_mjjred[4];
  TH1F* h_deleta[4];

  for(int i=0;i<2;i++) {
    for(int j=0;j<4;j++) {
      h_dsv[i][j] = new TH1F("","",20,-1,1);
      h_prmass[i][j] = new TH1F("","",25,90,140);
      h_pt[i][j] = new TH1F("","",40,200,1000);
      h_eta[i][j] = new TH1F("","",50,-2.5,2.5);
      h_tau21[i][j] = new TH1F("","",20,0,1);

      h_dsv[i][j]->Sumw2();
      h_prmass[i][j]->Sumw2();
      h_pt[i][j]->Sumw2();
      h_eta[i][j]->Sumw2();
      h_tau21[i][j]->Sumw2();
    }
  }

  for(int i=0;i<4;i++) {
    h_mjjred[i] = new TH1F("","",84,800,5000);
    h_mjj[i] = new TH1F("","",32,900,2500);
    h_deleta[i] = new TH1F("","",28,0,1.4);

    h_mjjred[i]->Sumw2();
    h_mjj[i]->Sumw2();
    h_deleta[i]->Sumw2();
  }

  //------------------method 1--------------------------------
  float fa[3] = {8.11492e-02,(8.11492e-02+1.12985e-02),(8.11492e-02-1.12985e-02)};
  float fb[3] = {3.40934e-04,(3.40934e-04+1.62748e-04),(3.40934e-04-1.62748e-04)};
  float fc[3] = {8.25267e-06,(8.25267e-06+4.98605e-06),(8.25267e-06-4.98605e-06)};
  //-----------------------------------------------------------

  //------------------method 2-----------------------------
  int nbin =75;
  double prunedm[nbin];
  double pfratio[3][nbin];
  //double pfup[nbin];
  //double pfdw[nbin];

  ifstream fin;
  fin.open("bkgestsysunc_binbybin.dat");
  for(int i=0;i<75;i++) {
    fin >> prunedm[i] >> pfratio[0][i] >> pfratio[1][i] >> pfratio[2][i];
  }
  fin.close();
  //-------------------------------------------------------
  
  TProfile* p_PRvW = new TProfile("","",25,90,140,0,1);
  TProfile* p_MRedvW = new TProfile("","",36,700,2500,0,1);
  /*
  float xbins1[] = {-75,-45,-20,10,40,75};
  const int nbins1=sizeof(xbins1)/sizeof(xbins1[0])-1;
  TH1F* hpass = new TH1F("","",nbins1,xbins1);
  TH1F* hfail = new TH1F("","",nbins1,xbins1);
  TH1F* hratio = new TH1F("","",nbins1,xbins1);
  hratio->Sumw2();
  hpass->Sumw2();
  hfail->Sumw2();
  */
  Long64_t nPass[20] ={0};
  total += data.GetEntriesFast();
  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);
    
    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCSV    = data.GetPtrFloat("FATjetCSV");
    Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
    Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Int_t*   FATnSubSDJet   = data.GetPtrInt("FATnSubSDJet");
    vector<float>* FATsubjetSDCSV       = data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);
    Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
    Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
    Float_t  mcWeight  = data.GetFloat("mcWeight");
    Float_t HT         = data.GetFloat("HT");
    //Float_t*  fatjet_doubleSV = data.GetPtrFloat("FATjet_DoubleSV");
    
    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doubleSV = data.GetPtrFloat("ADDjet_DoubleSV");
    

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    nPass[0]++;

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
      //if(fatjetPRmassL2L3Corr[ij]>105 && fatjetPRmassL2L3Corr[ij]<135)continue;
      if(fatjetPRmassL2L3Corr[ij]<105 || fatjetPRmassL2L3Corr[ij]>135)continue;
      
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

    Double_t leadtau21 = (fatjetTau2[aa]/fatjetTau1[aa]);
    Double_t subltau21 = (fatjetTau2[ee]/fatjetTau1[ee]);
    if(leadtau21>0.6)continue;
    if(subltau21>0.6)continue;

    nPass[5]++;
    
    int addJetIndex[2]={-1,-1}; 
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;


    Double_t weight;

     
    if((addjet_doubleSV[addJetIndex[0]]>0.6) && (addjet_doubleSV[addJetIndex[1]]>0.6)) {
      h_dsv[0][0]->Fill(addjet_doubleSV[addJetIndex[0]]);
      h_dsv[1][0]->Fill(addjet_doubleSV[addJetIndex[1]]);
      h_prmass[0][0]->Fill(fatjetPRmassL2L3Corr[aa]);
      h_prmass[1][0]->Fill(fatjetPRmassL2L3Corr[ee]);
      h_pt[0][0]->Fill(Jet1->Pt());
      h_pt[1][0]->Fill(Jet2->Pt());
      h_eta[0][0]->Fill(Jet1->Eta());
      h_eta[1][0]->Fill(Jet2->Eta());
      h_tau21[0][0]->Fill(leadtau21);
      h_tau21[1][0]->Fill(subltau21);
      h_mjj[0]->Fill(mff);
      h_mjjred[0]->Fill(msubt);
      h_deleta[0]->Fill(dEta);

      nPass[6]++;
    }
    //one anti-tag LFSP
    
    else if((addjet_doubleSV[addJetIndex[0]]<0.6) && (addjet_doubleSV[addJetIndex[1]]>0.6)) {
      //method2
      for(int i=0;i<3;i++) {
	for(int j=0;j<nbin;j++) {
	  if(h_prmass[0][i+1]->FindBin(fatjetPRmassL2L3Corr[aa]) == h_prmass[0][i+1]->FindBin(prunedm[j])) {
	    h_dsv[0][i+1]->Fill(addjet_doubleSV[addJetIndex[0]],pfratio[i][j]);
	    h_dsv[1][i+1]->Fill(addjet_doubleSV[addJetIndex[1]],pfratio[i][j]);
	    h_prmass[0][i+1]->Fill(fatjetPRmassL2L3Corr[aa],pfratio[i][j]);
	    h_prmass[1][i+1]->Fill(fatjetPRmassL2L3Corr[ee],pfratio[i][j]);
	    h_pt[0][i+1]->Fill(Jet1->Pt(),pfratio[i][j]);
	    h_pt[1][i+1]->Fill(Jet2->Pt(),pfratio[i][j]);
	    h_eta[0][i+1]->Fill(Jet1->Eta(),pfratio[i][j]);
	    h_eta[1][i+1]->Fill(Jet2->Eta(),pfratio[i][j]);
	    h_tau21[0][i+1]->Fill(leadtau21,pfratio[i][j]);
	    h_tau21[1][i+1]->Fill(subltau21,pfratio[i][j]);
	    h_mjj[i+1]->Fill(mff,pfratio[i][j]);
	    h_mjjred[i+1]->Fill(msubt,pfratio[i][j]);
	    h_deleta[i+1]->Fill(dEta,pfratio[i][j]);
	  }
	}
      }
      
      /* //method1
      for(int i=0;i<3;i++) {
	cout<<fa[i]<<" "<<fb[i]<<" "<<fc[i]<<endl;
	weight = (fa[i])+((fb[i])*(fatjetPRmassL2L3Corr[aa]-125))+((fc[i])*(fatjetPRmassL2L3Corr[aa]-125)*(fatjetPRmassL2L3Corr[aa]-125));
	cout<<weight<<endl;
	cout<<"2.1.1"<<endl;
	h_dsv[0][i+1]->Fill(addjet_doubleSV[addJetIndex[0]],weight);
	h_dsv[1][i+1]->Fill(addjet_doubleSV[addJetIndex[1]],weight);
	cout<<"2.1.2"<<endl;
	h_prmass[0][i+1]->Fill(fatjetPRmassL2L3Corr[aa],weight);
	h_prmass[1][i+1]->Fill(fatjetPRmassL2L3Corr[ee],weight);
	cout<<"2.1.3"<<endl;
	h_pt[0][i+1]->Fill(Jet1->Pt(),weight);
	h_pt[1][i+1]->Fill(Jet2->Pt(),weight);
	cout<<"2.1.4"<<endl;
	h_eta[0][i+1]->Fill(Jet1->Eta(),weight);
	h_eta[1][i+1]->Fill(Jet2->Eta(),weight);
	cout<<"2.1.5"<<endl;
	h_tau21[0][i+1]->Fill(leadtau21,weight);
	h_tau21[1][i+1]->Fill(subltau21,weight);
	cout<<"2.1.6"<<endl;
	h_mjj[i+1]->Fill(mff,weight);
	h_mjjred[i+1]->Fill(msubt,weight);
	h_deleta[i+1]->Fill(dEta,weight);
      }
      */
    }

  } //end of the event loop

  for(int i=0;i<20;i++) {
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;
  }
   //hratio->Divide(hpass,hfail); 
  TFile* outfile = new TFile("JetHT-Run2015_bkgest.10.root","recreate");
  for(int i=0;i<2;i++) {
    for(int j=0;j<4;j++) {
      h_dsv[i][j]->Write(Form("DBT%i_%i",i,j));
      h_prmass[i][j]->Write(Form("PR%i_%i",i,j));
      h_pt[i][j]->Write(Form("Pt%i_%i",i,j));
      h_eta[i][j]->Write(Form("Eta%i_%i",i,j));
      h_tau21[i][j]->Write(Form("Tau21%i_%i",i,j));
    }
  }

  for(int i=0;i<4;i++) {
    h_mjj[i]->Write(Form("Mjj_%i",i));
    h_mjjred[i]->Write(Form("MjjRed_%i",i));
    h_deleta[i]->Write(Form("DelEta_%i",i));
  }
   outfile->Write();
 
}
