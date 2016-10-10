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
void passfailratio(std::string inputFile, std::string name) {

  int total = 0;
  //read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
    
  TreeReader data(infiles);
  TH1F* h1=new TH1F("","",36,700,2500);    
  TH1F* h2=new TH1F("","",36,700,2500);
  TH1F* h3=new TH1F("","",36,700,2500);
  TH1F* h4=new TH1F("","",36,700,2500);
  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();
  h4->Sumw2();

  float xbins1[] = {-75,-45,-20,10,40,75};
  const int nbins1=sizeof(xbins1)/sizeof(xbins1[0])-1;
  TH1F* hpass = new TH1F("","",nbins1,xbins1);
  TH1F* hfail = new TH1F("","",nbins1,xbins1);
  TH1F* hratio = new TH1F("","",nbins1,xbins1);
  hratio->Sumw2();
  hpass->Sumw2();
  hfail->Sumw2();

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
      if(fatjetPRmassL2L3Corr[ij]<60 || fatjetPRmassL2L3Corr[ij]>200)continue;
      
      fatjet.push_back(ij);     
    }
    
    if(fatjet.size()<2)continue;
    
    
    for(unsigned int i=0; i<fatjet.size(); i++) {
      for(unsigned int j=0; j<i; j++) {
        int index_that = fatjet[i];
        int index_those = fatjet[j];
        TLorentzVector* thatJet  = (TLorentzVector*)fatjetP4->At(index_that);
        TLorentzVector* thoseJet = (TLorentzVector*)fatjetP4->At(index_those);
        float dEta = fabs(thatJet->Eta() - thoseJet->Eta());
        if(dEta>1.3)continue;
        

	Mjj.push_back(make_pair(index_that,index_those));

      }
    }

    
    if(Mjj.size()<1)continue;   

    int aa = Mjj[0].second;
    int ee = Mjj[0].first;
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa); 
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);
  
    
    Float_t mff=(*Jet1+*Jet2).M();

    Float_t msubt = mff-(fatjetPRmass[aa]-125)-(fatjetPRmass[ee]-125);
    if(msubt<800)continue;
    Double_t leadtau21 = (fatjetTau2[aa]/fatjetTau1[aa]);
    Double_t subltau21 = (fatjetTau2[ee]/fatjetTau1[ee]);

    if(leadtau21>0.6)continue;
    if(subltau21>0.6)continue;
    
    int addJetIndex[2]={-1,-1}; 
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;


    //-----------for  pass-fail ratio-----------
    Double_t pass_dbt = 0;
    Double_t fail_dbt = 0;
    Double_t ratio_dbt = 0;
    if(fatjetPRmass[aa]<105 || fatjetPRmass[aa]>135) {
      if(addjet_doubleSV[addJetIndex[0]]>0.6) {
	pass_dbt = fatjetPRmass[aa];
	hpass->Fill(fatjetPRmass[aa]-125);
      }
      if(addjet_doubleSV[addJetIndex[0]]<0.6) {
	fail_dbt = fatjetPRmass[aa];
	hfail->Fill(fatjetPRmass[aa]-125);
      }
    }
    ratio_dbt = pass_dbt/fail_dbt;
    //------------------------------------------


    if((fatjetPRmass[aa]<105 || fatjetPRmass[aa]>135) && (fatjetPRmass[ee]<105 || fatjetPRmass[ee]>135))continue;

    Double_t weight;
    
   
    //one anti-tag lead pass and subl fails
    if(addjet_doubleSV[addJetIndex[0]]>0.6 && addjet_doubleSV[addJetIndex[1]]<0.6) {
      h1->Fill(msubt);
    }
    //two anti-tag
    else if(addjet_doubleSV[addJetIndex[0]]<0.6 && addjet_doubleSV[addJetIndex[1]]<0.6) {
      h2->Fill(msubt);
      weight = (0.0702331)+(-0.000122783)*(fatjetPRmass[aa]-125)+(0.000001644981)*(fatjetPRmass[aa]-125)*(fatjetPRmass[aa]-125);
      //cout<<weight<<endl;
      h3->Fill(msubt,weight);
    }
    //one anti-tag lead fails and subl pass
    else if(addjet_doubleSV[addJetIndex[0]]<0.6 && addjet_doubleSV[addJetIndex[1]]>0.6) {
      h4->Fill(msubt);
    }

  } //end of the event loop

  hratio->Divide(hpass,hfail); 
  TFile* outfile = new TFile(Form("%s_bkgest.2.4.2.root",name.data()),"recreate");
  h1->Write("one_anti-tag_lead");
  h2->Write("two_anti-tag");
  h3->Write("two_anti-tag_weighted1");
  h4->Write("one_anti-tag_subl");
  hpass->Write("pass");
  hfail->Write("fail");
  hratio->Write("ratio");
  outfile->Write();

}
