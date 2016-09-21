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

  std::vector<string> infiles;

  readSample(inputFile, infiles);
    
  TreeReader data(infiles);

  float xbins1[] = {-75,-45,-20,10,40,75};
  const int nbins1=sizeof(xbins1)/sizeof(xbins1[0])-1;
  TH1F* h1 = new TH1F("","",nbins1,xbins1);
  TH1F* h2 = new TH1F("","",nbins1,xbins1);
  TH1F* h3 = new TH1F("","",nbins1,xbins1);
  TH1F* h4 = new TH1F("","",34,40,210);
  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();

  Long64_t counter[20] = {0};

  total += data.GetEntriesFast();
  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);

    int nFATJet         = data.GetInt("nFatjetAK08ungroomed");
    const int nFJets=nFATJet;
    Float_t*  FatjetMass = data.GetPtrFloat("FatjetAK08ungroomed_mass");
    Float_t*  FatjetPRmass = data.GetPtrFloat("FatjetAK08ungroomed_mpruned");
    Float_t*  FatjetPt = data.GetPtrFloat("FatjetAK08ungroomed_pt");
    Float_t*  FatjetEta = data.GetPtrFloat("FatjetAK08ungroomed_eta");
    Float_t*  FatjetIDTight = data.GetPtrFloat("FatjetAK08ungroomed_id_Tight");
    Float_t*  FatjetPRmassCorr = data.GetPtrFloat("FatjetAK08ungroomed_mprunedcorr");
    Float_t*  FatjetTau1 = data.GetPtrFloat("FatjetAK08ungroomed_tau1");
    Float_t*  FatjetTau2 = data.GetPtrFloat("FatjetAK08ungroomed_tau2");
    Float_t*  FatjetDBT = data.GetPtrFloat("FatjetAK08ungroomed_bbtag");
    //Int_t* HLT_PFHT800 = data.GetPtrInt("HLT_BIT_HLT_PFHT800_v");


    counter[0]++;

    //if(!HLT_PFHT800)continue;

    counter[1]++;

    vector<int> fatjet;
    for(int ij=0;ij<nFJets; ij++) {
      if(FatjetPt[ij]<200)continue;
      if(fabs(FatjetEta[ij])>2.4)continue;
      if(!FatjetIDTight[ij])continue;
      if(FatjetPRmassCorr[ij]<60 || FatjetPRmassCorr[ij]>200)continue;
  
      Double_t tau21 = (FatjetTau2[ij]/FatjetTau1[ij]);
      if(tau21>0.6)continue;

      fatjet.push_back(ij);
    }
  
    if(fatjet.size()<2)continue;
    counter[2]++;

    int lead = fatjet[0];
    int subl = fatjet[1];

    Float_t DelEta = fabs(FatjetEta[lead]-FatjetEta[subl]);
    if(DelEta>1.3)continue;
    Float_t mjj = FatjetPRmass[lead]+FatjetPRmass[subl];
    Float_t msubt = mjj-(FatjetPRmass[lead]-125)-(FatjetPRmass[subl]-125);
    if(msubt<800)continue;

    counter[3]++;


    h4->Fill(FatjetPRmass[lead]);

    if(FatjetPRmass[lead]>105 && FatjetPRmass[lead]<135)continue;
    if(FatjetDBT[lead]>0.6)h1->Fill(FatjetPRmass[lead]-125);
    else if(FatjetDBT[lead]<0.6)h2->Fill(FatjetPRmass[lead]-125);

  } //end of the event loop

  cout<<"entries="<<total<<endl;
  cout<<counter[0]<< " "<<counter[1]<<" "<<counter[2]<<" "<<counter[3]<<endl;

  h3->Divide(h1,h2);
  
  TFile* outfile = new TFile(Form("%s_passfailratio.3.4.root",name.data()),"recreate");
  h1->Write("pass");
  h2->Write("fail");
  h3->Write("pfratio");
  h4->Write("prmass");
  outfile->Write();

}
