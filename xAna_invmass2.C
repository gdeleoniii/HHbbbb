// example code to run Bulk Graviton->ZZ->ZlepZhad selections on electron-channel

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1D.h>
#include <TFile.h>
#include <TF1.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>

using namespace std;
void xAna_hh(std::string inputFile,char name) {

  //get TTree from file ...
  TreeReader data(inputFile.data());

  Long64_t nTotal=0;
  Long64_t nPass[20]={0};

  TH1F* h_SD=new TH1F("h_SD","",100,0,200);
  TH1F* h_FatMass=new TH1F("h_FatMass","",100,0,5000);
 
 
  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    nTotal++;

    //0. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
    nPass[0]++;

    int nFATJet         = data.GetInt("FATnJet");
    const int nJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
    Float_t*  fatjetSDmass = data.GetPtrFloat("FATjetSDmass");
    Int_t*   nSubSoftDropJet = data.GetPtrInt("FATnSubSDJet");
    vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);
    vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJet);
    vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJet);
    vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJet);
    vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDCE", nFATJet);
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));

    vector<int> fatty;
    int nFatBTag=0;
    for(int ij=0; ij<nJets; ij++)
      {
    	
     	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
    	if(thisJet->Pt()<300)continue;
    	if(fabs(thisJet->Eta())>2.5)continue;
    	if(fatjetSDmass[ij]<95 || fatjetSDmass[ij]>145)continue;
	if(!passFatJetLooseID[ij])continue;

    	if(fatjetCISVV2[ij] < 0.605)continue;
    	if(fatjetCISVV2[ij] > 1)continue;

    	nFatBTag++;
	fatty.push_back(ij);
      }
    
    if(nFatBTag>=2)nPass[1]++;   
    
    double MaxPt= 0;
    int leadjet1 = -99;
    vector<int> leader1;
    vector<int> leader2;
    for(unsigned int x=0; x<fatty.size(); x++) { //for the 1st leading jet
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(x);
      if(thisJet->Pt()>MaxPt) leadjet1 = fatty[x]; 
      leader1.push_back(x);
    }

    for(unsigned int y=0; y<fatty.size(); y++) { //for the 2nd leading jet
      if(fatty[y]==leadjet1)continue;
      leader2.push_back(y);
    }


    int nSubBTag[2]={0}; // check only the leading two fat jets 
    int nGoodFatJet=0;
    for(int ij=0; ij<nJets; ij++)
      {
    	
	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
	if(thisJet->Pt()<30)continue;
	if(fabs(thisJet->Eta())>2.5)continue;
	if(fatjetSDmass[ij]<95 || fatjetSDmass[ij]>145)continue;
	if(!passFatJetLooseID[ij])continue;

      	for(int is=0; is < nSubSoftDropJet[ij]; is++){
	  if(subjetSDCSV[ij][is] < 0.605)continue;
	  if(nGoodFatJet<2)
	  nSubBTag[nGoodFatJet]++;
	  
      	}

	nGoodFatJet++;
      }
  
    
    // if each fat jet has at least one subjet btag
    if(nSubBTag[0]>0 && nSubBTag[1]>0)nPass[2]++;

    // if one of the fat jets has at least two subjet btags
    if((nSubBTag[0]>1 && nSubBTag[1]>0) || 
       (nSubBTag[0]>0 && nSubBTag[1]>1))nPass[3]++;

    // if both fat jets have at least two subjet btags
    if(nSubBTag[0]>1 && nSubBTag[1]>1) nPass[4]++;

    vector<double> fatjetM;
    for(unsigned int g=0; g<fatty.size(); g++) {
      int gg = leader1[g];
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(gg);

      for(unsigned int h=0; h<g; h++) {
	int hh = leader2[h];
 	TLorentzVector* thatJet = (TLorentzVector*)fatjetP4->At(hh);

	Float_t mff = (*thisJet+*thatJet).M();
	fatjetM.push_back(mff);
      }
    }

    for(unsigned int g=0; g<fatjetM.size(); g++) {
      h_FatMass->Fill(fatjetM[g]);
    }

  } // end of loop over entries

  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;

  h_FatMass->Draw();
  
  TFile* outfile0 = new TFile(Form("fatmass_%d.root",name),"recreate");
  h_FatMass->Write("invmass");
  outfile0->Write();


}
