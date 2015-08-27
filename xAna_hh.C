// example code to run Bulk Graviton->ZZ->ZlepZhad selections on electron-channel

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>

using namespace std;
void xAna_hh(std::string inputFile){

  //get TTree from file ...
  TreeReader data(inputFile.data());

  Long64_t nTotal=0;
  Long64_t nPass[20]={0};

  TH1F* h_SD=new TH1F("h_SD","",100,0,200);
  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    nTotal++;

    //0. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
    nPass[0]++; //this part put the good vertices inside the nPass[0]

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
    
    int nFatBTag=0;
    for(int ij=0; ij<nJets; ij++)
      {
    	
     	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
    	if(thisJet->Pt()<30)continue;
    	if(fabs(thisJet->Eta())>2.5)continue; 
    	if(fatjetSDmass[ij]<95 || fatjetSDmass[ij]>145)continue; 
	//cuts for leading jets

    	if(fatjetCISVV2[ij] < 0.605)continue;
    	if(fatjetCISVV2[ij] > 1)continue; //this is the b-tagging for the fat jets

    	nFatBTag++;

      }
    
    if(nFatBTag>=2)nPass[1]++; //2 or greater fat jets(with cuts and b-tagging) that has b-tagging will be put inside the nPass[1] 


    int nSubBTag[2]={0}; // check only the leading two fat jets 
    int nGoodFatJet=0;
    for(int ij=0; ij<nJets; ij++)
      {
    	
	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
	if(thisJet->Pt()<30)continue; 
	if(fabs(thisJet->Eta())>2.5)continue; 
	if(fatjetSDmass[ij]<95 || fatjetSDmass[ij]>145)continue; 
	//cuts for leading jets

      	for(int is=0; is < nSubSoftDropJet[ij]; is++){
	  if(subjetSDCSV[ij][is] < 0.605)continue; //CSV distinguish between jets from b quarks and those form charm or light quarks and gluon
	  if(nGoodFatJet<2)
	    nSubBTag[nGoodFatJet]++;  
      	} //this loop is for the b-tagging on the subjet

	nGoodFatJet++;

      }
  
    // if each fat jet has at least one subjet btag
    if(nSubBTag[0]>0 && nSubBTag[1]>0)nPass[2]++;

    // if one of the fat jets has at least two subjet btags
    if((nSubBTag[0]>1 && nSubBTag[1]>0) || 
       (nSubBTag[0]>0 && nSubBTag[1]>1))nPass[3]++;

    // if both fat jets have at least two subjet btags
    if(nSubBTag[0]>1 && nSubBTag[1]>1) nPass[4]++;
        

  } // end of loop over entries

  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;
    
}
/*nPass[0] = for the good vertices
nPass[1] = leading jets with b-tagging (2 or more)
nPass[2] = both fat jets that has at least one subjet b-tag
nPass[3] = one of the fat jets that have at least two subjet b-tag
nPass[4] = both fat jets that have at least two subjet b-tag */
