#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1.h>
#include <TLine.h>
#include <TFrame.h>
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
void pdfscale(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());

  int nbin_mjj = 84;
  TH1F* h_Mjjred[9];
  for(int i=0;i<9;i++) {
    h_Mjjred[i] =  new TH1F("","",50,0,1);

  }
  Float_t nPass[20]={0};
  
  TH1F* h_pdf = new TH1F("","",50,0,1);

  Float_t scalea[9],scaleb[9],pdfa[100],pdfb[100];

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
    Float_t*  fatjetCorrUncUp = data.GetPtrFloat("FATjetCorrUncUp");
    Float_t*  fatjetCorrUncDown = data.GetPtrFloat("FATjetCorrUncDown");
    Float_t*  pdfscaleSysWeight = data.GetPtrFloat("pdfscaleSysWeights"); 
    vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));    

    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doubleSV = data.GetPtrFloat("ADDjet_DoubleSV");

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();


    for(int i=0;i<109;i++) {
      if(i<9) {
	scaleb[i] += pdfscaleSysWeight[i];
	if(pdfscaleSysWeight[i]<=0){
	  cout<<" scale before = "<<pdfscaleSysWeight[i]<<endl;
	}
      }
      else if(i>= 9) {
	pdfb[i-9] += pdfscaleSysWeight[i];
      }
    }

    if(nFJets<2)continue;
    int aa = 0; //Mjj[0].second;                                                                                                                           
    int ee = 1; //Mjj[0].first;                                                                                                                  
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa);
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
	bool results = trigResult[it];

	// std::cout << thisTrig << " : " << results << std::endl;
	
	if( (thisTrig.find("HLT_PFHT800")!= std::string::npos && results==1)
	    )
	  {
	    if(thisTrig != "HLT_PFHT800_v2") {cout<<thisTrig<<endl;}
	    passTrigger=true;
	    break;
	  }


      }


    if(!passTrigger)continue;

    nPass[1]++;   

    //3. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;

    if(!FATjetPassIDTight[aa])continue; 
    if(!FATjetPassIDTight[ee])continue;

    if(Jet1->Pt()<200)continue;
    if(Jet2->Pt()<200)continue;

    if(fabs(Jet1->Eta())>2.4)continue;
    if(fabs(Jet2->Eta())>2.4)continue;

    nPass[2]++;

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

    for(int i=0;i<109;i++) {
      if(i<9) {
        scalea[i] += pdfscaleSysWeight[i];
	//cout<<" scale after = "<<pdfscaleSysWeight[i]<<endl;
      }
      if(i>= 9) {
        pdfa[i-9] += pdfscaleSysWeight[i];
      }
    }
 
  } //end of the event loop

  //setNCUStyle();
  std::string bulkg_name[]={"1000","1200","1400","1600","1800","2000","2500","3000","4000","4500"};
  std::string radion_name[]={"1000","1200","1400","1600","1800","2000","2500","3000","3500","4500"};
  //std::string bulkg_name[]={"1000","1800","2000","2500","3000","4500"};
  //std::string radion_name[]={"1000","1800","2000","2500","3000","3500","4500"};
  float bmass[10]={1000,1200,1400,1600,1800,2000,2500,3000,4000,4500};
  float rmass[10]={1000,1200,1400,1600,1800,2000,2500,3000,3500,4500};

  //float bmass[4]={1000,1800,2000,2500};
  //float rmass[7]={1000,1800,2000,2500,3000,3500,4500};
  int nBM = sizeof(bmass)/sizeof(bmass[0]);
  int nRM = sizeof(rmass)/sizeof(rmass[0]);


  //--------------------xsec------------------------------------
  //float bxsec[6]={5.6657071801246945,0.18315617607742797,0.09062263600103764,0.01868308774736877,4.412202866735033e-3,8.933087336020322e-5};
  //float rxsec[7]={261.8915924170746,23.912714207366757,14.293231312611458,4.301706042824362,1.3843144368847916,0.46996032032601366,0.0534104307518241};

  //float bBR[6]={0.09869019214731996,0.09957353475983793,0.09961340921364503,0.09966492755238268,0.09968811468551785,0.09971268196476021};
  //float rBR = 0.236958;

  
  
  bool BulkGrav=(inputFile.find("BulkGrav")!= std::string::npos);
  //-------------------PDF------------------------------------------
  float pdfe[100];
  for(int i=0;i<100;i++) {
    pdfe[i] = pdfa[i]/pdfb[i];
    h_pdf->Fill(pdfe[i]);
    cout<<"pdf = "<<pdfe[i]<<" "<<pdfa[i]<<" "<<pdfb[i]<<endl;
  }
  
  
  TCanvas* c =  new TCanvas("c","c",0,0,600,600);
  c->cd();
  h_pdf->Draw();
  //cout<<h_pdf->GetMean()<<" "<<h_pdf->GetRMS()<<" "<<h_pdf->GetRMS()/h_pdf->GetMean()<<endl;
  //------------------------------------------------------------------
  

  TCanvas* c2 =  new TCanvas("c2","c2",0,0,600,600);
  h_Mjjred[0]->Draw();
  double integ[9],sysuncScale[9],mean[9],rms[9], scalee[9];
  for(int i= 0;i<9;i++) {
    //cout<<scalee[i]<<" "<<scalea[i]<<" "<<scaleb[i]<<endl; 
    scalee[i] = scalea[i]/scaleb[i];
    cout<<"scale = "<<scalee[i]<<" "<<scalea[i]<<" "<<scaleb[i]<<endl;
    h_Mjjred[i]->Fill(scalee[i]);
  

    h_Mjjred[i]->SetLineWidth(2);
    h_Mjjred[i]->SetLineColor(i+1);
    h_Mjjred[i]->SetLineStyle(i+1);
    h_Mjjred[i]->Draw("histsame"); 
    integ[i] = h_Mjjred[i]->Integral();
    mean[i] = h_Mjjred[i]->GetMean();
    rms[i] = h_Mjjred[i]->GetRMS();
  }
  h_Mjjred[0]->SetLineWidth(4);


  TLegend *leg1 = new TLegend(0.70, 0.66, 0.97, 0.88);                                                                                                              
  leg1->SetBorderSize(0);                                                                                                                                           
  leg1->SetFillColor(0);                                                                                                                                            
  leg1->SetFillStyle(0);                                                                                                                                            
  leg1->SetTextSize(0.035);
  for(int i= 0;i<9;i++) {
    leg1->AddEntry(h_Mjjred[i], Form("%i",i), "l");
  }
  leg1->Draw();



  for(int i= 0;i<9;i++) {
    sysuncScale[i] = fabs(integ[0] - integ[i])/integ[0];
    cout<<"scale uncertainty = "<<sysuncScale[i]<<endl;
  }
  for(int i= 0;i<9;i++) {
   
    cout<<"mean = "<<mean[i]<<endl;
  }
  for(int i= 0;i<9;i++) {
     cout<<"rms = "<<rms[i]<<endl;
  }  
 


}
