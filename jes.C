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
void jes(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());

  int nbin_mjj = 84;
  TH1F* h_Mjjred1 =  new TH1F("","",nbin_mjj,800,5000);
  TH1F* h_Mjjred2 =  new TH1F("","",nbin_mjj,800,5000);
  TH1F* h_Mjjred3 =  new TH1F("","",nbin_mjj,800,5000);


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
    Float_t*  fatjetCorrUncUp = data.GetPtrFloat("FATjetCorrUncUp");
    Float_t*  fatjetCorrUncDown = data.GetPtrFloat("FATjetCorrUncDown");
    vector<bool>    &FATjetPassIDTight = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));    

    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doubleSV = data.GetPtrFloat("ADDjet_DoubleSV");

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    //---------JES Uncertainty----------
    vector<int> fatjet;
    for(int ij=0; ij<nFJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      fatjet.push_back(ij);
    }

    if(fatjet.size()<2)continue;
    int aa = fatjet[0]; //Mjj[0].second;                                                                                                                           
    int ee = fatjet[1]; //Mjj[0].first;                                                                                                                  
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa);
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);
    TLorentzVector Jet1Up(0,0,0,0);
    TLorentzVector Jet2Up(0,0,0,0);
    TLorentzVector Jet1Dw(0,0,0,0);
    TLorentzVector Jet2Dw(0,0,0,0);

    Double_t pt1up = ((*Jet1)*(1+fatjetCorrUncUp[aa])).Pt();
    Double_t pt2up = ((*Jet2)*(1+fatjetCorrUncUp[ee])).Pt();
    Double_t pt1dw = ((*Jet1)*(1-fatjetCorrUncDown[aa])).Pt();
    Double_t pt2dw = ((*Jet2)*(1-fatjetCorrUncDown[ee])).Pt();

    Double_t eta1up = ((*Jet1)*(1+fatjetCorrUncUp[aa])).Eta();
    Double_t eta2up = ((*Jet2)*(1+fatjetCorrUncUp[ee])).Eta();
    Double_t eta1dw = ((*Jet1)*(1-fatjetCorrUncDown[aa])).Eta();
    Double_t eta2dw = ((*Jet2)*(1-fatjetCorrUncDown[ee])).Eta();

    Double_t phi1up = ((*Jet1)*(1+fatjetCorrUncUp[aa])).Phi();
    Double_t phi2up = ((*Jet2)*(1+fatjetCorrUncUp[ee])).Phi();
    Double_t phi1dw = ((*Jet1)*(1-fatjetCorrUncDown[aa])).Phi();
    Double_t phi2dw = ((*Jet2)*(1-fatjetCorrUncDown[ee])).Phi();

    Double_t m1up = ((*Jet1)*(1+fatjetCorrUncUp[aa])).M();
    Double_t m2up = ((*Jet2)*(1+fatjetCorrUncUp[ee])).M();
    Double_t m1dw = ((*Jet1)*(1-fatjetCorrUncDown[aa])).M();
    Double_t m2dw = ((*Jet2)*(1-fatjetCorrUncDown[ee])).M();

    Jet1Up.SetPtEtaPhiM(pt1up,eta1up,phi1up,m1up);
    Jet2Up.SetPtEtaPhiM(pt2up,eta2up,phi2up,m2up);
    Jet1Dw.SetPtEtaPhiM(pt1dw,eta1dw,phi1dw,m1dw);
    Jet2Dw.SetPtEtaPhiM(pt2dw,eta2dw,phi2dw,m2dw);

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
	bool results = trigResult[it];

	
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
    
    if(fatjet.size()<2)continue;*/
    if(!FATjetPassIDTight[aa])continue; 
    if(!FATjetPassIDTight[ee])continue;

    if(Jet1->Pt()<200)continue;
    if(Jet2->Pt()<200)continue;
    if(Jet1Up.Pt()<200)continue;
    if(Jet2Up.Pt()<200)continue; 
    if(Jet1Dw.Pt()<200)continue;
    if(Jet2Dw.Pt()<200)continue;

    if(fabs(Jet1->Eta())>2.4)continue;
    if(fabs(Jet2->Eta())>2.4)continue;
    if(fabs(Jet1Up.Eta())>2.4)continue;
    if(fabs(Jet2Up.Eta())>2.4)continue;
    if(fabs(Jet1Dw.Eta())>2.4)continue;
    if(fabs(Jet2Dw.Eta())>2.4)continue;

    nPass[2]++;
	  
    Double_t dEta = fabs(Jet1->Eta() - Jet2->Eta());
    if(dEta>1.3)continue;
    Double_t dEtaUp = fabs(Jet1Up.Eta() - Jet2Up.Eta());
    if(dEtaUp>1.3)continue;
    Double_t dEtaDw = fabs(Jet1Dw.Eta() - Jet2Dw.Eta());
    if(dEtaDw>1.3)continue;

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
    int addJetIndexUp[2]={-1,-1};
    int addJetIndexDw[2]={-1,-1};
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
      if(Jet1Up.DeltaR(*Jet3)<0.1 && addJetIndexUp[0] < 0) { addJetIndexUp[0]=ad;}
      if(Jet2Up.DeltaR(*Jet3)<0.1 && addJetIndexUp[1] < 0) { addJetIndexUp[1]=ad;}
      if(Jet1Dw.DeltaR(*Jet3)<0.1 && addJetIndexDw[0] < 0) { addJetIndexDw[0]=ad;}
      if(Jet2Dw.DeltaR(*Jet3)<0.1 && addJetIndexDw[1] < 0) { addJetIndexDw[1]=ad;}
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;
    if(addJetIndexUp[0]<0 || addJetIndexUp[1]<0)continue;
    if(addJetIndexDw[0]<0 || addJetIndexDw[1]<0)continue;

    if(addjet_doubleSV[addJetIndex[0]]<0.6)continue;
    if(addjet_doubleSV[addJetIndex[1]]<0.6)continue;
    if(addjet_doubleSV[addJetIndexUp[0]]<0.6)continue;
    if(addjet_doubleSV[addJetIndexUp[1]]<0.6)continue;
    if(addjet_doubleSV[addJetIndexDw[0]]<0.6)continue;
    if(addjet_doubleSV[addJetIndexDw[1]]<0.6)continue;

    nPass[7]++;
 
    Float_t msubtup = (Jet1Up+Jet2Up).M() - (fatjetPRmassL2L3Corr[aa]-125) - (fatjetPRmassL2L3Corr[ee]-125);
    Float_t msubtdw = (Jet1Dw+Jet2Dw).M() - (fatjetPRmassL2L3Corr[aa]-125) - (fatjetPRmassL2L3Corr[ee]-125);
    if(msubtup<800)continue;
    if(msubtdw<800)continue;

    h_Mjjred1->Fill(msubt);
    h_Mjjred2->Fill(msubtup);
    h_Mjjred3->Fill(msubtdw);




  } //end of the event loop

  setNCUStyle();

  double lowbin = h_Mjjred2->FindFirstBinAbove(0,1);
  double highbin = h_Mjjred2->FindLastBinAbove(0,1);



  //TFile* outfile = new TFile("dbtsf.root","update");

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

  h_Mjjred1->Scale((2.6837)/nPass[0]);
  h_Mjjred2->Scale((2.6837)/nPass[0]);
  h_Mjjred3->Scale((2.6837)/nPass[0]);

  bool BulkGrav=(inputFile.find("BulkGrav")!= std::string::npos);
  /*
  TFile* outfile = new TFile("jes.root","update");
  if(BulkGrav) {
    for(int i=0;i<nBM;i++){
      bool bulkmass=(inputFile.find(Form("%s",bulkg_name[i].data()))!= std::string::npos); 
      if(bulkmass) {
	h_Mjjred1->Write(Form("BulkGravM%s_JES",bulkg_name[i].data()));
        h_Mjjred2->Write(Form("BulkGravM%s_JESUp",bulkg_name[i].data()));
        h_Mjjred3->Write(Form("BulkGravM%s_JESDown",bulkg_name[i].data()));
	outfile->Write();
      }
    }
  }
  else {
    for(int i=0;i<nRM;i++){
      bool radionmass=(inputFile.find(Form("%s",radion_name[i].data()))!= std::string::npos);
      if(radionmass) {
	h_Mjjred1->Write(Form("RadionM%s_JES",radion_name[i].data()));
        h_Mjjred2->Write(Form("RadionM%s_JESUp",radion_name[i].data()));
        h_Mjjred3->Write(Form("RadionM%s_JESDown",radion_name[i].data()));
	outfile->Write();
      }
    }
  }
  */

  //---------------------------------------------------------------


  double binwidth = h_Mjjred2->GetBinWidth(lowbin);
  float xlow = (binwidth*lowbin)+(800-(50+100)); //+800 to get to the bin range, -50 to get to the first bin range; -100 for adjustment 
  float xhigh = (binwidth*highbin)+(800+400); //+800 to get to the bin range, +400 for adjustment
  float yhigh = 1.1*(h_Mjjred2->GetMaximum());
  float lowbincont = h_Mjjred2->GetBinContent(lowbin);

  cout<<xlow<<" "<<xhigh<<" "<<lowbin<<" "<<highbin<<" "<<binwidth<<endl;
  
  TCanvas* c2 =  new TCanvas("c2","c2",0,0,600,600);
  c2->cd();
 
  TH1F* hr = c2->DrawFrame(xlow,0,xhigh,yhigh,""); // new TH1F("","",nbin_mjj,xlow,xhigh);
  hr->GetYaxis()->SetRange(0,yhigh);
  hr->GetYaxis()->SetTitle("");
  hr->GetXaxis()->SetTitle("reduced dijet mass");
  hr->GetXaxis()->SetTitleSize(0.04);
  hr->GetYaxis()->SetTitleSize(0.04);
  hr->GetXaxis()->SetLabelSize(0.04);
  hr->GetYaxis()->SetLabelSize(0.04);
  hr->Draw();
  
  h_Mjjred1->SetLineWidth(2);
  h_Mjjred2->SetLineWidth(2);
  h_Mjjred3->SetLineWidth(2);
  h_Mjjred1->SetLineStyle(1);
  h_Mjjred2->SetLineStyle(9);
  h_Mjjred3->SetLineStyle(5);
  h_Mjjred2->Draw("histsame");
  //h_Mjjred2->GetYaxis()->SetTitle("");
  //h_Mjjred2->GetXaxis()->SetTitle("reduced dijet mass");
  //h_Mjjred2->GetXaxis()->SetTitleSize(0.04);
  //h_Mjjred2->GetYaxis()->SetTitleSize(0.04);
  h_Mjjred2->SetLineColor(kOrange+2);
  h_Mjjred1->SetLineColor(kViolet+2);
  h_Mjjred1->Draw("histsame");
  h_Mjjred3->SetLineColor(kGreen+2);
  h_Mjjred3->Draw("histsame");
  TLegend *leg1 = new TLegend(0.70, 0.72, 0.97, 0.88);
  leg1->SetBorderSize(0);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(h_Mjjred2, "up", "l");
  leg1->AddEntry(h_Mjjred1, "central", "l");
  leg1->AddEntry(h_Mjjred3, "down", "l");
  leg1->Draw();


  Double_t central1 = h_Mjjred1->Integral();
  Double_t up1      = h_Mjjred2->Integral();
  Double_t down1    = h_Mjjred3->Integral();

  float sysmjjred = abs(up1-down1)/(2*central1);
  float sys1up = abs(up1-central1)/central1;
  float sys1down = abs(down1-central1)/central1;

  cout<<sysmjjred<<" "<<sys1up<<" "<<sys1down<<" "<<endl;




  
  if(BulkGrav) {
    for(int i=0;i<nBM;i++){
      bool bulkmass=(inputFile.find(Form("%s",bulkg_name[i].data()))!= std::string::npos);
      if(bulkmass) {
	//c2->Print(Form("JES_Mjjred_SysUnc_BulkGrav_%s.pdf",bulkg_name[i].data()));
	//c->Print(Form("DBTSF_pT_BulkGrav_%s.pdf",bulkg_name[i].data()));
        //ofstream fout;
        //fout.open("bulkg_jes_sysunc.dat",ios::out | ios::app);
        //fout<<sys1up<<" "<<sys1down<<endl;
        //fout.close();
      }
    }
  }
  else {
    for(int i=0;i<nRM;i++){
      bool radionmass=(inputFile.find(Form("%s",radion_name[i].data()))!= std::string::npos);
      if(radionmass) {
	//c2->Print(Form("JES_Mjjred_SysUnc_Radion_%s.pdf",radion_name[i].data()));
	//c->Print(Form("DBTSF_pT_Radion_%s.pdf",radion_name[i].data()));
        //ofstream fout;
        //fout.open("radion_jes_sysunc.dat",ios::out | ios::app);
        //fout<<sys1up<<" "<<sys1down<<endl;
        //fout.close();
      }
    }
  }
  


}
