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
#include <TRandom.h>
#include <THnSparse.h>
#include <TStyle.h>
#include <TStyle.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "setNCUStyle.C"

double ApplyJERp4(double eta, int jerShift) {
  //76X 2015 DATA/MC SFs

  eta = abs(eta) ;
  if (eta >= 4.7) eta = 4.699 ;
  double jerscale(0) ; 

  if (eta >= 3.2 && eta < 5.0) {
    if (jerShift == 0) jerscale = 1.216  ;  
    else if (jerShift == 1) jerscale = 1.216 + 0.050 ;  
    else if (jerShift == -1) jerscale = 1.216 - 0.050 ;  
  }
  else if (eta >= 3.0 && eta < 3.2) {
    if (jerShift == 0) jerscale = 1.384  ;  
    else if (jerShift == 1) jerscale = 1.384 + 0.033;  
    else if (jerShift == -1) jerscale = 1.384 - 0.033;  
  }
  else if (eta >= 2.8 && eta < 3.0) {
    if (jerShift == 0) jerscale =  1.564 ;  
    else if (jerShift == 1) jerscale =  1.564 + 0.321;  
    else if (jerShift == -1) jerscale =  1.564 - 0.321;  
  }
  else if (eta >= 2.5 && eta < 2.8) {
    if (jerShift == 0) jerscale = 1.209  ;  
    else if (jerShift == 1) jerscale = 1.209 + 0.059;  
    else if (jerShift == -1) jerscale = 1.209 - 0.059;  
  }
  else if (eta >= 2.3 && eta < 2.5) {
    if (jerShift == 0) jerscale = 1.161 ;  
    else if (jerShift == 1) jerscale = 1.161 + 0.060;  
    else if (jerShift == -1) jerscale = 1.161 - 0.060;  
  }
  else if (eta >= 2.1 && eta < 2.3) {
    if (jerShift == 0) jerscale = 1.160 ;  
    else if (jerShift == 1) jerscale = 1.160 + 0.048;  
    else if (jerShift == -1) jerscale = 1.160 - 0.048 ;  
  }
  else if (eta >= 1.9 && eta < 2.1) {
    if (jerShift == 0) jerscale = 1.162  ;  
    else if (jerShift == 1) jerscale = 1.162 + 0.044;  
    else if (jerShift == -1) jerscale = 1.162 - 0.044;  
  }
  else if (eta >= 1.7 && eta < 1.9) {
    if (jerShift == 0) jerscale = 1.100  ;  
    else if (jerShift == 1) jerscale = 1.100 + 0.033;  
    else if (jerShift == -1) jerscale = 1.100 - 0.033 ;  
  }
  else if (eta >= 1.3 && eta < 1.7) {
    if (jerShift == 0) jerscale = 1.118  ;  
    else if (jerShift == 1) jerscale = 1.118 + 0.014;  
    else if (jerShift == -1) jerscale = 1.118 - 0.014 ;  
  }
  else if (eta >= 1.1 && eta < 1.3) {
    if (jerShift == 0) jerscale = 1.103  ;  
    else if (jerShift == 1) jerscale = 1.103 + 0.033;  
    else if (jerShift == -1) jerscale = 1.103 - 0.033;  
  }
  else if (eta >= 0.8 && eta < 1.1) {
    if (jerShift == 0) jerscale = 1.097  ;  
    else if (jerShift == 1) jerscale = 1.097 + 0.017;  
    else if (jerShift == -1) jerscale = 1.097 - 0.017;  
  }
  else if (eta >= 0.5 && eta < 0.8) {
    if (jerShift == 0) jerscale = 1.120  ;  
    else if (jerShift == 1) jerscale = 1.120 + 0.028;  
    else if (jerShift == -1) jerscale = 1.120 - 0.028;  
  }
  else if (eta >= 0.0 && eta < 0.5) {
    if (jerShift == 0) jerscale = 1.095  ;  
    else if (jerShift == 1) jerscale = 1.095 + 0.018;  
    else if (jerShift == -1) jerscale = 1.095 - 0.018;  
  }

  return jerscale ; 
}

using namespace std;
void jer(std::string inputFile, int mode) {
  if(mode == 0) {
    cout<<"jer is central"<<endl;
  }
  else if(mode == 1) {
    cout<<"jer is up"<<endl;
  }
  else if(mode == -1) {
    cout<<"jer is down"<<endl;
  }

  //get TTree from file ...
  TreeReader data(inputFile.data());

  int nbin_mjj = 84;
  TH1F* h_Mjjred1 =  new TH1F("","",nbin_mjj,800,5000);
  h_Mjjred1->Sumw2();

  Float_t nPass[20]={0};

  //------opening pt resolution text file--------------------
  int reso = 104;
  double eta1[reso],eta2[reso],rho1[reso],rho2[reso],num[reso],pt1[reso],pt2[reso],c0[reso],c1[reso],c2[reso],c3[reso];
  ifstream fin;
  fin.open("jer_ptresolution.txt");
  for(int i=0;i<reso;i++) {
    fin >> eta1[i] >> eta2[i] >> rho1[i] >> rho2[i] >> num[i] >> pt1[i] >> pt2[i] >> c0[i] >> c1[i] >> c2[i] >> c3[i];
  }
  fin.close();
  //---------------------------------------------------------

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
    Float_t fatjetRho = data.GetFloat("FATjetRho");

    TClonesArray* fatgenjetP4 = (TClonesArray*) data.GetPtrTObject("FATgenjetP4");
    Int_t nGenPar        = data.GetInt("nGenPar");
    Int_t* genParId      = data.GetPtrInt("genParId");
    Int_t* genParSt      = data.GetPtrInt("genParSt");
    Int_t* genMomParId   = data.GetPtrInt("genMomParId");
    Int_t* genDa1      = data.GetPtrInt("genDa1");
    Int_t* genDa2      = data.GetPtrInt("genDa2");


    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doubleSV = data.GetPtrFloat("ADDjet_DoubleSV");

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    //---------JER Uncertainty----------
   
    if(nFJets < 2)continue;

    int aa = 0; //Mjj[0].second;                                                                                                                           
    int ee = 1; //Mjj[0].first;                                                                                                                  
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa);
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);

    /*
    double ptsmearGlobal[3][2];
    for(int ij=0; ij<2; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      TLorentzVector* thisGenJet = (TLorentzVector*)fatgenjetP4->At(ij);

      if(thisGenJet->E() < 0 || abs(thisGenJet->E()) > 90000)continue;

      //cout<< ij <<" "<< thisGenJet->E() <<endl;

      double pt_gen = thisGenJet->Pt();
      double pt_reco   = thisJet->Pt();
      double eta_reco = thisJet->Eta();
      double ptreso;
      for(int g=0;g<reso;g++) {
	if((eta_reco>eta1[g] && eta_reco<eta2[g]) && (fatjetRho>rho1[g] && fatjetRho<rho2[g]) && (pt_reco>pt1[g] && pt_reco<pt2[g])) {
	  ptreso = sqrt(c0[g]*abs(c0[g])/(pt_reco*pt_reco)+c1[g]*c1[g]*pow(pt_reco,c3[g])+c2[g]*c2[g]);
	  //cout<<nPass[0]<<" "<<g<<endl;
	}
      }

      double jerscalep4[3];
      double ptsmear[3];
      jerscalep4[0]= ApplyJERp4(eta_reco,0);
      jerscalep4[1]= ApplyJERp4(eta_reco,1);
      jerscalep4[2]= ApplyJERp4(eta_reco,-1);
      for(int k=0;k<3;k++) { 
	if((thisJet->DeltaR(*thisGenJet)>0.4) && (abs(pt_reco - pt_gen) < 3*ptreso*pt_reco)) {
	    ptsmear[k] = 1 + (jerscalep4[k] -1)*((pt_reco -pt_gen)/pt_reco);
	  }
	    else {
	      TRandom* rand = new TRandom();
	      //double maxi = std::max(jerscalep4[k]*jerscalep4[k] -1,0.0);
	      ptsmear[k] = 1 + (rand->Gaus(0,ptreso))*sqrt(std::max(jerscalep4[k]*jerscalep4[k] -1,0.0));
	      delete rand;
	    }
	    ptsmearGlobal[k][ij]=ptsmear[k];
	    if(ptsmearGlobal[k][ij] < 0) {
	      cout<<ptsmearGlobal[k][ij]<<endl;
	      ptsmearGlobal[k][ij] = 0;
	    }
      }
    }


    if(mode ==1) {
      *Jet1 *= ptsmearGlobal[1][0];
      *Jet2 *= ptsmearGlobal[1][1];
    }
    else if(mode == 0) {
      *Jet1 *= ptsmearGlobal[0][0];
      *Jet2 *= ptsmearGlobal[0][1];
    }
    else if(mode == -1) {
      *Jet1 *= ptsmearGlobal[2][0];
      *Jet2 *= ptsmearGlobal[2][1];
    }
*/


    double ptsmearGlobal[3][2];
    bool leadingMatchleading=0;
    for(int ij=0; ij<2; ij++)
      {
	if( !FATjetPassIDTight[ij] )continue;
	
	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
	TLorentzVector* thisGenJet = (TLorentzVector*)fatgenjetP4->At(0);

	if(thisGenJet->E() < 0 || abs(thisGenJet->E()) > 90000)continue;

	cout<< ij <<" "<< thisGenJet->E() <<endl;

	if(ij==0){
	  leadingMatchleading=1;
	  if(thisJet->DeltaR(*thisGenJet)>0.4){
	    thisGenJet = (TLorentzVector*)fatgenjetP4->At(1);
	    if(ij==0)leadingMatchleading=0;
	  }
	}
	if(ij==1 && leadingMatchleading)thisGenJet = (TLorentzVector*)fatgenjetP4->At(1);
	double pt_gen = thisGenJet->Pt();
	double pt_reco   = thisJet->Pt() ;
	double eta_reco = thisJet->Eta() ; 
	double jerscalep4[3];
	double ptsmear[3];
	jerscalep4[0]= ApplyJERp4(eta_reco, 0) ; 
	jerscalep4[1]= ApplyJERp4(eta_reco, 1) ; 
	jerscalep4[2]= ApplyJERp4(eta_reco, -1) ;
	for(int k=0;k<3;k++){
	  if (pt_gen > 0.) ptsmear[k] = std::max( 0.0, pt_gen + jerscalep4[k]*(pt_reco - pt_gen) )/pt_reco ; 
	  else if (jerscalep4[k] > 1. && pt_reco > 0.) {
	    TRandom* rand = new TRandom();
	    ptsmear[k] = rand->Gaus(pt_reco, sqrt(jerscalep4[k]*jerscalep4[k] - 1)*0.2)/pt_reco ; //// Assuming 20% JER
	    delete rand; 
	  }
	  ptsmearGlobal[k][ij]=ptsmear[k] ;
	  if(ptsmearGlobal[k][ij] < 0) {
	    cout<<ptsmearGlobal[k][ij]<<endl;
	    ptsmearGlobal[k][ij] = 0;
            }
	  
	}
      }

    if(mode ==1) {
      *Jet1 *= ptsmearGlobal[1][0]/ptsmearGlobal[0][0];                                                                                                             
      *Jet2 *= ptsmearGlobal[1][1]/ptsmearGlobal[0][1];      
    }
    else if(mode == 0) {
      *Jet1 *= ptsmearGlobal[0][0]/ptsmearGlobal[0][0];
      *Jet2 *= ptsmearGlobal[0][1]/ptsmearGlobal[0][1];
    }
    else if(mode == -1) {
      *Jet1 *= ptsmearGlobal[2][0]/ptsmearGlobal[0][0];
      *Jet2 *= ptsmearGlobal[2][1]/ptsmearGlobal[0][1];
    }



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


    //cout<<msubt<<endl;
    h_Mjjred1->Fill(msubt);


  } //end of the event loop

  setNCUStyle();

  double lowbin = h_Mjjred1->FindFirstBinAbove(0,1);
  double highbin = h_Mjjred1->FindLastBinAbove(0,1);

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

  bool BulkGrav=(inputFile.find("BulkGrav")!= std::string::npos);
  
  TFile* outfile = new TFile("jer_chingwei.root","update");
  if(BulkGrav) {
    for(int i=0;i<nBM;i++){
      bool bulkmass=(inputFile.find(Form("%s",bulkg_name[i].data()))!= std::string::npos); 
      if(bulkmass) {
	if(mode == 0){
	  h_Mjjred1->Write(Form("BulkGravM%s_JER",bulkg_name[i].data()));
	}
	if(mode == 1) {
	    h_Mjjred1->Write(Form("BulkGravM%s_JERUp",bulkg_name[i].data()));
	  }
	if(mode == -1) {
        h_Mjjred1->Write(Form("BulkGravM%s_JERDown",bulkg_name[i].data()));
	}
	outfile->Write();
      }
    }
  }
  else {
    for(int i=0;i<nRM;i++){
      bool radionmass=(inputFile.find(Form("%s",radion_name[i].data()))!= std::string::npos);
      if(radionmass) {
	if(mode == 0) {
	  h_Mjjred1->Write(Form("RadionM%s_JER",radion_name[i].data()));
	}
	if(mode == 1) {
	  h_Mjjred1->Write(Form("RadionM%s_JERUp",radion_name[i].data()));
	}
	if(mode == -1) {
	  h_Mjjred1->Write(Form("RadionM%s_JERDown",radion_name[i].data()));
	}
	outfile->Write();
      }
    }
  }
  
  
  //---------------------------------------------------------------

  /*
  double binwidth = h_Mjjred1->GetBinWidth(lowbin);
  float xlow = (binwidth*lowbin)+(800-(50+100)); //+800 to get to the bin range, -50 to get to the first bin range; -100 for adjustment 
  float xhigh = (binwidth*highbin)+(800+400); //+800 to get to the bin range, +400 for adjustment
  float yhigh = 1.1*(h_Mjjred1->GetMaximum());
  float lowbincont = h_Mjjred1->GetBinContent(lowbin);

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
  //h_Mjjred2->Draw("histsame");
  //h_Mjjred2->GetYaxis()->SetTitle("");
  //h_Mjjred2->GetXaxis()->SetTitle("reduced dijet mass");
  //h_Mjjred2->GetXaxis()->SetTitleSize(0.04);
  //h_Mjjred2->GetYaxis()->SetTitleSize(0.04);
  h_Mjjred2->SetLineColor(kOrange+2);
  h_Mjjred1->SetLineColor(kViolet+2);
  h_Mjjred1->Draw("histsame");
  h_Mjjred3->SetLineColor(kGreen+2);
  // h_Mjjred3->Draw("histsame");
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
  
  */

}
