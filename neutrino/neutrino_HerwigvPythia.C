// example code to run Bulk Graviton->ZZ->ZlepZhad selections on electron-channel

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TPad.h>
#include <TStyle.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "setNCUStyle.C"
//#include "twopads.C"


using namespace std;
void neutrino(std::string inputFile,std::string name) {
  bool isHerwigpp=(inputFile.find("herwigpp")!= std::string::npos);
  cout << endl;
  if(isHerwigpp)cout << "Note!! This is a herwigpp MC sample " << endl;
  else
    cout << "This is a pythia8 MC sample" << endl;
  cout << endl;

  //get TTree from file ...
  TreeReader data(inputFile.data());
  //float xbins1[]={300,325,350,375,400,425,450,475,500,700};
  //float xbins1[]={300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1400};
  //float xbins1[]={300,375,450,525,600,675,750,825,900,975,1050,1125,1200,1275,1350,1425,1500,2000};
  float xbins1[]={300,400,500,600,700,800,900,1000,1200,1300,1400,1500,1600,1700,1800,1900,2000,2750};
  
  const int nbins1=sizeof(xbins1)/sizeof(xbins1[0])-1;

  float xbins2[]={40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,130,140,150,160,170,180,190,200};
  const int nbins2=sizeof(xbins2)/sizeof(xbins2[0])-1;

  TH2F *h_leadJet_neuP4 = new TH2F("","",75,0,3600,45,0,900);
  TProfile* p1 = new TProfile("","",nbins1,xbins1);
  TProfile* p2 = new TProfile("","",nbins1,xbins1);
  TH2F *h_leadJet_nNeu  = new TH2F("","",36,0,3600,17,0,17);
  TH2F *h_sublJet_neuP4 = new TH2F("","",75,0,3600,45,0,900);
  TH2F *h_sublJet_nNeu  = new TH2F("","",36,0,3600,17,0,17);
  TH2F *h_leadPR_neuP4  = new TH2F("","",50,40,200,45,0,900);//
  TProfile* p3 = new TProfile("","",nbins2,xbins2);
  TProfile* p4 = new TProfile("","",nbins2,xbins2);
  TH2F *h_leadPR_nNeu   = new TH2F("","",32,40,200,17,0,17);//
  TH2F *h_sublPR_neuP4  = new TH2F("","",50,40,200,45,0,900);//
  TH2F *h_sublPR_nNeu   = new TH2F("","",32,40,200,17,0,17);//

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    
    int nFATJet         = data.GetInt("FATnJet");
    const int nJets=nFATJet;
    TClonesArray* genParP4    = (TClonesArray*) data.GetPtrTObject("genParP4");
    TClonesArray* fatjetP4    = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCISVV2    = data.GetPtrFloat("FATjetCISVV2");
    Float_t*  fatjetPRmass    = data.GetPtrFloat("FATjetPRmass");
    Int_t     nGenPar         = data.GetInt("nGenPar");
    Int_t*    genParId        = data.GetPtrInt("genParId");
    Int_t*    genParSt        = data.GetPtrInt("genParSt");
    Int_t*    genMomParId     = data.GetPtrInt("genMomParId");
    Int_t*    genGMomParId     = data.GetPtrInt("genGMomParId");
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));

    vector<int> fatjet;
    for(int ij=0; ij<nJets; ij++)
      {

        TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
        if(thisJet->Pt()<300)continue;
        if(fabs(thisJet->Eta())>2.4)continue;
	if(!passFatJetLooseID[ij])continue;
	if(fatjetPRmass[ij]<40)continue;
        if(fatjetCISVV2[ij] < 0.605)continue;
        if(fatjetCISVV2[ij] > 1)continue;

        fatjet.push_back(ij);
      }
    

    if(fatjet.size()<2)continue;

    Int_t lead = fatjet[0];
    Int_t subl = fatjet[1];

    TLorentzVector *leadjet = (TLorentzVector*)fatjetP4->At(lead);
    TLorentzVector *subljet = (TLorentzVector*)fatjetP4->At(subl);
    
    int n_neutrinos_lead=0;
    TLorentzVector neutrinos_lead_p4;
    
    int n_neutrinos_subl=0;
    TLorentzVector neutrinos_subl_p4;    
    
    for(int ig=0; ig < nGenPar; ig++){
      
      int pid =abs(genParId[ig]);
      int mom_pid = abs(genMomParId[ig]);
      int grandmom_pid = abs(genGMomParId[ig]);
      if(genParSt[ig]!=1)continue;

      if(pid!=12 && pid!=14 && pid!=16)continue;
      

      // look for daughters of charm and b hadrons
      // bhadrons can decay to l nu + charm hadrons 
      // charm hadrons could further decay to l nu + other hadrons
      
      if(isHerwigpp) {
        if ( (mom_pid < 400 || mom_pid>600) &&
             (grandmom_pid<400 || grandmom_pid>600) )continue;
      }
      else {
	if(mom_pid < 400 || mom_pid > 600)continue;
      }
      
      TLorentzVector* thisGen = (TLorentzVector*)genParP4->At(ig);
      
      if(thisGen->DeltaR(*leadjet)<0.8) { 
	n_neutrinos_lead++;
	neutrinos_lead_p4 += *thisGen; 
      }
      else if(thisGen->DeltaR(*subljet)<0.8) { 
	n_neutrinos_subl++;
	neutrinos_subl_p4 += *thisGen;
      }
      else continue;
      
    }

    if(n_neutrinos_lead > 0) {
      h_leadJet_neuP4->Fill(leadjet->Pt(),neutrinos_lead_p4.Pt());
      p1->Fill(leadjet->Pt(),neutrinos_lead_p4.Pt());
      h_leadPR_neuP4->Fill(fatjetPRmass[lead],neutrinos_lead_p4.Pt());
      p3->Fill(fatjetPRmass[lead],neutrinos_lead_p4.Pt());
    } 
    h_leadJet_nNeu->Fill(leadjet->Pt(),n_neutrinos_lead);
    h_leadPR_nNeu->Fill(fatjetPRmass[lead],n_neutrinos_lead);

    if(n_neutrinos_subl > 0) {
      h_sublJet_neuP4->Fill(subljet->Pt(),neutrinos_subl_p4.Pt());
      p2->Fill(subljet->Pt(),neutrinos_subl_p4.Pt());
      h_sublPR_neuP4->Fill(fatjetPRmass[subl],neutrinos_subl_p4.Pt());
      p4->Fill(fatjetPRmass[subl],neutrinos_subl_p4.Pt());
    }
    h_sublJet_nNeu->Fill(subljet->Pt(),n_neutrinos_subl);
    h_sublPR_nNeu->Fill(fatjetPRmass[subl],n_neutrinos_subl);
    

  } // end of loop over entries
  
  TFile* outfile = new TFile(Form("%s_profile.root",name.data()),"recreate");
    p1->Write("lead_pt");
    p2->Write("subl_pt");
    p3->Write("lead_pr");
    p4->Write("subl_pr");
    outfile->Write();
  
    /*
  setNCUStyle();

  
  TCanvas* c0 = new TCanvas("c0","c0",0,0,750,550); 
  //p1->SetLineWidth(2);
  //p1->SetLineColor(kRed-2);//pyhtia
  p1->SetLineColor(kOrange+7);//herwig
  p2->SetXTitle("Jet p_{T} [GeV]");
  p2->SetYTitle("Average of vector-summed neutrino p_{T} [GeV]");
  p2->SetMarkerStyle(21);
  p2->SetMarkerSize(0.8);
  //p1->SetMarkerColor(kRed-2);
  p1->SetMarkerColor(kOrange+7);
  p1->SetMarkerStyle(22);
  p1->SetMarkerSize(0.8);
  //p2->SetMarkerColor(kGreen-3); //pythia
  p2->SetMarkerColor(kTeal-6); //herwig;
  p2->Draw();
  //p2->SetLineWidth(2);
  //p2->SetLineColor(kGreen-3); //pythia
  p2->SetLineColor(kTeal-6); //herwig
  p1->Draw("same");

  TLegend *legend1 = new TLegend(0.66, 0.80, 0.79, 0.88);
  //legend1->SetHeader("M_{BulkG} = 1 TeV");
  legend1->SetBorderSize(0);
  legend1->SetFillColor(0);
  legend1->SetFillStyle(0);
  legend1->SetTextSize(0.03);
  legend1->AddEntry(p1,"Leading Jet","lep");
  legend1->AddEntry(p2,"Subleading Jet", "lep");
  legend1->Draw();

  TCanvas* a = new TCanvas("a","a",0,0,750,550);
  //p3->SetLineWidth(2);
  p3->SetLineColor(kOrange+2);//pythia
  //p3->SetLineColor(kOrange-2);//herwig
  p4->SetXTitle("Jet pruned mass [GeV]");
  p4->SetYTitle("Average of vector-summed neutrino p_{T} [GeV]");
  p3->SetMarkerStyle(23);
  p3->SetMarkerSize(0.8);
  p3->SetMarkerColor(kOrange+2); //pythia
  //p3->SetMarkerColor(kOrange-2); //herwig
  p4->SetMarkerStyle(20);
  p4->SetMarkerSize(0.8);
  p4->SetMarkerColor(kBlue-2);//pythia
  //p4->SetMarkerColor(kAzure+5);//herwig
  p4->Draw();
  //p4->SetLineWidth(2);
  p4->SetLineColor(kBlue-2);//pythia
  //p4->SetLineColor(kAzure+5);//herwig
  p3->Draw("same");

  TLegend *legend2 = new TLegend(0.66, 0.80, 0.79, 0.88);
  legend2->SetBorderSize(0);
  legend2->SetFillColor(0);
  legend2->SetFillStyle(0);
  legend2->SetTextSize(0.03);
    //legend2->SetHeader("M_{BulkG} = 3.5 TeV");
  legend2->AddEntry(p3,"Leading Jet","lep");
  legend2->AddEntry(p4,"Subleading Jet", "lep");
  legend2->Draw();
    */

  /*
  TStyle *gStyle;

  TAxis *axis1 = h_leadJet_neuP4->GetYaxis();
  TAxis *axis2 = h_sublJet_neuP4->GetYaxis();
  TAxis *axis3 = h_leadPR_neuP4->GetYaxis();
  TAxis *axis4 = h_sublPR_neuP4->GetYaxis();
  
  if(isHerwigpp) gStyle->SetPalette(55);
  else gStyle->SetPalette(1);

  TCanvas* c1 = new TCanvas("c1","c1",0,0,750,750);
  c1->cd();
  c1->SetLeftMargin(0.105);
  c1->SetRightMargin(0.125);
  h_leadJet_neuP4->Draw("colz");
  h_leadJet_neuP4->SetStats(0);
  h_leadJet_neuP4->SetXTitle(" Leading Jet p_{T} [GeV]");
  h_leadJet_neuP4->SetYTitle("Neutrino p_{T} [GeV]");
  axis1->SetTitleOffset(1.4);

  c1->Print(Form("%s_leadPT_neupt.pdf",name.data()));

  TCanvas* c2 = new TCanvas("c2","c2",0,0,750,750);
  c2->cd();
  c2->SetLeftMargin(0.105);
  c2->SetRightMargin(0.125);
  h_leadJet_nNeu->Draw("colz");
  h_leadJet_nNeu->SetStats(0);
  h_leadJet_nNeu->SetXTitle("Leading Jet p_{T} [GeV]");
  h_leadJet_nNeu->SetYTitle("No. of neutrino");

  c2->Print(Form("%s_leadPT_nneu.pdf",name.data()));

  TCanvas* c3 =  new TCanvas("c3","c3",0,0,750,750);
  c3->cd();
  c3->SetLeftMargin(0.105);
  c3->SetRightMargin(0.125);
  h_sublJet_neuP4->Draw("colz");
  h_sublJet_neuP4->SetStats(0);
  h_sublJet_neuP4->SetXTitle("Subleading Jet p_{T} [GeV]");
  h_sublJet_neuP4->SetYTitle("Neutrino p_{T} [GeV]");
  axis2->SetTitleOffset(1.4);

  c3->Print(Form("%s_sublPT_neupt.pdf",name.data()));

  TCanvas* c4 = new TCanvas("c4","c4",0,0,750,750);
  c4->cd();
  c4->SetLeftMargin(0.105);
  c4->SetRightMargin(0.125);
  h_sublJet_nNeu->Draw("colz");
  h_sublJet_nNeu->SetStats(0);
  h_sublJet_nNeu->SetXTitle("Subleading Jet p_{T} [GeV]");
  h_sublJet_nNeu->SetYTitle("No. of neutrino");

  c4->Print(Form("%s_sublPT_nneu.pdf",name.data()));

  TCanvas* c5 = new TCanvas("c5","c5",0,0,750,750);
  c5->cd();
  c5->SetLeftMargin(0.105);
  c5->SetRightMargin(0.125);
  h_leadPR_neuP4->Draw("colz");
  h_leadPR_neuP4->SetStats(0);
  h_leadPR_neuP4->SetXTitle("Leading Jet pruned mass [GeV]");
  h_leadPR_neuP4->SetYTitle("Neutrino p_{T} [GeV]");
  axis3->SetTitleOffset(1.4);

  c5->Print(Form("%s_leadPR_neupt.pdf",name.data()));
  
  TCanvas* c6 = new TCanvas("c6","c6",0,0,750,750);
  c6->cd();
  c6->SetLeftMargin(0.105);
  c6->SetRightMargin(0.125);
  h_leadPR_nNeu->Draw("colz");
  h_leadPR_nNeu->SetStats(0);
  h_leadPR_nNeu->SetXTitle("Leading Jet pruned mass [GeV]");
  h_leadPR_nNeu->SetYTitle("No. of neutrino");

  c6->Print(Form("%s_leadPR_nneu.pdf",name.data()));
  
  TCanvas* c7 = new TCanvas("c7","c7",0,0,750,750);                                                                                               
  c7->cd();              
  c7->SetLeftMargin(0.105);
  c7->SetRightMargin(0.125);                                                                                                                                   
  h_sublPR_neuP4->Draw("colz");                                                                                                                             
  h_sublPR_neuP4->SetStats(0);                                                                                                                              
  h_sublPR_neuP4->SetXTitle("Subleading pruned mass [GeV]");                                                                                
  h_sublPR_neuP4->SetYTitle("Neutrino p_{T} [GeV]");
  axis4->SetTitleOffset(1.4);                                                       

  c7->Print(Form("%s_sublPR_neupt.pdf",name.data())); 
  
  TCanvas* c8 = new TCanvas("c8","c8",0,0,750,750);                                                                                                  
  c8->cd();            
  c8->SetLeftMargin(0.105);                                                                                                                                     
  c8->SetRightMargin(0.125);
  h_sublPR_nNeu->Draw("colz");                                                                                                                            
  h_sublPR_nNeu->SetStats(0);                                                                                                                               
  h_sublPR_nNeu->SetXTitle("Subleading Jet pruned mass [GeV]");                                                  
  h_sublPR_nNeu->SetYTitle("No. of neutrino");

  c8->Print(Form("%s_sublPR_nneu.pdf",name.data()));
  */
}
