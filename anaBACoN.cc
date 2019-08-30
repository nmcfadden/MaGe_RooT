#include "MGTMCEventSteps.hh"
#include "MGTMCStepData.hh"
#include "MGTMCRun.hh"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TProof.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TObject.h"
#include "TString.h"
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <getopt.h>
#include <TDatime.h>
#include <TNtuple.h>
#include <unordered_map>
#include <map>

using namespace std;

inline bool fileExist (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


int main(){

  TDatime time;
  //time 12:36:26 133626
  ///date 24/12/1997 19971224
  TString fileName = TString("BACoNAnaOutput")+TString(".root");
  TFile *outFileMaGe = new TFile(fileName,"recreate");
  
  TNtuple *ntpleLAr = new TNtuple("primaryLAr","primaryLAr","scint:wls:ceren:edep:nDetected:nDetectedOpticalMap:XenonDoped:uPanel1:uPanel2:muonTrigger");
  TNtuple *ntupleStep = new TNtuple("step","step","parentID:trackID:particleID:edep:KE:stepLength:x:y:z");
  TNtuple *ntpleWLS = new TNtuple("wls","wls","px:py:pz");
  TNtuple *ntplePrimary = new TNtuple("primary","primary","x:y:z:r:theta:px:py:pz:KE:muonTrigger");

  //TH2D * hRZTrackWeight = new TH2D("RZTrackWeights","RZTrackWeight",
  Double_t totalEvents = 0;
 
  TH3D *hMap;
  TH3D *hMapXe;
  TH2D *hYZMap;
  //TString mapDir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
  TString mapDir = "/home/nmcfadden/BACoNSims/root/";
  //TString mapFileName = "OpticalMapBACON.1e9.5mm";
  //TString mapFileName = "OpticalMapBACONArgon.1e10.5mm";
  //TString mapFileName = "OpticalMapBACONXenon.1e10.5mm";
  //TString mapFileName = "OpticalMapLEGEND200.4e9_10mm";
  //TString mapFileName = "Detect.113005-20181019";
  TString mapFileName = "BACoN-heat.PMTs.10000.CorrectedPMTQE";
  TFile* mapFile = TFile::Open(mapDir+mapFileName+TString(".root"));
  if(mapFile != NULL){
    mapFile->GetObject("OpticalMap",hMap);
    mapFile->GetObject("2DOpticalMap_YZ",hYZMap);
  }
  mapFileName = "BACoN-heat.PMTs.10000.CorrectedPMTQE.XeDoped";
  mapFile = TFile::Open(mapDir+mapFileName+TString(".root"));
  if(mapFile != NULL){
    mapFile->GetObject("OpticalMap",hMapXe);
  }

  cout<<"starting run"<<endl; 
  Int_t pmtHitCounter = 0;
  for(int k = 0; k < 19; k++){
   
    //TString dir = "/home/nmcfadden/XenonDoping/bin/Linux-g++/";
    //TString dir = "/home/nmcfadden/RooT/MaGe_RooT/root/";
    TString dir = "/home/nmcfadden/BACoNSims/cosmicMuons/";
    //TString dir = "/home/nmcfadden/BACoNSims/Co60/";
    TString fileName = "";
    //fileName = "0.5MeVGammaEvents.29.0.0.cm";
    //fileName = "pmtVUV";
    //fileName = "outputGamma";
    fileName = "cosmicMuon"+to_string(k);
    //fileName = "Co60"+to_string(k);
    ///*
    if(!fileExist(string(dir+fileName+TString(".root")))){
      cout<<"processed "<<k<<" files"<<endl;
      break;
    }
    //*/
    TFile* infile = TFile::Open(dir+fileName+TString(".root"));
    string neventsString;
    if(infile->Get("NumberOfEvents") == NULL){
      cout<<"Number of event Tag is null"<<endl;
      neventsString = "0";
    }
    else{
      neventsString = infile->Get("NumberOfEvents")->GetTitle();
    }
    totalEvents += std::stoi(neventsString,nullptr,10);
    infile->Close();
    outFileMaGe->cd();
    
    //TFile *infile = new TFile(dir+fileName,"READONLY");
    TChain *fTree = new TChain("fTree");
    fTree->Add(dir+fileName+TString(".root"));
    
    Long64_t nentries = (Long64_t)fTree->GetEntries();
    
    cout<<"File location at: "<<dir+fileName+TString(".root")<<". Processing "<<nentries<<" entries"<<endl;
    MGTMCEventSteps *eventSteps = 0;
    MGTMCEventSteps *eventPrimaries = 0;
    if(fTree != NULL){
      fTree->SetBranchAddress("eventSteps",&eventSteps);
      fTree->SetBranchAddress("eventPrimaries",&eventPrimaries);
    }
    else{
      cout<<"NULL fTree"<<endl;
      break;
    }
       
    const MGTMCStepData *step,*primaries;
   
    
    for(Int_t i = 0; i < nentries ; i++){
      //if((i+1)%100 == 0 || i == nentries - 1 ) cout<<"\tprocessed "<<i+1<<" events"<<endl;
      if((i+1)%(int(nentries*.1)) == 0 || i == nentries - 1 ) cout<<"\tprocessed "<<i+1<<" events"<<endl;
      fTree->GetEntry(i);
      TString physName;
      primaries = eventPrimaries->GetStep(0);
      if(primaries == NULL){
        cout<<"null primary"<<endl;
        continue;
      }
      ///*
      Double_t x = primaries->GetX(),y = primaries->GetY(),z = primaries->GetZ();//,time = primaries->GetT();
      Double_t r = sqrt(x*x+y*y);
      Double_t theta = std::acos(x/r);
      if( y < 0) theta += 3.14159265359;
      //*/
      Double_t nPhotons = 0,nPhotonsXe = 0,eDepLAr = 0,eDepTotal = 0,eDepPanel1 = 0,eDepPanel2 = 0;
      Int_t nDetected = 0,scintNum = 0,wlsNum = 0,cereNum = 0,muonTrigger = -1;
      bool hitLAr = false,hitPMT = false,muonFlag1 = false,muonFlag2 = false;
      //472,630,711,841,887,997
      for (Int_t j = 0; j < eventSteps->GetNSteps();j++){
        step = eventSteps->GetStep(j);
        Double_t edep = step->GetEdep();
        eDepTotal += edep;
        physName = step->GetPhysVolName();
        TString procName = step->GetProcessName();
        Int_t trackWeight = step->GetTrackWeight();
        Int_t trackID = step->GetTrackID();
        Int_t stepNum = step->GetStepNumber();
        Int_t parentID = step->GetParentTrackID();
        Int_t particleID = step->GetParticleID();
        Double_t KE = step->GetKineticE();
        Double_t stepX = step->GetX(),stepY = step->GetY(),stepZ = step->GetY();
        Double_t stepLength = step->GetStepLength();

        if(physName == "Detector") ntupleStep->Fill(parentID,trackID,particleID,edep,KE,stepLength,stepX,stepY,stepZ);
        //cout<<"Event "<<i<<", trackID "<<trackID<<", stepNumber "<<stepNum<<", KE "<<KE<<", eDep "<<edep<<
        //", "<<physName<<", processName "<<procName<<" ("<<step->GetX()<<","<<step->GetY()<<","<<step->GetZ()<<")"<<endl;
        //Count optical photons
        //Make sure write all steps is on
        //Scint 1, OpWLS 2, Cerenkov 3
        if(trackWeight == 1 && stepNum == 1) {
          scintNum++;
        }
        if(trackWeight == 2 && stepNum == 1) {
          wlsNum++;
          ntpleWLS->Fill(step->GetPx(),step->GetPy(),step->GetPz());
          //cout<<"trackID "<<trackID<<", stepNumber "<<stepNum<<", KE "<<KE<<", eDep "<<edep<<endl;
        }
        if(trackWeight == 3 ) {
          cereNum++;
        }
        //count particles that are not optical photons (trackWeight == 0) in the detector
        //Note "Detector" is LAr
        if(physName == "Detector" && trackWeight == 0){
          eDepLAr += edep;
          hitLAr = true;
        }
        //argon scint can deposit energy in the SiPM though that is not real detection
        if(physName.Contains("physicalPMT") && trackWeight == 2 && edep > 0){
          nDetected++;
          hitPMT = true;
        }
        if(physName == "ScintPanel1" && edep > 0 && parentID == 1){
          muonFlag1 = true;
          eDepPanel1 += edep;
        }
        if(physName == "ScintPanel2" && edep > 0 && parentID == 1){
          muonFlag2 = true;
          eDepPanel2 += edep;
        }

        if(mapFile!= NULL && trackWeight == 0){
          Int_t bin = hMap->FindBin(step->GetX(),step->GetY(),step->GetZ());
          Double_t mapProb = hMap->GetBinContent(bin);
          Double_t scintYeild = 42370;
          Double_t eThresh = mapProb*scintYeild;//*(1200000.)*(1200000.)/(1e8*120);
          if(mapProb == 0) continue;
          //average photon yeild accord to map
          nPhotons += edep*eThresh;
          //cout<<nPhotons<<"...nPhotons = edep*(mapProb*scintYield)..."<<edep*eThresh<<" = "<<edep<<"*("<<mapProb<<"*"<<scintYeild<<")"<<endl;
          //Xenon Doped Map
          mapProb = hMapXe->GetBinContent(bin);
          eThresh = mapProb*scintYeild;
          nPhotonsXe += edep*eThresh;
          //cout<<"xenon "<<nPhotonsXe<<"...nPhotons = edep*(mapProb*scintYield)..."<<edep*eThresh<<" = "<<edep<<"*("<<mapProb<<"*"<<scintYeild<<")"<<endl;
        }

      }
      if(hitPMT) pmtHitCounter++;
      if(muonFlag1 && muonFlag2 )muonTrigger = 1;
      ntplePrimary->Fill(x,y,z,r,theta,primaries->GetPx(),primaries->GetPy(),primaries->GetPz(),primaries->GetKineticE(),muonTrigger);
      if(eDepLAr > 0) ntpleLAr->Fill(scintNum,wlsNum,cereNum,eDepLAr,nDetected,nPhotons,nPhotonsXe,eDepPanel1,eDepPanel2,muonTrigger);
      
    }
  }

  outFileMaGe->Write();
  outFileMaGe->Close();

  cout<<"Total events = "<<totalEvents<<", pmt hits = "<<pmtHitCounter<<endl;
  cout<<"root -l "<<fileName<<endl;
  return 0;
}
