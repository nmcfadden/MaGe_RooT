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
  TString fileName = TString("MaGeAnaOutput")+TString(".root");
  TFile *outFileMaGe = new TFile(fileName,"recreate");
  
  TNtuple *ntpleLAr = new TNtuple("primaryLAr","primaryLAr","scint:wls:ceren:edep:nDetected:nDetectedOpticalMap");
  //TH2D * hRZTrackWeight = new TH2D("RZTrackWeights","RZTrackWeight",
  Double_t totalEvents = 0;
  Int_t noHitcounter = 0;
 
  TH3D *hMap;
  TH2D *hYZMap;
  //TString mapDir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
  TString mapDir = "";//"/home/nmcfadden/RooT/MaGe_RooT/root/";
  //TString mapFileName = "OpticalMapBACON.1e9.5mm";
  //TString mapFileName = "OpticalMapBACONArgon.1e10.5mm";
  //TString mapFileName = "OpticalMapBACONXenon.1e10.5mm";
  //TString mapFileName = "OpticalMapLEGEND200.4e9_10mm";
  //TString mapFileName = "Detect.113005-20181019";
  TString mapFileName = "BACoN-heat.PMTs.10000";
  TFile* mapFile = TFile::Open(mapDir+mapFileName+TString(".root"));
  if(mapFile != NULL){
    mapFile->GetObject("OpticalMap",hMap);
    mapFile->GetObject("2DOpticalMap_YZ",hYZMap);
  }

  cout<<"starting run"<<endl; 
  
  for(int k = 0; k < 1 ; k++){
   
    //TString dir = "/home/nmcfadden/XenonDoping/bin/Linux-g++/";
    TString dir = "/home/nmcfadden/MaGe/bin/Linux-g++/";
    TString fileName = "";
    //fileName = "100keVe-";
    fileName = "1MeValpha";
    ///*
    if(!fileExist(string(dir+fileName+TString(".root")))){
      cout<<"processed "<<k<<" files"<<endl;
      break;
    }
    //*/
    TFile* infile = TFile::Open(dir+fileName+TString(".root"));
    string neventsString = infile->Get("NumberOfEvents")->GetTitle();
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
      Double_t x = primaries->GetX(),y = primaries->GetY(),z = primaries->GetZ();//,time = primaries->GetT();
      Double_t r = sqrt(x*x+y*y);
      Double_t theta = std::acos(x/r);
      Double_t nPhotons = 0,eDepLAr = 0,eDepTotal = 0;
      Int_t stepCounter = 0,nDetected = 0,scintNum = 0,wlsNum = 0,cereNum = 0;
      bool hitLAr = false,hitGe = false,hitGAr = false;
      if( y < 0) theta += 3.14159265359;
      for (Int_t j = 0; j < eventSteps->GetNSteps();j++){
        Double_t edep = step->GetEdep();
        step = eventSteps->GetStep(j);
        physName = step->GetPhysVolName();
        TString procName = step->GetProcessName();
        Int_t trackWeight = step->GetTrackWeight();
        Int_t trackID = step->GetTrackID();
        //Count optical photons
        //Scint 1, OpWLS 2, Cerenkov 3
        if(trackWeight == 1 && step->GetStepNumber() == 1) {
          scintNum++;
        }
        if(trackWeight == 2 ) {
          wlsNum++;
        }
        if(trackWeight == 3 ) {
          cereNum++;
        }
        //count particles that are not optical photons (trackWeight == 0) in the detector
        //Note "Detector" is LAr
        if(physName == "Detector" && trackWeight == 0){
          eDepLAr += edep;
        }
        //argon scint can deposit energy in the SiPM though that is not real detection
        if(physName.Contains("SiPM") && trackWeight == 2){
          nDetected++;
        }

        if(mapFile!= NULL){
          Int_t bin = hMap->FindBin(step->GetX(),step->GetY(),step->GetZ());
          Double_t mapProb = hMap->GetBinContent(bin);
          Double_t scintYeild = 40.;//40 photons/keV
          Double_t eThresh = mapProb*scintYeild;
          //average photon yeild accord to map
          nPhotons += edep*eThresh;
        }

      }
      ntpleLAr->Fill(scintNum,wlsNum,cereNum,eDepLAr,nDetected,nPhotons);
      
    }
  }

  outFileMaGe->Write();
  outFileMaGe->Close();

  cout<<"Total events = "<<totalEvents<<", no LAr hits counter = "<<noHitcounter<<endl;
  cout<<"root -l "<<fileName<<endl;
  return 0;
}
