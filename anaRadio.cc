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
  TString fileName = TString("Radio.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root");
  TFile *outFileMaGe = new TFile(fileName,"recreate");
  
  TNtuple *ntpleStepLAr = new TNtuple("stepLAR","LArCore","event:nstep:edep:ke:x:y:z:r:t");
  TNtuple *ntpleLAr = new TNtuple("primaryLAr","primaryLAr","x:y:z:r:theta:px:py:pz:nPhotons:edep:mapThresh:pmtCounter:totalPhotons");
  TNtuple *ntpleGAr = new TNtuple("primaryGAr","primaryGAr","x:y:z:r:theta:px:py:pz:nPhotons:edep:mapThresh:pmtCounter:totalPhotons");
  TNtuple *ntpleGe = new TNtuple("primaryGe","primaryGe","x:y:z:r:theta:px:py:pz:nPhotons:edep:mapThresh:pmtCounter:totalPhotons");
  //TH2D * hRZTrackWeight = new TH2D("RZTrackWeights","RZTrackWeight",
  Double_t totalEvents = 0;
  Int_t noHitcounter = 0;
 
  TH3D *hMap;
  TH2D *hYZMap;
  //TString mapDir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
  TString mapDir = "/home/nmcfadden/RooT/MaGe_RooT/";
  //TString mapFileName = "OpticalMapBACON.1e9.5mm";
  //TString mapFileName = "OpticalMapBACONArgon.1e10.5mm";
  //TString mapFileName = "OpticalMapBACONXenon.1e10.5mm";
  //TString mapFileName = "OpticalMapLEGEND200.4e9_10mm";
  //TString mapFileName = "Detect.113005-20181019";
  TString mapFileName = "OpticalMapBACoN.1e10";
  TFile* mapFile = TFile::Open(mapDir+mapFileName+TString(".root"));
  mapFile->GetObject("EnergyMap",hMap);
  mapFile->GetObject("2DOpticalMap_YZ",hYZMap);

  cout<<"starting run"<<endl; 
  
  for(int k = 0; k < 1 ; k++){
   
    TString dir = "/home/nmcfadden/XenonDoping/bin/Linux-g++/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/MaGe/bin/Linux-g++/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/BACoN/bin/Linux-g++/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/arrayGenerator/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/array/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/LGND_200Orig100M1.1mAttenuation/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/arrBACON/";
    TString fileName = "OpticalRun";
    //TString fileName = "SensitiveVolumesLGND_200Alt1" +to_string(k);
    //TString fileName = "SensitiveVolumesLGND_200Orig1.1mAttenuation" +to_string(k);
    //TString fileName = "SensitiveVolumesLGND_200Orig" +to_string(k);
    //TString fileName = "SensitiveVolumes" +to_string(k);
    //TString fileName = "RadaioDecaySensitiveVolume";
    //TString fileName = "PMTCountingTest";
    //TString fileName = "1MeVBeta100000Events";
    //TString fileName = "1MeVAlphaEvents";
    //TString fileName = "1MeVBetaEvents";
    //TString fileName = "1MeVGammaEvents";
    //TString fileName = "RDMiso224.88";
    //TString fileName = "2MeVGamma.1e4";
    //TString fileName = "SensitiveVolumesLGND_200Orig"+to_string(k);
    //    TString fileName = "RadioChecks"+to_string(k);
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
      if((i+1)%100 == 0 || i == nentries - 1 ) cout<<"\tprocessed "<<i+1<<" events"<<endl;
      fTree->GetEntry(i);
      TString physName;
      primaries = eventPrimaries->GetStep(0);
      if(primaries == NULL){
        cout<<"null primary"<<endl;
        continue;
      }

      Double_t x = primaries->GetX(),y = primaries->GetY(),z = primaries->GetZ();//,time = primaries->GetT();
      Double_t px = primaries->GetPx(),py = primaries->GetPy(),pz = primaries->GetPz();
      Double_t r = sqrt(x*x+y*y);
      Double_t theta = std::acos(x/r);
      Double_t nPhotons = 0,eDepLAr = 0,eDepGe = 0,mapThresh = 0.,eDepGAr = 0;
      Int_t stepCounter = 0,pmtCounter = 0,nTotalPhotons = 0;
      bool hitLAr = false,hitGe = false,hitGAr = false;
      if( y < 0) theta += 3.14159265359;
      for (Int_t j = 0; j < eventSteps->GetNSteps();j++){
        step = eventSteps->GetStep(j);
        physName = step->GetPhysVolName();
        TString procName = step->GetProcessName();
        if(step->GetTrackWeight() == 1 && step->GetStepNumber() == 1) {
          nTotalPhotons++;
        }
        //eDep +=step->GetEdep();
        if(physName.Contains("physicalPMT")&& procName.Contains("WLS") && step->GetParticleID() == 0 && step->GetEdep() != 0) pmtCounter++;
        //if(physName.Contains("SiPM")&& step->GetParticleID() == 0 && step->GetEdep() != 0) pmtCounter++;
        if(physName =="Detector" || physName == "argonGasPhysical" ){
          //if(step->GetParticleID() == 0) continue;
          ///*
          //if(sqrt(step->GetX()*step->GetX() +step->GetY()*step->GetY()) > 254.) continue;
          //if(fabs(step->GetZ()) > 23.*2.54*10/2.) continue;
          Int_t bin = hMap->FindBin(step->GetX(),step->GetY(),step->GetZ());
          //Int_t bin = hYZMap->FindBin(step->GetY(),step->GetZ());
          Double_t eThresh = (hMap->GetBinContent(bin)/1000.); //map in keV, Geant4 is in MeV
          //Double_t eThresh = hYZMap->GetBinContent(bin)/1000.; //map in keV, Geant4 is in MeV
          //cout<<"bin "<<bin<<", threshold "<<eThresh<<", eDep "<<step->GetEdep()<<endl;
          if(eThresh == 0){
            eThresh = 100.; //100 MeV as a large threshold
            //eThresh = hMap->Interpolate(step->GetX(),step->GetY(),step->GetZ())/1000.;
            //cout<<eThresh<<", ("<<step->GetY()<<","<<step->GetZ()<<") ,"<<endl;
            //continue;
          }
          nPhotons += step->GetEdep()/eThresh;
          mapThresh+= eThresh;
          stepCounter++;
          //*/
          eDepLAr +=step->GetEdep();
          hitLAr = true;
          ntpleStepLAr->Fill(i,j,step->GetEdep(),step->GetKineticE(),step->GetX(),step->GetY(),step->GetZ(),sqrt(step->GetX()*step->GetX() +step->GetY()*step->GetY()),step->GetT());
        }
        if(physName == "argonGasPhysical"){
          Int_t bin = hMap->FindBin(step->GetX(),step->GetY(),step->GetZ());
          //Int_t bin = hYZMap->FindBin(step->GetY(),step->GetZ());
          Double_t eThresh = (hMap->GetBinContent(bin)/1000.); //map in keV, Geant4 is in MeV
          //Double_t eThresh = hYZMap->GetBinContent(bin)/1000.; //map in keV, Geant4 is in MeV
          //cout<<"bin "<<bin<<", threshold "<<eThresh<<", eDep "<<step->GetEdep()<<endl;
          if(eThresh == 0){
            eThresh = 100.; //100 MeV as a large threshold
          }
          nPhotons += step->GetEdep()/eThresh;
          mapThresh+= eThresh;
          //stepCounter++;
          //*/
          eDepGAr +=step->GetEdep();
          hitGAr = true;
        }
        if(physName.Contains("ActiveDet")){
          eDepGe += step->GetEdep();
          hitGe = true;
        }
        /*
        if(i%1000 == 0) 
          cout<<"Ethresh "<<eThresh<<", ("<<step->GetY()<<","<<step->GetZ()<<") ,"<<" eDep "<<eDep<<", mapThresh "<<mapThresh/stepCounter<<endl;
        */

      }
      /*
      cout<<"*********************"<<endl;
      cout<<" Edep "<<eDep<<" nPhotons "<<nPhotons<<endl;
      cout<<"*********************"<<endl;
      */
      if(hitLAr){
        ntpleLAr->Fill(x,y,z,r,theta,px,py,pz,nPhotons,eDepLAr,mapThresh/stepCounter,pmtCounter,nTotalPhotons);
      }
      if(hitGe){ 
        ntpleGe->Fill(x,y,z,r,theta,px,py,pz,nPhotons,eDepGe,mapThresh/stepCounter,pmtCounter,nTotalPhotons);
      }
      if(hitGAr){
        ntpleGAr->Fill(x,y,z,r,theta,px,py,pz,nPhotons,eDepGAr,mapThresh/stepCounter,pmtCounter,nTotalPhotons);
      }
      if(!hitGe && !hitLAr){
        noHitcounter++;
        //cout<<physName<<endl;
      }
    }
  }
 //  outFileMaGe->cd();
  //ntpleStep->Write();
  //ntplePMT->Write();
  //ntpleLAr->Write();
  outFileMaGe->Write();
  outFileMaGe->Close();
/*
  outFileOpticalMap->cd();
  hMap->Write();
  hMapUnscaled->Write();
  outFileOpticalMap->Close();
*/
  cout<<"Total events = "<<totalEvents<<", no Ger or LAr hit counter = "<<noHitcounter<<endl;
  cout<<"root -l "<<fileName<<endl;
  return 0;
}
