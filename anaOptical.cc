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
  TString fileName = TString("MaGe.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root");
  TFile *outFileMaGe = new TFile(fileName,"recreate");
//  TFile *outFileOpticalMap = new TFile(TString("OpticalMap.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root"),"recreate");
  
  TNtuple *ntpleStepSiPM = new TNtuple("stepSiPM","stepSiPM","event:nstep:edep:ke:x:y:z:r:t");
  TNtuple *ntpleSiPM = new TNtuple("primarySiPM","primarySiPM","x:y:z:r:theta:px:py:pz:nsteps");
  TNtuple *ntpleStepFiberCore = new TNtuple("stepFireCore","stepFireCore","event:nstep:edep:ke:x:y:z:t:t");
  TNtuple *ntpleFiberCore = new TNtuple("primaryFiber","primaryFiber","x:y:z:r:theta:px:py:pz:nsteps");

  std::map<std::vector<Int_t>,Double_t> prob_map;
  //Double_t gridSpacing = 2;//*mm
  Double_t gridSpacing = 5;//*mm
  //Double_t maxX = 320.,maxY=320.0,maxZ=600,minZ = -200.,minR = 0,maxR=320;
  bool twoD = false;//true;
  //Double_t maxX = 1000.,maxY=1000,maxZ=1200,minZ = -800.,minR = 320., maxR = 1000.;
  //Double_t maxX = 300.,minX = 0,maxY=300.,minY = 0,maxZ=925,minZ = -425.,minR = 0., maxR = 300.;
  Double_t maxX = 700.,minX = 0.,maxY=700.,minY = 0,maxZ=925,minZ = -425.,minR = 300., maxR = 700.;
  Int_t nbinsX = 2*(maxX-minX)/gridSpacing,nbinsY = 2*(maxY-minY)/gridSpacing,nbinsZ = (maxZ-minZ)/gridSpacing;
  TH3D* hMap = new TH3D("OpticalMap","OpticalMap",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  TH3D* hMapUnscaled = new TH3D("OpticalMap_unScaled","OpticalMap_unScaled",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);

  TH2D* h2DMapRZ = new TH2D("2DOpticalMap_RZ","2DOpticalMap_RZ",nbinsX/2.,minR,maxR,nbinsZ,minZ,maxZ);
  TH2D* h2DMapRZUnscaled = new TH2D("2DOpticalMap_RZUnscaled","2DOpticalMap_RZUnscaled",nbinsX/2.,minR,maxR,nbinsZ,minZ,maxZ);

  TH2D* h2DMapXY = new TH2D("2DOpticalMap_XY","2DOpticalMap_XY",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);
  TH2D* h2DMapXYUnscaled = new TH2D("2DOpticalMap_XYUnscaled","2DOpticalMap_XYUnscaled",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);

  TH2D* h2DMapYZ = new TH2D("2DOpticalMap_YZ","2DOpticalMap_YZ",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  TH2D* h2DMapYZUnscaled = new TH2D("2DOpticalMap_YZUnscaled","2DOpticalMap_YZUnscaled",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);  
  Double_t totalEvents = 0;
  Int_t n3D = 0, nXY = 0, nYZ = 0, nRZ = 0;
  cout<<"Histogram has "<<nbinsX<<" "<<nbinsY<<" "<<nbinsZ<<" "<<nbinsX*nbinsY*nbinsZ <<" bins"<<endl;

  cout<<"starting run"<<endl; 
  TFile * infile;
  for(int k =0; k <= 100 ; k++){
    //TODO--change input dir
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/MaGe/bin/Linux-g++/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/BACoN/bin/Linux-g++/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden//LEGEND_200_One_FIBER_Array/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/LGND_200Orig100M1.1mAttenuation/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/array/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/arrOpticalDist/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/LGND19stringNoHV/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/ExteriorMap/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/opticalArray/";
    TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/19StringRadon700/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/14StringRadon700/";
    //TString fileName = "SensitiveVolumesLGND_200Alt1" +to_string(k);
    //TString fileName = "SensitiveVolumesLGND_200Orig1.1mAttenuation" +to_string(k);
    //TString fileName = "SensitiveVolumesLGND_200Orig" +to_string(k);
    //TString fileName = "SensitiveVolumes" +to_string(k);
    //TString fileName = "SensitiveVolumesLGND_200Orig"+to_string(k);
    //TString fileName = "RDMiso224.88.Optical";
    //TODO -- change fileName
    //TString fileName = "ExteriorMap_"+to_string(k);
    //TString fileName = "OpticalRunExterior"+to_string(k);
    //TString fileName = "OpticalRun14String"+to_string(k); 
    TString fileName = "OpticalRun19StringExterior"+to_string(k); 
    //TString fileName = "OpticalRun14StringExterior"+to_string(k); 
    //TString fileName = "OpticalRun"+to_string(k);
    //TString fileName = "SpeedTest";
    cout<<"File location at: "<<dir+fileName+TString(".root")<<endl;;
    ///*
    if(!fileExist(string(dir+fileName+TString(".root")))){
      cout<<"processed "<<k<<" files"<<endl;
      continue;
      //break;
    }
    //*/
    //TFile* infile = TFile::Open(dir+fileName+TString(".root"));
    infile = new TFile(dir+fileName+TString(".root"));
    if(infile->IsZombie()){
      cout<<"File "<<k<<" is Zombie"<<endl;
      continue;
    }
    else
      infile = TFile::Open(dir+fileName+TString(".root"));
    if(infile->Get("NumberOfEvents") == NULL) {
      cout<<"file has not closed yet"<<endl;
      continue;
    }
    string neventsString = infile->Get("NumberOfEvents")->GetTitle();
    totalEvents += std::stoi(neventsString,nullptr,10);
    infile->Close();
    delete infile;
    outFileMaGe->cd();
    
    //TFile *infile = new TFile(dir+fileName,"READONLY");
    TChain *fTree = new TChain("fTree");
    fTree->Add(dir+fileName+TString(".root"));
    
    Long64_t nentries = (Long64_t)fTree->GetEntries();
    
    //cout<<"File location at: "<<dir+fileName+TString(".root")<<". Processing "<<nentries<<" entries"<<endl;
    cout<<". Processing "<<nentries<<" entries"<<endl;
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
      if((i+1)%10000 == 0 || i == nentries - 1 ) cout<<"\tprocessed "<<i+1<<" events"<<endl;
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
      bool hitSiPM = false,hitFiber = false;
      Int_t hitSiPMCounter = 0,hitFiberCounter = 0;
      if( y < 0) theta += 3.14159265359;
      for (Int_t j = 0; j < eventSteps->GetNSteps();j++){
        step = eventSteps->GetStep(j);
        physName = step->GetPhysVolName();
        //cout<<physName<<endl;
        if(physName.Contains("SiPM") && step->GetEdep() > 0 ){
          hitSiPM = true;
          hitSiPMCounter++;
          ntpleStepSiPM->Fill(i,j,step->GetEdep(),step->GetKineticE(),step->GetX(),step->GetY(),step->GetZ(),sqrt(step->GetX()*step->GetX() +step->GetY()*step->GetY()),step->GetT());
          ntpleSiPM->Fill(x,y,z,r,theta,px,py,pz,hitSiPMCounter);
          hMapUnscaled->Fill(x,y,z);
          h2DMapRZUnscaled->Fill(r,z);
          h2DMapXYUnscaled->Fill(x,y);
          h2DMapYZUnscaled->Fill(y,z);
        }
        else if(physName.Contains("FiberCore")){
          hitFiber = true;
          hitFiberCounter++;
          ntpleStepFiberCore->Fill(i,j,step->GetEdep(),step->GetKineticE(),step->GetX(),step->GetY(),step->GetZ(),sqrt(step->GetX()*step->GetX() +step->GetY()*step->GetY()),step->GetT());
        }
      }
      /*
      if(hitSiPM){
        ntpleSiPM->Fill(x,y,z,r,theta,px,py,pz,hitSiPMCounter);
        hMapUnscaled->Fill(x,y,z);
        h2DMapRZUnscaled->Fill(r,z);
        h2DMapXYUnscaled->Fill(x,y);
        h2DMapYZUnscaled->Fill(y,z);
      }
      */
      if(hitFiber){
        ntpleFiberCore->Fill(x,y,z,r,theta,px,py,pz,hitFiberCounter);
      }
    }
    delete step;
    delete primaries;
    delete fTree;
  }
  if(!twoD){
  ///*
  //weight 3Dhistogram
  for(int i = 0; i <= hMap->GetNbinsX();i++){
    for(int j = 0; j <= hMap->GetNbinsY();j++){
      for(int k = 0; k <= hMap->GetNbinsZ();k++){
        Double_t binVal = hMapUnscaled->GetBinContent(i,j,k);
        Double_t x = gridSpacing*(i)-maxX;
        Double_t y = gridSpacing*(j)-maxY;
        Double_t z = gridSpacing*(k)+minZ;
        if(binVal == 0){
            //binVal = hMapUnscaled->Interpolate(x,y,z);
            n3D++;
        }
        hMap->SetBinContent(i,j,k,binVal/totalEvents);
      }
    }
  }
  cout<<"finished XYZ Histogram"<<endl;
  for(int i = 0; i <= h2DMapXY->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapXY->GetNbinsY();j++){
      Double_t binVal = h2DMapXYUnscaled->GetBinContent(i,j);
      Double_t x = gridSpacing*(i)-maxX;
      Double_t y = gridSpacing*(j)-maxY;
      if(binVal == 0){
         //binVal = h2DMapXYUnscaled->Interpolate(x,y);
         nXY++;
      }
      h2DMapXY->SetBinContent(i,j,binVal/totalEvents);
    }
  }
  cout<<"finished XY Histogram"<<endl;
  for(int i = 0; i <= h2DMapYZ->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapYZ->GetNbinsY();j++){
      Double_t binVal = h2DMapYZUnscaled->GetBinContent(i,j);
      Double_t y = gridSpacing*(i)-maxY;
      Double_t z = gridSpacing*(j)+minZ;
      if(binVal == 0){
        //binVal = h2DMapYZUnscaled->Interpolate(y,z);
        nYZ++;
      }
      h2DMapYZ->SetBinContent(i,j,binVal/totalEvents);
    }
  }
  }
  for(int i = 0; i <= h2DMapRZ->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapRZ->GetNbinsY();j++){
      Double_t binVal = h2DMapRZUnscaled->GetBinContent(i,j);
      Double_t r = gridSpacing*(i)+minR;
      Double_t z = gridSpacing*(j)+minZ;
      if(binVal == 0){
        //binVal = h2DMapRZUnscaled->Interpolate(r,z);
        nRZ++; 
      }
    h2DMapRZ->SetBinContent(i,j,binVal/totalEvents);
    }
  }
  cout<<"finished RZ Histogram"<<endl;

  //*/
  cout<<"finished YZ Histogram"<<endl;
 //  outFileMaGe->cd();
  //ntpleStep->Write();
  //ntpleSiPM->Write();
  //ntpleFiberCore->Write();
  outFileMaGe->Write();
  outFileMaGe->Close();
/*
  outFileOpticalMap->cd();
  hMap->Write();
  hMapUnscaled->Write();
  outFileOpticalMap->Close();
*/
  cout<<"Total events = "<<totalEvents<<endl;
  cout<<"Histogram has "<<nbinsX<<" "<<nbinsY<<" "<<nbinsZ<<" "<<nbinsX*nbinsY*nbinsZ <<" bins"<<endl;
  cout<<"Extropalted "<<100.*n3D/(nbinsX*nbinsY*nbinsZ)<<" 3D ... "<<100.*nXY/(nbinsX*nbinsY)<<" XY... "<<100.*nYZ/(nbinsY*nbinsZ)<<" YZ... "<<100.*nRZ/(nbinsZ*nbinsX/2.)<<" RZ "<<endl;
  cout<<"root -l "<<fileName<<endl;
  return 0;
}
