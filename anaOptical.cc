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
  
  TNtuple *ntpleStepSiPM = new TNtuple("stepSiPM","stepSiPM","event:edep:ke:x:y:z:r:primX:primY:primZ:primR:primTheta");

  std::map<std::vector<Int_t>,Double_t> prob_map;
  //Double_t gridSpacing = 2;//*mm
  Double_t gridSpacing = 5;//*mm
  //Double_t maxX = 320.,maxY=320.0,maxZ=600,minZ = -200.,minR = 0,maxR=320;
  //Double_t maxX = 1000.,maxY=1000,maxZ=1200,minZ = -800.,minR = 320., maxR = 1000.;
  //Double_t maxX = 300.,minX = 0,maxY=300.,minY = 0,maxZ=925,minZ = -425.,minR = 0., maxR = 300.;
  Double_t maxX = 300.,minX = -300.,maxY=300.,minY = -300,maxZ=860,minZ = -860,minR = 0, maxR = 300.;
  Int_t nbinsX = (maxX-minX)/gridSpacing,nbinsY = (maxY-minY)/gridSpacing,nbinsZ = (maxZ-minZ)/gridSpacing;
  TH3D* hMap = new TH3D("OpticalMap","OpticalMap",nbinsX,minX,maxX,nbinsY,minY,maxY,nbinsZ,minZ,maxZ);
  TH3D* hMapUnscaled = new TH3D("OpticalMap_unScaled","OpticalMap_unScaled",nbinsX,minX,maxX,nbinsY,minY,maxY,nbinsZ,minZ,maxZ);

  TH2D* h2DMapRZ = new TH2D("2DOpticalMap_RZ","2DOpticalMap_RZ",nbinsX/2.,minR,maxR,nbinsZ,minZ,maxZ);
  TH2D* h2DMapRZUnscaled = new TH2D("2DOpticalMap_RZUnscaled","2DOpticalMap_RZUnscaled",nbinsX/2.,minR,maxR,nbinsZ,minZ,maxZ);

  TH2D* h2DMapXY = new TH2D("2DOpticalMap_XY","2DOpticalMap_XY",nbinsX,minX,maxX,nbinsY,minY,maxY);
  TH2D* h2DMapXYUnscaled = new TH2D("2DOpticalMap_XYUnscaled","2DOpticalMap_XYUnscaled",nbinsX,minX,maxX,nbinsY,minY,maxY);

  TH2D* h2DMapYZ = new TH2D("2DOpticalMap_YZ","2DOpticalMap_YZ",nbinsY,minY,maxY,nbinsZ,minZ,maxZ);
  TH2D* h2DMapYZUnscaled = new TH2D("2DOpticalMap_YZUnscaled","2DOpticalMap_YZUnscaled",nbinsY,minY,maxY,nbinsZ,minZ,maxZ);  
  Double_t totalEvents = 0,totalPrimaries;
  Int_t n3D = 0, nXY = 0, nYZ = 0, nRZ = 0;
  cout<<"Histogram has "<<nbinsX<<" "<<nbinsY<<" "<<nbinsZ<<" "<<nbinsX*nbinsY*nbinsZ <<" bins"<<endl;

  cout<<"starting run"<<endl; 
  TFile * infile;
  for(int k =0; k < 1 ; k++){
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
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/19StringRadon700/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/14StringRadon700/";
    TString dir = "/home/nmcfadden/MaGe/bin/Linux-g++/";
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
    //TString fileName = "OpticalRun19StringExterior"+to_string(k); 
    //TString fileName = "OpticalRun14StringExterior"+to_string(k); 
    //TString fileName = "OpticalRun"+to_string(k);
    TString fileName = "SpeedTest";
    cout<<"File location at: "<<dir+fileName+TString(".root")<<endl;;
    if(!fileExist(string(dir+fileName+TString(".root")))){
      cout<<"processed "<<k<<" files"<<endl;
      continue;
      //break;
    }
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

    string nprimariesString = infile->Get("NumberOfPrimaries")->GetTitle();
    totalPrimaries += std::stoi(nprimariesString,nullptr,10);

    infile->Close();
    delete infile;
    outFileMaGe->cd();
    
    TChain *fTree = new TChain("fTree");
    fTree->Add(dir+fileName+TString(".root"));
    
    Long64_t nentries = (Long64_t)fTree->GetEntries();
    
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
      /*
      Double_t x = primaries->GetX(),y = primaries->GetY(),z = primaries->GetZ();//,time = primaries->GetT();
      Double_t px = primaries->GetPx(),py = primaries->GetPy(),pz = primaries->GetPz();
      Double_t r = sqrt(x*x+y*y);
      Double_t theta = std::acos(x/r);
      bool hitSiPM = false,hitFiber = false;
      Int_t hitSiPMCounter = 0,hitFiberCounter = 0;
      if( y < 0) theta += 3.14159265359;
      */
      Int_t pastTrackID = -1;
      for (Int_t j = 0; j < eventSteps->GetNSteps();j++){
        step = eventSteps->GetStep(j);
        Int_t trackID = step->GetTrackID();
        physName = step->GetPhysVolName();
        Double_t x = step->GetLocalX(),y = step->GetLocalY(),z = step->GetLocalZ();//,time = step->GetT();
        Double_t r = sqrt(x*x+y*y);
        Double_t theta = std::acos(x/r);
        //cout<<"x "<<x<<", y "<<y<<" z "<<z<<" volName "<< physName<<", Edep "<<step->GetEdep()<<" process "<<step->GetTrackWeight()
        //  <<" trackID "<<trackID<<" pastTrackID "<<pastTrackID<<" parentID "<<step->GetParentTrackID()<<endl;
        if( y < 0) theta += 3.14159265359;
        //Track weight == 1 is for scintillation, track weight == 2 is for OpWLS
        if(physName.Contains("SiPM") && step->GetEdep() > 0 && step->GetTrackWeight() == 2 && pastTrackID != trackID){
          ntpleStepSiPM->Fill(i,step->GetEdep(),step->GetKineticE(),step->GetX(),step->GetY(),step->GetZ(),sqrt(step->GetX()*step->GetX() +step->GetY()*step->GetY()),x,y,z,r,theta);
          //cout<<"\tx "<<x<<", y "<<y<<" z "<<z<<" volName "<< physName<<endl;
          hMapUnscaled->Fill(x,y,z);
          h2DMapRZUnscaled->Fill(r,z);
          h2DMapXYUnscaled->Fill(x,y);
          h2DMapYZUnscaled->Fill(y,z);
        }
        //only fill when a new track is being looped over
        if(pastTrackID != trackID) pastTrackID = trackID;
      }
    }
    delete step;
    delete primaries;
    delete fTree;
  }
  //weight 3Dhistogram
  cout<<hMap->GetNbinsX()<<" "<<hMap->GetNbinsY()<<" "<<hMap->GetNbinsZ()<<endl;
  for(int i = 0; i < hMap->GetNbinsX();i++){
    for(int j = 0; j < hMap->GetNbinsY();j++){
      for(int k = 0; k < hMap->GetNbinsZ();k++){
        Double_t binVal = hMapUnscaled->GetBinContent(i,j,k);
        Double_t x = gridSpacing*(i)-maxX;
        Double_t y = gridSpacing*(j)-maxY;
        Double_t z = gridSpacing*(k)+minZ;
        if(binVal == 0){
            //binVal = hMapUnscaled->Interpolate(x,y,z);
            //cout<<"("<<x<<","<<y<<","<<z<<") ... bin center value"<<endl;
            n3D++;
        }
        hMap->SetBinContent(i,j,k,binVal/(totalPrimaries/(nbinsX*nbinsY*nbinsZ)));//totalEvents);
      }
    }
  }
  cout<<"finished XYZ Histogram"<<endl;
  for(int i = 0; i < h2DMapXY->GetNbinsX();i++){
    for(int j = 0; j < h2DMapXY->GetNbinsY();j++){
      Double_t binVal = h2DMapXYUnscaled->GetBinContent(i,j);
      Double_t x = gridSpacing*(i)-maxX;
      Double_t y = gridSpacing*(j)-maxY;
      if(binVal == 0){
         //binVal = h2DMapXYUnscaled->Interpolate(x,y);
         nXY++;
      }
      h2DMapXY->SetBinContent(i,j,binVal/(totalPrimaries/(nbinsX*nbinsY)));//totalEvents);
    }
  }
  cout<<"finished XY Histogram"<<endl;
  for(int i = 0; i < h2DMapYZ->GetNbinsX();i++){
    for(int j = 0; j < h2DMapYZ->GetNbinsY();j++){
      Double_t binVal = h2DMapYZUnscaled->GetBinContent(i,j);
      Double_t y = gridSpacing*(i)-maxY;
      Double_t z = gridSpacing*(j)+minZ;
      if(binVal == 0){
        //binVal = h2DMapYZUnscaled->Interpolate(y,z);
        nYZ++;
      }
      h2DMapYZ->SetBinContent(i,j,binVal/(totalPrimaries/(nbinsX*nbinsY)));//totalEvents);
    }
  }
  cout<<"finished YZ Histogram"<<endl;
  for(int i = 0; i < h2DMapRZ->GetNbinsX();i++){
    for(int j = 0; j < h2DMapRZ->GetNbinsY();j++){
      Double_t binVal = h2DMapRZUnscaled->GetBinContent(i,j);
      Double_t r = gridSpacing*(i)+minR;
      Double_t z = gridSpacing*(j)+minZ;
      if(binVal == 0){
        //binVal = h2DMapRZUnscaled->Interpolate(r,z);
        nRZ++; 
      }
    h2DMapRZ->SetBinContent(i,j,binVal/(totalPrimaries/(nbinsX*nbinsZ/2.)));//totalEvents);
    }
  }
  cout<<"finished RZ Histogram"<<endl;

  outFileMaGe->Write();
  outFileMaGe->Close();
  cout<<"Total events = "<<totalEvents<<", total primaries "<<totalPrimaries<<endl;
  cout<<"Events per bin..."<<totalPrimaries/(nbinsX*nbinsY)<<" XY, "<<totalPrimaries/(nbinsY*nbinsZ)
    <<" YZ, "<<totalPrimaries/(nbinsX*nbinsZ/2.)<<" RZ, "<<totalPrimaries/(nbinsX*nbinsY*nbinsZ)<<" XYZ"<<endl;
  cout<<"Histogram has "<<nbinsX*nbinsY<<" XY bins, "<<nbinsY*nbinsZ<<" YZ bins, "<<nbinsZ*nbinsX/2<<" RZ bins, "<<nbinsX*nbinsY*nbinsZ <<", xyz bins"<<endl;
  cout<<"Zero bins found..."<<nXY<<" XY, "<<nYZ<<" YZ, "<<nRZ<<" RZ bins, "<<n3D <<", xyz bins"<<endl;
  cout<<"root -l "<<fileName<<endl;
  return 0;
}
