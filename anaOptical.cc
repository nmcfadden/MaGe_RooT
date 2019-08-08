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
  Double_t QE = 0.164;
  Double_t gridSpacing = 5;//*mm
  //Double_t maxX = 320.,maxY=320.0,maxZ=600,minZ = -200.,minR = 0,maxR=320;
  //Double_t maxX = 1000.,maxY=1000,maxZ=1200,minZ = -800.,minR = 320., maxR = 1000.;
  //Double_t maxX = 300.,minX = 0,maxY=300.,minY = 0,maxZ=925,minZ = -425.,minR = 0., maxR = 300.;
  Double_t maxX = 700.,minX = -700.,maxY=700.,minY = -700,maxZ= 845,minZ = -850,minR = 0, maxR = 700.;
  //Double_t maxX = 700.,minX = -700.,maxY=700.,minY = -700,maxZ= 850,minZ = -850,minR = 300, maxR = 700.;
  Int_t nbinsX = (maxX-minX)/gridSpacing,nbinsY = (maxY-minY)/gridSpacing,nbinsZ = (maxZ-minZ)/gridSpacing, nbinsR = (maxR-minR)/gridSpacing;
  TH3D* hMap = new TH3D("OpticalMap","OpticalMap",nbinsX,minX,maxX,nbinsY,minY,maxY,nbinsZ,minZ,maxZ);
  TH3D* hMapRaw = new TH3D("OpticalMap_Raw","OpticalMap_Raw",nbinsX,minX,maxX,nbinsY,minY,maxY,nbinsZ,minZ,maxZ);
  TH3D* hMapScaled = new TH3D("OpticalMap_Scaled","OpticalMap_Scaled",nbinsX,minX,maxX,nbinsY,minY,maxY,nbinsZ,minZ,maxZ);

  TH2D* h2DMapRZ = new TH2D("2DOpticalMap_RZ","2DOpticalMap_RZ",nbinsR,minR,maxR,nbinsZ,minZ,maxZ);
  TH2D* h2DMapRZRaw = new TH2D("2DOpticalMap_RZRaw","2DOpticalMap_RZRaw",nbinsR,minR,maxR,nbinsZ,minZ,maxZ);
  TH2D* h2DMapRZScaled = new TH2D("2DOpticalMap_RZScaled","2DOpticalMap_RZScaled",nbinsR,minR,maxR,nbinsZ,minZ,maxZ);

  TH2D* h2DMapXY = new TH2D("2DOpticalMap_XY","2DOpticalMap_XY",nbinsX,minX,maxX,nbinsY,minY,maxY);
  TH2D* h2DMapXYRaw = new TH2D("2DOpticalMap_XYRaw","2DOpticalMap_XYRaw",nbinsX,minX,maxX,nbinsY,minY,maxY);
  TH2D* h2DMapXYScaled = new TH2D("2DOpticalMap_XYScaled","2DOpticalMap_XYScaled",nbinsX,minX,maxX,nbinsY,minY,maxY);

  TH2D* h2DMapYZ = new TH2D("2DOpticalMap_YZ","2DOpticalMap_YZ",nbinsY,minY,maxY,nbinsZ,minZ,maxZ);
  TH2D* h2DMapYZRaw = new TH2D("2DOpticalMap_YZRaw","2DOpticalMap_YZRaw",nbinsY,minY,maxY,nbinsZ,minZ,maxZ);  
  TH2D* h2DMapYZScaled = new TH2D("2DOpticalMap_YZScaled","2DOpticalMap_YZScaled",nbinsY,minY,maxY,nbinsZ,minZ,maxZ);  
  
  Double_t totalEvents = 0,totalPrimaries;
  Int_t n3D = 0, nXY = 0, nYZ = 0, nRZ = 0;
  cout<<"Histogram has "<<nbinsX<<" "<<nbinsY<<" "<<nbinsZ<<" "<<nbinsX*nbinsY*nbinsZ <<" bins"<<endl;

  cout<<"starting run"<<endl; 
  TFile * infile;
  for(int k =0; k <= 339; k++){
    //TODO--change input dir
    TString dir = "";
    //dir = "/data4/legend/sim/l200sims/LGND200_14String_Optical/default/";
    dir = "/data4/legend/sim/l200sims/LGND200_14String_Optical/realistic/";
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
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/14StringBaseline/";
    //TString dir = "/home/nmcfadden/MaGe/bin/Linux-g++/";
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
    TString fileName = "14StringBaselineInterior"+to_string(k);
    //TString fileName = "14StringBaselineExterior"+to_string(k);
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
      cout<<"NumberOfEvents is NULL... maybe file has not closed yet"<<endl;
      continue;
    }
    string neventsString = infile->Get("NumberOfEvents")->GetTitle();
    totalEvents += std::stoi(neventsString,nullptr,10);

    string nprimariesString = infile->Get("NumberOfPrimaries")->GetTitle();
    if(fileName.Contains("Exterior")) totalPrimaries += 78400;//std::stoi(nprimariesString,nullptr,10);
    else totalPrimaries += std::stoi(nprimariesString,nullptr,10);
    cout<<"adding to total Primaries "<<std::stoi(nprimariesString,nullptr,10)<<"... Total Primaries is "<<totalPrimaries<<endl;
    infile->Close();
    delete infile;
    outFileMaGe->cd();
    
    TChain *fTree = new TChain("fTree");
    fTree->Add(dir+fileName+TString(".root"));
    if(fTree == NULL) break; 
    Long64_t nentries = (Long64_t)fTree->GetEntries();
    
    cout<<". Processing "<<nentries<<" entries"<<endl;
    MGTMCEventSteps *eventSteps = 0;
    MGTMCEventSteps *eventPrimaries = 0;
    
    if(fTree->GetBranchStatus("eventSteps")){
      fTree->SetBranchAddress("eventSteps",&eventSteps);
    }
    else{
      cout<<"eventSteps is NULL"<<endl;
      continue;
    }

    if(fTree->GetBranchStatus("eventPrimaries")){
      fTree->SetBranchAddress("eventPrimaries",&eventPrimaries);
    }
    else{
      cout<<"eventSteps is NULL"<<endl;
      continue;
    }  
    
    const MGTMCStepData *step,*primaries;
    
    for(Int_t i = 0; i < nentries ; i++){
      if((i+1)%10000 == 0 || i == nentries - 1 ) cout<<"\tprocessed "<<i+1<<" events"<<endl;
      //if(i==0) break;
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
          //cout<<"\tx "<<x<<", y "<<y<<" z "<<z<<" r "<<r<<" volName "<< physName<<endl;
          hMapRaw->Fill(x,y,z);
          h2DMapRZRaw->Fill(r,z);
          h2DMapXYRaw->Fill(x,y);
          h2DMapYZRaw->Fill(y,z);
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
  cout<<hMap->GetNbinsX()<<" "<<hMap->GetNbinsY()<<" "<<hMap->GetNbinsZ()<<" rBins "<<h2DMapRZRaw->GetNbinsX()<<endl;
  for(int i = 0; i < hMap->GetNbinsX();i++){
    for(int j = 0; j < hMap->GetNbinsY();j++){
      for(int k = 0; k < hMap->GetNbinsZ();k++){
        Double_t binVal = hMapRaw->GetBinContent(i+1,j+1,k+1);
        Double_t x = gridSpacing*(i)-maxX;
        Double_t y = gridSpacing*(j)-maxY;
        Double_t z = gridSpacing*(k)+minZ;
        if(binVal == 0){
            //binVal = hMapRaw->Interpolate(x,y,z);
            //cout<<"("<<x<<","<<y<<","<<z<<") ... bin center value"<<endl;
            n3D++;
        }
        //hMap->SetBinContent(i+1,j+1,k+1,binVal/(totalPrimaries/(nbinsX*nbinsY*nbinsZ)));//totalEvents);
        hMap->SetBinContent(i+1,j+1,k+1,binVal/(500));//totalPrimaries/(nbinsX*nbinsY*nbinsZ)));//totalEvents);
	hMapScaled->SetBinContent(i+1,j+1,k+1,QE*binVal/(500));
      }
    }
  }
  cout<<"finished XYZ Histogram"<<endl;
  for(int i = 0; i < h2DMapXY->GetNbinsX();i++){
    for(int j = 0; j < h2DMapXY->GetNbinsY();j++){
      Double_t binVal = h2DMapXYRaw->GetBinContent(i+1,j+1);
      Double_t x = gridSpacing*(i)-maxX;
      Double_t y = gridSpacing*(j)-maxY;
      if(binVal == 0){
         //binVal = h2DMapXYRaw->Interpolate(x,y);
         nXY++;
      }
      h2DMapXY->SetBinContent(i+1,j+1,binVal/(totalPrimaries/(nbinsX*nbinsY)));//totalEvents);
      h2DMapXYScaled->SetBinContent(i+1,j+1,QE*binVal/(totalPrimaries/(nbinsX*nbinsY)));//totalEvents);
    }
  }
  cout<<"finished XY Histogram"<<endl;
  for(int i = 0; i < h2DMapYZ->GetNbinsX();i++){
    for(int j = 0; j < h2DMapYZ->GetNbinsY();j++){
      Double_t binVal = h2DMapYZRaw->GetBinContent(i+1,j+1);
      Double_t y = gridSpacing*(i+1)-maxY;
      Double_t z = gridSpacing*(j+1)+minZ;
      if(binVal == 0){
        //binVal = h2DMapYZRaw->Interpolate(y,z);
        nYZ++;
      }
      h2DMapYZ->SetBinContent(i+1,j+1,binVal/(totalPrimaries/(nbinsX*nbinsY)));//totalEvents);
      h2DMapYZScaled->SetBinContent(i+1,j+1,QE*binVal/(totalPrimaries/(nbinsX*nbinsY)));//totalEvents);
    }
  }
  cout<<"finished YZ Histogram"<<endl;
TH1D * hWeight = new TH1D("rWeight","rWeight",h2DMapRZ->GetNbinsX(),minR,maxR);
  for(int i = 0; i < h2DMapRZ->GetNbinsX();i++){
    for(int j = 0; j < h2DMapRZ->GetNbinsY();j++){
      Double_t binVal = h2DMapRZRaw->GetBinContent(i+1,j+1);
      Double_t r = gridSpacing*(i+1)+minR;
      Double_t z = gridSpacing*(j+1)+minZ;
      Double_t weight = ((r+gridSpacing)*(r+gridSpacing)-r*r)/(maxR*maxR-minR*minR);
      binVal /= weight;
      if(binVal == 0){
        //binVal = h2DMapRZRaw->Interpolate(r,z);
        nRZ++; 
      }
      hWeight->SetBinContent(i+1,weight);
    //h2DMapRZ->SetBinContent(i+1,j+1,binVal/(totalPrimaries/(nbinsR*nbinsZ)));//totalEvents);
    h2DMapRZ->SetBinContent(i+1,j+1,binVal/(totalPrimaries/(nbinsZ)));//totalEvents);
    h2DMapRZScaled->SetBinContent(i+1,j+1,QE*binVal/(totalPrimaries/(nbinsZ)));//totalEvents);
    }
  }
  cout<<"finished RZ Histogram"<<endl;

  outFileMaGe->Write();
  outFileMaGe->Close();
  cout<<"Total events = "<<totalEvents<<", total primaries "<<totalPrimaries<<endl;
  cout<<"Events per bin..."<<totalPrimaries/(nbinsX*nbinsY)<<" XY, "<<totalPrimaries/(nbinsY*nbinsZ)
    <<" YZ, "<<totalPrimaries/(nbinsR*nbinsZ)<<" RZ, "<<totalPrimaries/(nbinsX*nbinsY*nbinsZ)<<" XYZ"<<endl;
  cout<<"Histogram has "<<nbinsX*nbinsY<<" XY bins, "<<nbinsY*nbinsZ<<" YZ bins, "<<nbinsZ*nbinsR<<" RZ bins, "<<nbinsX*nbinsY*nbinsZ <<", xyz bins"<<endl;
  cout<<"Zero bins found..."<<nXY<<" XY, "<<nYZ<<" YZ, "<<nRZ<<" RZ bins, "<<n3D <<", xyz bins"<<endl;
  cout<<"root -l "<<fileName<<endl;
  return 0;
}
