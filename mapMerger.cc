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
  TString fileName = "OpticalMapL200.14String.5mm.2e9.root";

  TFile *outFileMaGe = new TFile(fileName,"recreate");
  
  //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
  TString dir = "/home/nmcfadden/RooT/mage/";
  
  //TString inFileNameTH3D = "OpticalMapLEGEND200.4e9_2mm.TwoFiber.5Rows_RawProb.root";
  //TString inFileNameTH3D = "OpticalMapLEGEND200.2.032e+09.19String.twoFiber.NO_HV.5mm.root";
  TString inFileNameTH3D = "OpticalMapLEGEND200.2e9.14String.root";
  //TString inFileNameTH3D = "OpticalMapLEGEND200.4e8.14String.25mm.root";
  //TString inFileNameTH2D = "OpticalMapLEGEND200.4e9_50mm.1e7.Exterior.root";
  TString inFileNameTH2D = "OpticalMapLEGEND200.Exterior.25mm.root";

  TH3D * h3D;
  TH2D* h2D;
  TFile* fileTH3D = TFile::Open(dir+inFileNameTH3D);
  fileTH3D->GetObject("EnergyMapUnscaled",h3D);
  h3D->SetNameTitle("ProbMapInterior","ProbMapInterior");
  TFile* fileTH2D = TFile::Open(dir+inFileNameTH2D);
  fileTH2D->GetObject("2DOpticalMap_RZ_Unscaled",h2D);
  h2D->SetNameTitle("ProbMapExterior","ProbMapExterior");
  outFileMaGe->cd();
  h3D->Write();
  h2D->Write();
  outFileMaGe->Write();
  outFileMaGe->Close();
  cout<<"root -l "<<fileName<<endl;
}

