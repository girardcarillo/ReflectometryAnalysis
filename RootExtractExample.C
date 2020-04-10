// RootPulsesMap.C
// Author: Cloé Girard-Carillo <girardcarillo@lal.in2p3.fr>

// .x RootPulsesMap.C+("103","fr","back")

#include <iostream>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TLegend.h>

#include "RootDisplay.h"

using namespace std;

//----------------------------//


double c=2.99e1;// cm/ps
double theoretical_celerity = 0.69*c;

void RootPulsesMap_LengthTests ()
{
  gROOT->Reset ();
  gStyle->SetOptStat (0);
  gStyle->SetPaintTextFormat("1.2f");
  gStyle->SetLineStyleString(11,"100 20");

  TCanvas *canvas = new TCanvas ("canvas", "canvas", 200, 200, 1100, 800);

  // Root histos:




  // Merging Root files:

  TList *list = new TList;

  for (int cable_counter = 0; cable_counter < 3; ++cable_counter) {
    stringstream ss_cable;
    ss_cable << cable_counter;
    string str_cable = ss_cable.str();

    string filename_ = "root_files_LengthTests/Run105_LengthTests_"+str_cable+".root";
    const char* c_filename = filename_.c_str();

    TFile *file_adress = TFile::Open(c_filename);
    TTree *tree_adress = (TTree*)file_adress->Get("pulses");
    list->Add(tree_adress);
  }

  TFile* outputfile = TFile::Open("output.root", "recreate");
  TTree *TPulses = TTree::MergeTrees(list);

  cout << "The TTree has " << TPulses->GetEntries () << " total entries."<< endl;

  // Declaration of leaf types:

  Int_t pulse_slot;
  Int_t pulse_channel;
  Double_t pulse_amplitude_first_peak;
  Double_t pulse_charge_first_peak;
  Double_t pulse_charge_second_peak;
  Double_t pulse_baseline;
  Double_t pulse_amplitude_second_peak;
  Double_t pulse_attenuation;
  Double_t pulse_time_difference_CFD;
  Double_t pulse_time_difference_max;
  Int_t pulse_time_second_pic;

  // Branches:

  TPulses->SetBranchAddress ("pulse_slot", &pulse_slot);
  TPulses->SetBranchAddress ("pulse_channel", &pulse_channel);
  TPulses->SetBranchAddress ("pulse_amplitude_first_peak", &pulse_amplitude_first_peak);
  TPulses->SetBranchAddress ("pulse_amplitude_second_peak", &pulse_amplitude_second_peak);
  TPulses->SetBranchAddress ("pulse_charge_first_peak", &pulse_charge_first_peak);
  TPulses->SetBranchAddress ("pulse_charge_second_peak", &pulse_charge_second_peak);
  TPulses->SetBranchAddress ("pulse_baseline", &pulse_baseline);
  TPulses->SetBranchAddress ("pulse_attenuation", &pulse_attenuation);
  TPulses->SetBranchAddress ("pulse_time_difference_CFD", &pulse_time_difference_CFD);
  TPulses->SetBranchAddress ("pulse_time_difference_max", &pulse_time_difference_max);
  TPulses->SetBranchAddress ("pulse_time_second_pic", &pulse_time_second_pic);

  for (unsigned int i_pulse = 0; i_pulse < TPulses->GetEntries(); i_pulse++){
    TPulses->GetEntry (i_pulse);

  }

  TPulses->Write();
  outputfile->Close();
  delete outputfile;


}//end RootPulsesMap_LengthTests
