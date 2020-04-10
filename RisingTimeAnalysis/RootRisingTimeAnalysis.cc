// RootPulsesMap.C
// Author: Cloé Girard-Carillo <girardcarillo@lal.in2p3.fr>

// This macro is a simplified version of RootPulsesMap.C (created on 19/02/2020 for editing results for thesis manuscript)
// .x RootTimingAnalysis.cc+("103","fr","back")

#include <iostream>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include <TROOT.h>
#include <TLine.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TMathText.h>

#include "/home/girardcarillo/Workdir/SNPlot/RootDisplay.h"

using namespace std;

//----------------------------//

double c=2.99e1 ;// cm/ns
const double amp_bin = 0.6 ; //mV
const int channel_tot_number = 13 ;
const double measured_celerity = 0.7*c ;

void RootRisingTimeAnalysis (bool enable_drawing = 0) {

  // Root histos:


  TH1F *hrisingTime = new TH1F ("RisingTimes","", 100, 0, 0) ;
  TProfile *prisingtime = new TProfile ("prisingtime","", 10, 9, 14, 2, 20) ;


  // Merging .root files
  TList *list = new TList ;

  string filename1_ = "root_files/RisingTime_slot0.root" ;
  string filename2_ = "root_files/RisingTime_slot9.root" ;
  string filename3_ = "root_files/RisingTime_slot20.root" ;

  TFile *file_adress1 = TFile::Open(filename1_.c_str()) ;
  TTree *tree_adress1 = (TTree*)file_adress1->Get("pulses") ;
  list->Add(tree_adress1) ;
  TFile *file_adress2 = TFile::Open(filename2_.c_str()) ;
  TTree *tree_adress2 = (TTree*)file_adress2->Get("pulses") ;
  list->Add(tree_adress2) ;
  TFile *file_adress3 = TFile::Open(filename3_.c_str()) ;
  TTree *tree_adress3 = (TTree*)file_adress3->Get("pulses") ;
  list->Add(tree_adress3) ;

  TFile* outputfile = TFile::Open("output.root", "recreate") ;
  TTree *TPulses = TTree::MergeTrees(list) ;

  cout << "The TTree has " << TPulses->GetEntries () << " total entries."<< endl ;

  // Declaration of leaf types
  Int_t pulse_slot ;
  Int_t pulse_channel ;
  Double_t pulse_rising_time ;
  Double_t pulse_amplitude_first_peak ;
  Double_t pulse_charge_first_peak ;
  Double_t pulse_charge_second_peak ;
  Double_t pulse_baseline ;
  Double_t pulse_amplitude_second_peak ;
  Double_t pulse_attenuation ;
  Double_t pulse_time_difference_CFD ;
  Double_t pulse_time_difference_max ;
  Int_t pulse_time_second_pic ;

  // Branches:
  TPulses->SetBranchAddress ("pulse_slot", &pulse_slot) ;
  TPulses->SetBranchAddress ("pulse_channel", &pulse_channel) ;
  TPulses->SetBranchAddress ("pulse_rising_time", &pulse_rising_time) ;
  TPulses->SetBranchAddress ("pulse_amplitude_first_peak", &pulse_amplitude_first_peak) ;
  TPulses->SetBranchAddress ("pulse_amplitude_second_peak", &pulse_amplitude_second_peak) ;
  TPulses->SetBranchAddress ("pulse_charge_first_peak", &pulse_charge_first_peak) ;
  TPulses->SetBranchAddress ("pulse_charge_second_peak", &pulse_charge_second_peak) ;
  TPulses->SetBranchAddress ("pulse_baseline", &pulse_baseline) ;
  TPulses->SetBranchAddress ("pulse_attenuation", &pulse_attenuation) ;
  TPulses->SetBranchAddress ("pulse_time_difference_CFD", &pulse_time_difference_CFD) ;
  TPulses->SetBranchAddress ("pulse_time_difference_max", &pulse_time_difference_max) ;
  TPulses->SetBranchAddress ("pulse_time_second_pic", &pulse_time_second_pic) ;

  vector < vector <double> > risingtime(channel_tot_number) ;
  vector < vector <double> > time_differences(channel_tot_number) ;
  vector < double > mean_time_difference(channel_tot_number) ;

  for (unsigned int i_pulse = 0 ; i_pulse < TPulses->GetEntries() ; i_pulse++){
    TPulses->GetEntry (i_pulse) ;

    if (pulse_channel <= channel_tot_number) {

      if (pulse_slot == 0) {

        time_differences[pulse_channel].push_back(pulse_time_difference_CFD) ;
        risingtime[pulse_channel].push_back(pulse_rising_time) ;

      }
    }
  }

  double length[channel_tot_number] ;

  for (int i = 0 ; i < channel_tot_number ; ++i) {

    mean_time_difference[i] = TMath::Mean(time_differences[i].begin(),time_differences[i].end()) ;
    length[i] = (measured_celerity*mean_time_difference[i])/2./100. ; // Convert cm in m

    for (int k = 0 ; k < risingtime[i].size() ; ++k) {

      prisingtime->Fill(length[i],risingtime[i].at(k)) ;

    }

  }

  if (enable_drawing) {
    gStyle->SetOptStat(0) ;
    config_profile(prisingtime, "", "l (m)", "Rising time (ns)","",2,kCyan+2) ;
  }


}//end macro
