// RootPulsesMap.C
// Author: Cloé Girard-Carillo <girardcarillo@lal.in2p3.fr>

// A macro to determine the signal velocity in coaxial cables, using three cables measured @LSM and
// using reflectometry
// Macro derived from RootPulsesMap.C

// root 'RootPulsesMap_LengthTests.C'

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

#include "/home/girardcarillo/Workdir/SNPlot/RootDisplay.h"

using namespace std;

//----------------------------//

/*
  RunData_02212019_17h40min31s_reflecto_slot7_cable_ext_Ascii --> 698 +/- 1 cm
  RunData_02212019_17h46min07s_reflecto_slot7_ext_et_court_Ascii --> X-1.0.0.0 345 +/- 1 cm
  RunData_02212019_17h51min01s_reflecto_slot7_ext_et_long_Ascii --> X-1.1.1.15 1327 +/- 1 cm
*/

const double c=2.997e8;// m/s
const double theoretical_celerity = 0.69*c;
const double tdc_to_ns = 0.390625;
const double adc_to_mv = 0.610351563;
const int channel_tot_number = 16;


void RootPulsesMap_LengthTests ()
{
  gROOT->Reset ();
  gStyle->SetOptStat (1);
  gStyle->SetPaintTextFormat("1.2f");
  gStyle->SetLineStyleString(11,"100 20");


  double tab_length[3];
  tab_length[0]=2*6.98;tab_length[1]=2*(6.98+3.45);tab_length[2]=2*(6.98+13.27);
  double tab_length_error[3] = {1,1,1};
  double tab_time[3] = {0,0,0};



  // Merging .root files
  for (int cable_counter = 0; cable_counter < 3; ++cable_counter) {
    int counts = 0 ;

    stringstream ss_cable;
    ss_cable << cable_counter;
    string str_cable = ss_cable.str();

    string filename_ = "root_files_LengthTests/Run105_LengthTests_"+str_cable+".root";

    TFile *theInFile = new TFile(filename_.c_str(),"READ") ;
    TTree *theTree = nullptr ;
    theInFile->GetObject("T",theTree) ;

    if (theInFile->IsOpen()) {
      cout << "File " << filename_ << " opened sucessfully" << endl ;
    }

    theTree = (TTree*)theInFile->Get("pulses") ;

    // cout << "The TTree has " << theTree->GetEntries () << " total entries."<< endl;

    // Declaration of leaf types:
    Double_t pulse_time_difference_CFD;

    // // Branches:
    theTree->SetBranchAddress ("pulse_time_difference_CFD", &pulse_time_difference_CFD);

    // Root histos:


    for (unsigned int i_pulse = 0; i_pulse < theTree->GetEntries(); i_pulse++){
      theTree->GetEntry (i_pulse);

      tab_time[cable_counter] += pulse_time_difference_CFD ;
      counts++ ;

    }

    tab_time[cable_counter] /= counts;
    tab_time[cable_counter] *= 1e-9;

  }

  TCanvas *canvas = new TCanvas ("canvas", "canvas", 10,10,2000,1000);

  TGraph *gr_celerity = new TGraph(3, tab_time, tab_length);

  config_graph(gr_celerity,"","t (s)","l (m)","AP",1.6,20,kRed+2) ;
  TF1 *fit = new TF1("linfit", "pol1");
  gr_celerity->Fit("linfit") ;
  gr_celerity->GetFunction("linfit")->SetLineColor(kOrange+7) ;
  gr_celerity->GetFunction("linfit")->SetLineWidth(1) ;

  TLegend* legend1 = new TLegend(0.178,0.755,0.515,0.990);

  float chi2lin = gr_celerity->GetFunction("linfit")->GetChisquare() ;
  float ndflin = fit->GetNDF() ;

  legend1->AddEntry(gr_celerity,"Measured cable","lep");
  legend1->AddEntry(gr_celerity->GetFunction("linfit"),"Linear fit","l");
  legend1->AddEntry((TObject*)0,Form("#chi^{2}/ndf  %1.4f/%1.f",chi2lin,ndflin),"");
  legend1->AddEntry((TObject*)0,Form("v_{p} /c    %1.3f #pm %1.4f",fit->GetParameter(1)/c,fit->GetParError(1)/c),"");

  legend1->Draw() ;

  canvas->SaveAs("plots/celerity.eps") ;

}//end RootPulsesMap_LengthTests
