// RootPulsesAnalysis.C
// Author: Mathieu Bongrand <bongrand@lal.in2p3.fr>
// Copyright: 2018 (C) SuperNEMO - LAL (IN2P3/CNRS)

// root 'RootPulsesAnalysis.C("data/test1/root_pulses_test1_run0.root", 1 , 0)'
// root 'RootPulsesAnalysis.C("data/root_run100_0_9.root", 1 , 0)'
// root 'RootPulsesAnalysis.C("data/root_run104_0_64.root", 6, 12)'

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include <TROOT.h>
#include <TStyle.h>

#include <TFile.h>
#include <TTree.h>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>

using namespace std;

//----------------------------//



void RootPulsesAnalysis (string filename_,int slot_, int channel_)
{
  gROOT->Reset ();
  // gStyle->SetHistFillColor (kAzure-8);
  gStyle->SetOptStat (0000);

  TH2D *h_counts = new TH2D ("counts","counts; Y; Z", 22, -1, 21, 16, -1, 14);

  // Output tree:
  TFile *rootfile = new TFile (filename_.c_str (), "read");
  TTree *TPulses = (TTree*) rootfile->Get("pulses");

  cout << TPulses->GetEntries () << endl;

  // Declaration of leaf types
  Int_t pulse_slot;
  Int_t pulse_channel;
  Double_t pulse_time;
  Double_t pulse_amplitude;
  Double_t pulse_charge;
  Double_t pulse_baseline;

  // Branches:
  TPulses->SetBranchAddress ("pulse_slot", &pulse_slot);
  TPulses->SetBranchAddress ("pulse_channel", &pulse_channel);
  TPulses->SetBranchAddress ("pulse_time", &pulse_time);
  TPulses->SetBranchAddress ("pulse_amplitude", &pulse_amplitude);
  TPulses->SetBranchAddress ("pulse_charge", &pulse_charge);
  TPulses->SetBranchAddress ("pulse_baseline", &pulse_baseline);

  // Root histos:
  ostringstream iss;
  iss << "OM: C" << slot_ <<"-R"<< channel_ << "; charge [V.nS]";

  TH1D *h_amplitude = new TH1D ("amplitude","; amplitude [mV]", 40, 0, 800);
  TH1D *h_charge = new TH1D ("charge",iss.str().c_str (), 50, 0, 100);
  TH2D *h_charge_amplitude = new TH2D ("charge vs amplitude","; amplitude [mV]; charge [V.nS]", 80, 0, 800, 50, 0, 100);

  for ( unsigned int i_pulse = 0; i_pulse < TPulses->GetEntries(); i_pulse++ )
    {
      TPulses->GetEntry (i_pulse);

      if (pulse_slot != slot_) continue;
      if (pulse_channel != channel_) continue;

      //if (pulse_amplitude > -100.0) continue;

      //if (i_pulse < 1000) TPulses->Show (i_pulse);

      double amplitude = fabs (pulse_amplitude);
      double charge = fabs (pulse_charge);

      h_amplitude->Fill (amplitude);
      h_charge->Fill (charge);
      h_charge_amplitude->Fill (amplitude,charge);
    }

  TCanvas *canvas = new TCanvas ("canvas", "canvas", 200, 200, 1400, 700);
  canvas->Divide(3,1);
  canvas->cd (1);
  gPad->SetLogy ();
  h_amplitude->Draw ();
  canvas->cd (2);
  gPad->SetLogy ();
  h_charge->Draw ();
  canvas->cd (3);
  h_charge_amplitude->Draw ("colz");

  ostringstream iss_canvas;
  iss_canvas << "plot-run104-c" << slot_ <<"-r"<< channel_ << ".pdf";
  canvas->Print(iss_canvas.str().c_str ());
}
