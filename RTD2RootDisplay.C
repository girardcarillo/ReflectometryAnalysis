// RTD2RootDisplay.C
// Author: Mathieu Bongrand <bongrand@lal.in2p3.fr>
// Copyright: 2019 (C) SuperNEMO - LAL (IN2P3/CNRS)

// root 'RTD2RootDisplay.C("/sps/nemo/snemo/snemo_data/raw_data/RTD/root/run_104/snemo_run-104_rtd.root",1)'
// root 'RTD2RootDisplay.C("/sps/nemo/snemo/snemo_data/raw_data/RTD/root/run_103/RunData_12042018_15h51min12s_Slot1_Ascii/snemo_run-104_rtd.root",1)'
// root 'RTD2RootDisplay.C("rtd2root/snemo_run-103_slot-1_rtd.root",0)'

#include "CnRoot.h"
#include "RTD2RootPulsesFunctions.h"

using namespace std;

void RTD2RootDisplay (std::string filename_, int entry_)
{
  // Settings:
  bool zoom = false;

  // Root settings:
  gROOT->Reset ();
  //gStyle->SetHistFillColor (kAzure-8);
  gStyle->SetOptStat (00000);
  //gStyle->SetOptFit (1111);
  // gStyle->SetStatX (0.89); // 48
  // gStyle->SetStatY (0.89);
  gStyle->SetPadGridX (true);
  gStyle->SetPadGridY (true);

  // Constants:
  static const uint16_t MAX_CALO_HITS    = 800;
  static const uint16_t MAX_WAVEFORM_SAMPLES = 1024;
  const double tdc_to_ns = 0.390625;
  const double adc_to_mv = 0.610351563;

  // TTrees:

  const char * MyTFile = filename_.c_str ();

  TChain * TRTD = new TChain("RTD");
  TRTD->Add(MyTFile);
  TRTD->Show(entry_);

  int16_t calo_ch0_baseline[MAX_CALO_HITS];
  int16_t calo_ch1_baseline[MAX_CALO_HITS];
  TRTD->SetBranchAddress("calo_ch0_baseline", &calo_ch0_baseline);
  TRTD->SetBranchAddress("calo_ch1_baseline", &calo_ch1_baseline);

  int16_t  calo_ch0_waveform[MAX_CALO_HITS][MAX_WAVEFORM_SAMPLES];
  int16_t  calo_ch1_waveform[MAX_CALO_HITS][MAX_WAVEFORM_SAMPLES];
  TRTD->SetBranchAddress("calo_ch0_waveform", &calo_ch0_waveform);
  TRTD->SetBranchAddress("calo_ch1_waveform", &calo_ch1_waveform);

  int16_t  calo_board_num[MAX_CALO_HITS];
  TRTD->SetBranchAddress("calo_board_num", &calo_board_num);
  int16_t  calo_chip_num[MAX_CALO_HITS];
  TRTD->SetBranchAddress("calo_chip_num", &calo_chip_num);

  bool calo_ch0_lt;
  TRTD->SetBranchAddress("calo_ch0_lt", &calo_ch0_lt);
  bool calo_ch1_lt;
  TRTD->SetBranchAddress("calo_ch1_lt", &calo_ch1_lt);

  bool calo_ch0_ht;
  TRTD->SetBranchAddress("calo_ch0_ht", &calo_ch0_ht);
  bool calo_ch1_ht;
  TRTD->SetBranchAddress("calo_ch1_ht", &calo_ch1_ht);

  for (int i = entry_ ; i < TRTD->GetEntries(); i++)
    {
      if (i % 10000 == 0) cout << "Read : " << i << endl;

      TRTD->GetEntry(i);
      if (calo_board_num[0] == 14 && calo_chip_num[0] == 4 && calo_ch1_lt == 1 && calo_ch1_ht == 0)
	{
	  TRTD->Show(i);
	  break;
	}
    }

  TRTD->GetEntry(entry_);

  TH1F *hpulse_0 = new TH1F ("hpulse_0", "; time [ns] ; amplitude [mV]", 1024, 0, 399);
  TH1F *hderiv_0 = new TH1F ("hderiv_0", "; time [ns] ; amplitude [mV]", 1024, 0, 399);
  TH1F *hpulse_1 = new TH1F ("hpulse_1", "; time [ns] ; amplitude [mV]", 1024, 0, 399);
  TH1F *hderiv_1 = new TH1F ("hderiv_1", "; time [ns] ; amplitude [mV]", 1024, 0, 399);

  // compute parameters:

  double baseline_0 = compute_baseline (calo_ch0_waveform[0], 150);
  double baseline_1 = compute_baseline (calo_ch1_waveform[0], 150);

  cout << "----------------" << endl;
  cout << "Baselines: " << baseline_0 <<' '<< baseline_1 << endl;

  double amplitude_0 = compute_amplitude (calo_ch0_waveform[0], 250, 600, +1.0, baseline_0);
  double amplitude_1 = compute_amplitude (calo_ch1_waveform[0], 250, 600, +1.0, baseline_1);

  cout << "Amplitudes: " << amplitude_0 <<' '<< amplitude_1 << endl;

  double cfd_time_0 = compute_cfd_time (calo_ch0_waveform[0], 250, 600, +1.0, baseline_0, amplitude_0);
  double cfd_time_1 = compute_cfd_time (calo_ch1_waveform[0], 250, 600, +1.0, baseline_1, amplitude_1);

  cout << "Cfd_Times: " << cfd_time_0 <<' '<< cfd_time_1 << endl;
  cout << "Debug: " << compute_cfd_time (calo_ch0_waveform[0], 650, 1000, +1.0, baseline_0, amplitude_0) << endl;

  double charge_0 = 0.001 * compute_charge (calo_ch0_waveform[0], 250, 600, baseline_0);
  double charge_1 = 0.001 * compute_charge (calo_ch1_waveform[0], 250, 600, baseline_1);

  cout << "Charges: " << charge_0 <<' '<< charge_1 << endl;
  cout << "----------------" << endl;

  // fill pulse histos:

  double previous_ampl_0 = 0.0;
  double previous_ampl_1 = 0.0;

  for (size_t i = 0; i < MAX_WAVEFORM_SAMPLES; i++)
    {
      double iamplitude_0 = (calo_ch0_waveform[0][i] - baseline_0) * ADC_TO_MV;
      double iamplitude_1 = (calo_ch1_waveform[0][i] - baseline_1) * ADC_TO_MV;

      hpulse_0->SetBinContent (i+1, iamplitude_0);
      hpulse_1->SetBinContent (i+1, iamplitude_1);

      if (i>0)
	{
	  hderiv_0->SetBinContent (i+1, (iamplitude_0 - previous_ampl_0) / TDC_TO_NS);
	  hderiv_1->SetBinContent (i+1, (iamplitude_1 - previous_ampl_1) / TDC_TO_NS);
	}
      else
	{
	  hderiv_0->SetBinContent (i+1, 0);
	  hderiv_1->SetBinContent (i+1, 0);
	}

      previous_ampl_0 = iamplitude_0;
      previous_ampl_1 = iamplitude_1;
    }

  TCanvas *c0 = new TCanvas ("c0", "c0", 100, 100, 1300, 900);
  c0->Divide (1,2);
  c0->cd (1);
  hpulse_0->Draw ();
  if (zoom) hpulse_0->GetXaxis ()->SetRangeUser (150,300);

  c0->cd (2);
  hderiv_0->Draw ();
  if (zoom) hderiv_0->GetXaxis ()->SetRangeUser (150,300);

  TCanvas *c1 = new TCanvas ("c1", "c1", 150, 150, 1300, 900);
  c1->Divide (1,2);
  c1->cd (1);
  hpulse_1->Draw ();
  if (zoom) hpulse_1->GetXaxis ()->SetRangeUser (150,300);

  c1->cd (2);
  hderiv_1->Draw ();
  if (zoom) hderiv_1->GetXaxis ()->SetRangeUser (150,300);

}
