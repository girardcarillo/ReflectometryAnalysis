// RTD2RootTimeDelay1Channel.C
// Author: Mathieu Bongrand <bongrand@lal.in2p3.fr>
// Copyright: 2019 (C) SuperNEMO - LAL (IN2P3/CNRS)

// root 'RTD2RootTimeDelay1Channel.C("rtd2root/run_105/snemo_run-105_345cm_rtd.root", 14)'

// root -b 'RTD2RootTimeDelay1Channel.C("rtd2root/run_105/snemo_run-105_345cm_rtd.root", 14)' -q

#include "RTD2RootPulsesFunctions.h"

using namespace std;

//----------------------------//

void RTD2RootTimeDelay1Channel (std::string filename_, int channel_)
{
  // Root settings:
  gROOT->Reset ();
  gStyle->SetOptStat (0000);

  // Constants:
  static const uint16_t MAX_CALO_HITS    = 800;
  static const uint16_t MAX_WAVEFORM_SAMPLES = 1024;
  const double tdc_to_ns = 0.390625;
  const double adc_to_mv = 0.610351563;

  // TTrees:

  const char * MyTFile = filename_.c_str ();

  TChain * TRTD = new TChain("RTD");
  TRTD->Add(MyTFile);
  int ttree_entries = TRTD->GetEntries();

  cout << "Nb of events in tree: " << ttree_entries << endl;

  int16_t  calo_board_num[MAX_CALO_HITS];
  int16_t  calo_chip_num[MAX_CALO_HITS];
  TRTD->SetBranchAddress("calo_board_num", &calo_board_num);
  TRTD->SetBranchAddress("calo_chip_num", &calo_chip_num);

  int16_t  calo_ch0_waveform[MAX_CALO_HITS][MAX_WAVEFORM_SAMPLES];
  int16_t  calo_ch1_waveform[MAX_CALO_HITS][MAX_WAVEFORM_SAMPLES];
  TRTD->SetBranchAddress("calo_ch0_waveform", &calo_ch0_waveform);
  TRTD->SetBranchAddress("calo_ch1_waveform", &calo_ch1_waveform);


  const size_t nchannels = 14;
  double hbins [nchannels][MAX_WAVEFORM_SAMPLES];
  TH1D *hChannel = new TH1D ("Cable length 345 cm", "Cable length 345 cm", 1024, 0, 399);;
  hChannel->GetXaxis ()->SetTitle ("Time [ns]");
  hChannel->GetYaxis ()->SetTitle ("Amplitude [mV]");
  hChannel->SetLineColor (kAzure+1);

  TRTD->GetEntry(0);
  int slot = calo_board_num [0];

  for (size_t i = 0 ; i < ttree_entries; i ++)
    {
      TRTD->GetEntry (i);
      TRTD->Show (i);

      int channel_0 = 2 * calo_chip_num [0];
      int channel_1 = 2 * calo_chip_num [0] + 1;

      if (channel_0 != 14) continue;


      double baseline_0 = compute_baseline (calo_ch0_waveform[0], 150);
      double baseline_1 = compute_baseline (calo_ch1_waveform[0], 150);


       for (size_t i = 0; i < MAX_WAVEFORM_SAMPLES; i++)
	{
	  hChannel->SetBinContent (i+1, calo_ch0_waveform[0][i] - baseline_0);
	}

      if (channel_0 == 14)  break;
    }

  TCanvas *canvas = new TCanvas ("canvas", "canvas", 100, 100, 1200, 400);
  canvas->cd ();
  hChannel->Draw ();
  hChannel->GetXaxis ()->SetRangeUser (50,200);

}
