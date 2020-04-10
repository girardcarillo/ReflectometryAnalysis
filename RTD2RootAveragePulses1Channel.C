// RTD2RootAveragePulses1Channel.C
// Author: Mathieu Bongrand <bongrand@lal.in2p3.fr>
// Copyright: 2019 (C) SuperNEMO - LAL (IN2P3/CNRS)

// root 'RTD2RootAveragePulses1Channel.C("rtd2root/snemo_run-103_slot-1_rtd.root")'

// root -b 'RTD2RootAveragePulses1Channel.C("rtd2root/snemo_run-103_slot-0_rtd.root")' -q

// root -b 'RTD2RootAveragePulses1Channel.C("rtd2root/snemo_run-105_wall-it_slot-0_rtd.root")' -q

#include "CnRoot.h"
#include "RTD2RootPulsesFunctions.h"

using namespace std;

//----------------------------//

void RTD2RootAveragePulses1Channel (std::string filename_, int channel_)
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
  TH1D *hChannel [nchannels];
  unsigned int channel_counts [nchannels];

  TRTD->GetEntry(0);
  int slot = calo_board_num [0];

  for (int i = 0; i < nchannels; i++)
    {
      std::string hname = "Board " + std::to_string (slot) + " - Channel " + std::to_string (i);
      hChannel [i] = new TH1D (hname.c_str (), hname.c_str (), 1024, 0, 399);
      hChannel [i]->GetXaxis ()->SetTitle ("Time [ns]");
      hChannel [i]->GetYaxis ()->SetTitle ("Amplitude [mV]");

      //hChannel [i]->SetLineWidth (2);
      hChannel [i]->SetLineColor (kAzure+i);

      channel_counts [i] = 0;

       for (size_t j = 0; j < MAX_WAVEFORM_SAMPLES; j++)
	{
	  hbins [i][j] = 0.0;
	}
    }

  for (size_t i = 0 ; i < ttree_entries; i ++)
    {
      TRTD->GetEntry(i);

      //if (i>100) break;

      double baseline_0 = compute_baseline (calo_ch0_waveform[0], 150);
      double baseline_1 = compute_baseline (calo_ch1_waveform[0], 150);

      int channel_0 = 2 * calo_chip_num [0];
      int channel_1 = 2 * calo_chip_num [0] + 1;

      channel_counts [channel_0]++;
      channel_counts [channel_1]++;

       for (size_t i = 0; i < MAX_WAVEFORM_SAMPLES; i++)
	{
	  hbins [channel_0][i] += calo_ch0_waveform[0][i] - baseline_0;
	  hbins [channel_1][i] += calo_ch1_waveform[0][i] - baseline_1;
	}

    }

  for (int i = 0; i < nchannels; i++)
    {
      for (size_t j = 0; j < MAX_WAVEFORM_SAMPLES; j++)
	{
	  double bin = ADC_TO_MV * hbins [i][j] / (double) channel_counts [i];
	  hbins [i][j] = bin;

	  hChannel[i]->SetBinContent (j+1, bin);
	}
    }

  TCanvas *canvas = new TCanvas ("canvas", "canvas", 100, 100, 1200, 400);
  canvas->cd ();
  hChannel [channel_]->Draw ();
  //hChannel [channel_]->GetXaxis ()->SetRangeUser (100,300);
  hChannel [channel_]->GetXaxis ()->SetRangeUser (200,400);

  ostringstream iss_canvas;
  iss_canvas << "plots/pdf/reflecto-it-board-" << slot << "-channel-" << channel_ << ".pdf";
  canvas->Print(iss_canvas.str().c_str ());
}
