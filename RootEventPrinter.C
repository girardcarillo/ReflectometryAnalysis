// RootEventPrinter.C
// Author: Mathieu Bongrand <bongrand@lal.in2p3.fr>
// Copyright: 2018 (C) SuperNEMO - LAL (IN2P3/CNRS)

// root 'RootEventPrinter.C("/sps/nemo/snemo/snemo_data/raw_data/CRD/run_104/RunData_12052018_12h58min53s_test2_Ascii.dat",0,0,1000)' -b -q

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"

using namespace std;

//----------------------------//

void RootEventPrinter (std::string filename_, int slot_, int channel_, int max_pulses_)
{
  const double tdc_to_ns = 0.390625;
  const double adc_to_mv = 0.610351563;

  // Read data file :

  ifstream file;

  file.open (filename_);

  if (!file.is_open ())
    {
      cerr << "ERROR: file isn't open !" << endl;
      return;
    }

  bool read_pulses = false;

  size_t pulse_read = 0;

  while (!file.eof ())
    {
      if (pulse_read % 10000 == 0) cout << "Read pulse: " << pulse_read << endl;

      if (pulse_read > max_pulses_) break;

      string file_line;
      getline (file, file_line);

      // read the event:

      string hit_str;
      int hit;
      size_t hit_found = file_line.find ("= HIT");
      if (hit_found != std::string::npos)
	{
	  istringstream hit_line (file_line);
	  hit_line >> hit_str >> hit_str >> hit;
	}

      ostringstream oss_slot_channel;
      oss_slot_channel << "Slot ";

      size_t event_found = file_line.find (oss_slot_channel.str ());
      if (event_found != std::string::npos)
	{
	  string slot_str, channel_str;
	  int slot, channel;

	  istringstream iss_line (file_line);
	  iss_line >> slot_str >> slot >> channel_str >> channel;

	  if (slot != slot_) continue;
	  if (channel != channel_) continue;

	  string pulse_line;
	  getline (file, pulse_line);

	  istringstream iss_pulse (pulse_line);

	  int isample = 0;
	  double baseline = 0.0;
	  double adc = 0.0;

	  vector<double> vec_adc;

	  while (iss_pulse >> adc >> ws)
	    {
	      vec_adc.push_back(adc);
	      if (isample < 150) baseline += adc;
	      isample++;
	    }

	  baseline /= 150.0;

	  TH1F *hpulse = new TH1F ("hpulse", "SuperNEMO run 104 - 5/12/2018; time [ns] ; amplitude [mV]", 1024, 0, 150);

	  size_t ibin = 1;
	  for (vector<double>::iterator it = vec_adc.begin(); it != vec_adc.end(); ++it)
	    {
	      double iamplitude = (*it - baseline)  * adc_to_mv;

	      hpulse->SetBinContent (ibin, iamplitude);
	      ibin++;
	    }

	  TCanvas *canvas = new TCanvas ("canvas", "canvas", 100, 100, 1300, 900);
	  canvas->cd ();
	  hpulse->Draw ();

	  ostringstream iss_canvas;
	  iss_canvas << "pulses/pulses-run104-" << hit << "-c" << slot_ <<"-r"<< channel_ << ".png";
	  canvas->Print(iss_canvas.str().c_str ());

	  }

      pulse_read++;
    }
  file.close();

}
