// DisplayPulses.C
// Author: Cloe Girard-Carillo <girardcarillo@lal.in2p3.fr>
// Copyright: 2018 (C) SuperNEMO - LAL (IN2P3/CNRS)

// example :
// .x DisplayPulse.C("/sps/nemo/snemo/snemo_data/raw_data/CRD/run_103/RunData_12042018_15h51min12s_Slot1_Ascii/RunData_12042018_15h51min12s_Slot1_Ascii.dat")

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLine.h"
#include "TROOT.h"
#include "TGraph.h"

#include "/home/girardcarillo/Workdir/SNPlot/RootDisplay.h"

using namespace std;

TH1F *pulse(string filename, int pulse_number);
void config_histo(TH1F *histo,int linewidth,Color_t color,const char *draw_type);

void DisplayPulse (string filename)
{

  TH1F* hpulse = pulse(filename,1);
  // TH1F* hpulse1 = pulse(filename,28);
  // TH1F* hpulse2 = pulse(filename,44);
  // TH1F* hpulse3 = pulse(filename,62);
  // TH1F* hpulse4 = pulse(filename,74);
  // TH1F* hpulse5 = pulse(filename,90);
  // TH1F* hpulse6 = pulse(filename,102);

  TCanvas *c = new TCanvas();
  config_histo1D(hpulse,"","Time (ns)","Amplitude (mV)",1,1,1) ;
  // config_histo(hpulse1,2,2,"SAME");
  // config_histo(hpulse2,2,3,"SAME");
  // config_histo(hpulse3,2,4,"SAME");
  // config_histo(hpulse4,2,5,"SAME");
  // config_histo(hpulse5,2,6,"SAME");
  // config_histo(hpulse6,2,7,"SAME");

  // auto legend = new TLegend(0.2,0.7,0.55,0.9);
  // legend->AddEntry(hpulse,hpulse->GetName(),"l");
  // legend->AddEntry(hpulse1,hpulse1->GetName(),"l");
  // legend->AddEntry(hpulse2,hpulse2->GetName(),"l");
  // legend->AddEntry(hpulse3,hpulse3->GetName(),"l");
  // legend->AddEntry(hpulse4,hpulse4->GetName(),"l");
  // legend->AddEntry(hpulse5,hpulse5->GetName(),"l");
  // legend->AddEntry(hpulse6,hpulse6->GetName(),"l");

  // legend->Draw();

  c->SaveAs("plots/pulse.pdf");

}


TH1F *pulse(string filename, int pulse_number){

  const double tdc_to_ns = 0.390625;
  const double adc_to_mv = 0.610351563;

  // Prepare histos :

  gROOT->Reset ();
  gStyle->SetOptStat (0000);

  stringstream ss;
  ss << pulse_number;
  string s_histo_name = string("hit")+ss.str();
  const char *histo_name = s_histo_name.c_str();

  TH1F *hpulse = new TH1F (histo_name, "", 1024, 0, 400);

  // Read data file :

  ifstream file;

  file.open (filename);

  if (!file.is_open ())
    {
      cerr << "ERROR: file isn't open !" << endl;
    }

  vector<double> vec_adc;

  bool read_pulses = false;

  ostringstream oss0;
  oss0 << pulse_number;
  string str_pulse_nb = "= HIT " + oss0.str ();

  while (read_pulses == false)
    {
      string hit_line, slot_line;
      getline(file , hit_line);

      size_t found1 = hit_line.find(str_pulse_nb);
      if (found1 != string::npos)
	{
	  read_pulses = true;
	}
      getline(file , slot_line);

      size_t found2 = hit_line.find("Slot ");
      if (found2 != string::npos)
	{

	  string tmp1, tmp2;
	  int slot, channel;
	  istringstream iss(slot_line);
	  iss >> tmp1 >> slot >> tmp2 >> channel;

          cout << "slot " << found2 << endl ;

	}
      size_t found3 = hit_line.find("Ch ");
      if (found3 != string::npos)
	{
          cout << "channel " << found3 << endl ;
        }
    }

  // Read the pulse

  string line_to_read;

  getline(file , line_to_read);
  istringstream iss(line_to_read);

  int nb_samples = 0;
  double time = 0.0, charge = 0.0;
  double amplitude = 0.0, baseline = 0.0;
  double adc = 0.0;

  while (iss >> adc >> ws)
    {
      vec_adc.push_back(adc);
      if (nb_samples < 150) baseline += adc;
      nb_samples++;
    }
  baseline /= 150.0;

  size_t ibin = 1;
  for (vector<double>::iterator it = vec_adc.begin(); it != vec_adc.end(); ++it)
    {
      double iamplitude = (*it - baseline)  * adc_to_mv;

      if (iamplitude < amplitude)
	{
	  amplitude = iamplitude;
	  time = ibin * tdc_to_ns;
	}

      if (ibin > 150) charge += iamplitude;

      hpulse->SetBinContent (ibin, iamplitude);

      ibin++;
    }

  cout << "Baseline: " << baseline << endl;
  cout << "Time: " << time << " ns" << endl;
  cout << "Amplitude: " << amplitude << " mV" << endl;
  cout << "Charge: " << charge / 1000.0 << " V.ns" << endl;
  cout << "Q/A: " << charge / 1000.0 / amplitude << " us" << endl;

  return hpulse;
}
