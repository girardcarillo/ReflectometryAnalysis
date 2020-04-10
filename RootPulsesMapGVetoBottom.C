// RootPulsesMap.C
// Author: Cloé Girard-Carillo <girardcarillo@lal.in2p3.fr>

// .x RootPulsesMapGVeto.C+("118","france","tunnel")

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

using namespace std;

//----------------------------//


void config_histo1D(TH1F *histo,const char * draw_type,const char * XaxisTitle=0, const char *YaxisTitle=0,int linewidth=0,int linestyle=0,Color_t linecolor=0);
void config_histo2D(TH2D *histo, const char *histoTitle, const char * XaxisTitle, const char *YaxisTitle,const char * draw_type);
bool is_labelled_view(string wall,string view);
TH2D* from_mean(vector <vector <double> > value,vector <vector <double> > lengths,vector <vector <double> > counts, const char* title, const char* axis);

double c=2.99e1;// cm/ps
double theoretical_celerity = 0.69*c;
int channel_tot_number=16;
int slot_tot_number=21;

void RootPulsesMapGVetoBottom (string run, string side)
{
  gROOT->Reset ();
  gStyle->SetOptStat (0);
  gStyle->SetPaintTextFormat("1.2f");
  gStyle->SetLineStyleString(11,"100 20");

  ofstream corrected_times;
  corrected_times.open("corrected_times.txt");
  if (!corrected_times.is_open()) cout << "File corrected_times not open!" << endl;

  int w=0;int s;
  if (side=="france") s=1;
  else if (side=="italy") s=0;
  else cout << "Wrong side value!" << endl;

  TCanvas *canvas = new TCanvas ("canvas", "canvas", 200, 200, 1100, 800);

  // Root histos:
  TH2D *h2_counts = new TH2D ("counts","counts; Y; Z",  slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_amplitude = new TH2D ("h2_amplitude_first_peak","Amplitude_First_Peak; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_charge_first_peak = new TH2D ("h2_charge_first_peak","Charge_First_Peak; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_charge_second_peak = new TH2D ("h2_charge_second_peak","Charge_Second_Peak; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_amplitude_second_peak = new TH2D ("h2_amplitude_second_peak","Amplitude_Second_Peak; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_attenuation_amplitude = new TH2D ("attenuation_amplitude","Attenuation_Amplitude; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_attenuation_charge = new TH2D ("attenuation_amplitude","Attenuation_Amplitude; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_time_difference_CFD = new TH2D ("time_diff_CFD","time_diff_CFD; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_time_difference_max = new TH2D ("time_diff_max","time_diff_max; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_rms_attenuation = new TH2D ("rms_attenuation","Rms_Attenuation; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_rms_time_difference_CFD = new TH2D ("rms_time_difference_CFD","Rms_Time_Difference_CFD; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_rms_time_difference_max = new TH2D ("rms_time_difference_max","Rms_Time_Difference_max; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_expected_length = new TH2D ("expected_length","expected_length; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_time_difference = new TH2D ("expected_timing","Expected_Timing; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_length_difference = new TH2D ("length_difference","Length_Difference; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);
  TH2D *h2_real_length = new TH2D ("real_length","Real_Length; Y; Z", slot_tot_number+1, -1, slot_tot_number, channel_tot_number+2, -1, channel_tot_number+1);

  TH1F *h_time_test = new TH1F ("test","channel 0 slot 0",100,98,99);

  TH1F *h_time_max_slot = new TH1F ("time_CFD_slot","Time difference (max)",22,-1,21);
  TH1F *h_time_CFD_slot = new TH1F ("time_max_slot","Time difference (CFD)",22,-1,21);
  TH1F *h_time_max_ch = new TH1F ("time_CFD_ch","Time difference (max)",17,-1,16);
  TH1F *h_time_CFD_ch = new TH1F("time_max_ch","Time difference (CFD)",17,-1,16);

  TH1F *h_amplitude_second_peak = new TH1F("amplitude_first_peak","Amplitude",100,200,450);
  TH1F *h_amplitude_first_peak = new TH1F("amplitude_second_peak","Amplitude",100,200,450);
  TH1F *h_charge_first_peak_slot = new TH1F("charge_first_peak","Charge first peak",100,9,12);
  TH1F *h_charge_second_peak_slot = new TH1F("charge_second_peak","Charge second peak",100,9,12);


  TH1F** h_corrected_times = new TH1F*[channel_tot_number]; for (int i=0;i<channel_tot_number;i++) { h_corrected_times[i] = new TH1F(Form("h_corrected_times_%d",i),"Corrected times",22,-1,21); }
  TH1F** h_attenuation_amplitude_slot = new TH1F*[channel_tot_number]; for (int i=0;i<channel_tot_number;i++) { h_attenuation_amplitude_slot[i] = new TH1F(Form("h_attenuation_amplitude_slot_%d",i),"Amplitude attenuation",22,-1,21); }
  TH1F** h_attenuation_charge_slot = new TH1F*[channel_tot_number]; for (int i=0;i<channel_tot_number;i++) { h_attenuation_charge_slot[i] = new TH1F(Form("h_attenuation_charge_slot_%d",i),"Charge attenuation",22,-1,21); }
  TH1F** h_amplitude_slot = new TH1F*[channel_tot_number]; for (int i=0;i<channel_tot_number;i++) { h_amplitude_slot[i] = new TH1F(Form("h_amplitude_slot_%d",i),"Amplitude attenuation",22,-1,21); }
  TH1F** h_charge_slot = new TH1F*[channel_tot_number]; for (int i=0;i<channel_tot_number;i++) { h_charge_slot[i] = new TH1F(Form("h_charge_slot_%d",i),"Charge second peak",22,-1,21); }
  TH1F** h_time_CFD = new TH1F*[channel_tot_number]; for (int i=0;i<channel_tot_number;i++) { h_time_CFD[i] = new TH1F(Form("h_time_CFD_%d",i),"Time difference between 2 peaks",22,-1,21); }

  TList *list = new TList;

  stringstream ss_side,ss_wall;
  ss_side << s;ss_wall << w;

  string filename_ = string("root_files_run")+run+string("/Run")+run+string("_s")+ss_side.str()+string("_w")+ss_wall.str()+string(".root");

  cout << filename_ << endl;
  const char* c_filename = filename_.c_str();

  TFile *file_adress = TFile::Open(c_filename);
  TTree *tree_adress = (TTree*)file_adress->Get("pulses");
  list->Add(tree_adress);


  TFile* outputfile = TFile::Open("output.root", "recreate");
  TTree *TPulses = TTree::MergeTrees(list);

  cout << "The TTree has " << TPulses->GetEntries () << " total entries."<< endl;

  // Declaration of leaf types
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

  vector <vector <double> > count_map(channel_tot_number);
  vector <vector <double> > amplitude_map(channel_tot_number);
  vector <vector <double> > charge_first_peak_map(channel_tot_number);
  vector <vector <double> > charge_second_peak_map(channel_tot_number);
  vector <vector <double> > amplitude_first_peak_map(channel_tot_number);
  vector <vector <double> > amplitude_second_peak_map(channel_tot_number);
  vector <vector <double> > attenuation_amplitude_map(channel_tot_number);
  vector <vector <double> > attenuation_charge_map(channel_tot_number);
  vector <vector <double> > attenuation_squared_map(channel_tot_number);
  vector <vector <double> > time_difference_CFD_squared_map(channel_tot_number);
  vector <vector <double> > time_difference_CFD_map_moy(channel_tot_number);//[ns]
  vector <vector <double> > time_difference_CFD_map(channel_tot_number);//[ns]
  vector <vector <double> > time_difference_max_map(channel_tot_number);//[ns]
  vector <vector <double> > tab_length(channel_tot_number);//[cm]
  vector <vector <double> > tab_tmp(channel_tot_number);//[cm]
  vector <double> tab_first_column(channel_tot_number);
  int cable_length = 325+700;//[cm]


  for (unsigned i = 0; i < tab_length.size(); ++i) {
    tab_length[i]=vector<double>(slot_tot_number);
    tab_tmp[i]=vector<double>(slot_tot_number);
    count_map[i]=vector<double>(slot_tot_number);
    amplitude_map[i]=vector<double>(slot_tot_number);
    charge_first_peak_map[i]=vector<double>(slot_tot_number);
    charge_second_peak_map[i]=vector<double>(slot_tot_number);
    amplitude_first_peak_map[i]=vector<double>(slot_tot_number);
    amplitude_second_peak_map[i]=vector<double>(slot_tot_number);
    attenuation_amplitude_map[i]=vector<double>(slot_tot_number);
    attenuation_charge_map[i]=vector<double>(slot_tot_number);
    time_difference_CFD_squared_map[i]=vector<double>(slot_tot_number);
    attenuation_squared_map[i]=vector<double>(slot_tot_number);
    time_difference_CFD_map_moy[i]=vector<double>(slot_tot_number);
    time_difference_CFD_map[i]=vector<double>(slot_tot_number);
    time_difference_max_map[i]=vector<double>(slot_tot_number);

    tab_first_column[i]=cable_length;
    if (i%2!=0) {
      cable_length+=50;
    }
  }


  for (unsigned i = 0; i<tab_length.size(); i++){
    for (unsigned j = 0; j<tab_length[i].size(); j++){
      count_map[i][j]=0;
      amplitude_map[i][j]=0;
      charge_first_peak_map[i][j]=0;
      charge_second_peak_map[i][j]=0;
      amplitude_first_peak_map[i][j]=0;
      amplitude_second_peak_map[i][j]=0;
      attenuation_amplitude_map[i][j]=0;
      attenuation_charge_map[i][j]=0;
      time_difference_CFD_squared_map[i][j]=0;
      attenuation_squared_map[i][j]=0;
      time_difference_CFD_map_moy[i][j]=0;
      time_difference_CFD_map[i][j]=0;
      time_difference_max_map[i][j]=0;
      if (j==0) {
        tab_length[i][j]=tab_first_column[i];
      }
      else {
        tab_length[i][j]=tab_length[i][j-1]+25;
      }
      tab_tmp[i][j]=tab_length[i][j];

    }
  }

  // bool test_label_view = is_labelled_view(wall,view);


  // if (!test_label_view) {
  //   for (unsigned i = 0; i<tab_length.size(); i++){
  //     for (unsigned j = 0; j<tab_length[i].size(); j++){
  //       tab_length[i][j]=tab_tmp[i][slot_tot_number-j-1];
  //     }
  //   }
  // }

  int time_second_pic=0;

  for (unsigned int i_pulse = 0; i_pulse < TPulses->GetEntries(); i_pulse++){
    TPulses->GetEntry (i_pulse);

    // cout << i_pulse << " " << pulse_channel << " " << pulse_slot << endl;

    if (pulse_channel<channel_tot_number) {
      if (pulse_slot==9) {

        count_map[pulse_channel][pulse_slot]++;
        amplitude_map [pulse_channel][pulse_slot] += pulse_amplitude_first_peak;
        charge_first_peak_map [pulse_channel][pulse_slot] += pulse_charge_first_peak;
        charge_second_peak_map [pulse_channel][pulse_slot] += pulse_charge_second_peak;
        amplitude_first_peak_map [pulse_channel][pulse_slot] += pulse_amplitude_first_peak;
        amplitude_second_peak_map [pulse_channel][pulse_slot] += pulse_amplitude_second_peak;
        // attenuation_amplitude_map [pulse_channel][pulse_slot] += pulse_attenuation;
        attenuation_squared_map [pulse_channel][pulse_slot] += pow(pulse_attenuation,2);
        time_difference_CFD_squared_map [pulse_channel][pulse_slot] += pow(pulse_time_difference_CFD,2);
        time_difference_CFD_map [pulse_channel][pulse_slot] += pulse_time_difference_CFD;
        time_difference_max_map [pulse_channel][pulse_slot] += pulse_time_difference_max;

        if (pulse_channel==9&&pulse_slot==0) {
          time_second_pic+=pulse_time_second_pic;
        }

        if (pulse_slot==0&&pulse_channel==9) {
          h_amplitude_first_peak->Fill(pulse_amplitude_first_peak);
          h_amplitude_second_peak->Fill(pulse_amplitude_second_peak);
          h_charge_first_peak_slot->Fill(pulse_charge_first_peak);
          h_charge_second_peak_slot->Fill(pulse_charge_second_peak);
          h_time_test->Fill(pulse_time_difference_CFD);
        }

        //h2_charge_first_peak->Fill (-1.0 * pulse_charge_first_peak);
        if (pulse_amplitude_first_peak > -100.0) continue;

      }
    }
  }
  TPulses->Write();
  outputfile->Close();
  delete outputfile;


  double v_amplitude[20],v_charge[20];

  h2_rms_time_difference_CFD->SetEntries(1);
  for (unsigned i = 0; i<tab_length.size(); i++){//16 channel
    for (unsigned j = 0; j<tab_length[i].size(); j++){//21 slot
      if (j==9) {
        h2_counts->SetBinContent (j+2, i+2,count_map [i][j]);
        h2_amplitude->SetBinContent(j+2,i+2,amplitude_map [i][j]/count_map [i][j]);
        h2_charge_first_peak->SetBinContent(j+2,i+2,charge_first_peak_map [i][j]/count_map [i][j]);
        h2_charge_second_peak->SetBinContent(j+2,i+2,charge_second_peak_map [i][j]/count_map [i][j]);
        h2_amplitude_second_peak->SetBinContent(j+2,i+2,amplitude_second_peak_map [i][j]/count_map [i][j]);
        //h2_attenuation_amplitude->SetBinContent(j+2,i+2,attenuation_amplitude_map [i][j]/count_map[i][j]);
        h2_attenuation_amplitude->SetBinContent(j+2,i+2,(amplitude_first_peak_map[i][j]/amplitude_second_peak_map[i][j]));
        h2_attenuation_charge->SetBinContent(j+2,i+2,(charge_second_peak_map[i][j]/charge_first_peak_map[i][j])/0.83);
        h2_time_difference_CFD->SetBinContent(j+2,i+2,time_difference_CFD_map [i][j]/count_map [i][j]);
        h2_time_difference_max->SetBinContent(j+2,i+2,time_difference_max_map [i][j]/count_map [i][j]);
        h2_time_difference->SetBinContent(j+2,i+2,(2*tab_length[i][j])/theoretical_celerity-time_difference_CFD_map[i][j]/count_map[i][j]);
        h2_length_difference->SetBinContent(j+2,i+2,((theoretical_celerity*(time_difference_CFD_map[i][j]/count_map[i][j]))/2)-tab_length[i][j]);
        h2_real_length->SetBinContent(j+2,i+2,2*((theoretical_celerity*(time_difference_CFD_map[i][j]/count_map[i][j]))/2));
        h2_expected_length->SetBinContent(j+2,i+2,2*tab_length[i][j]);

        // cout << i << " " << j << " " << count_map[i][j] << endl;

        if (j==9) {
          h_time_CFD_ch->SetBinContent(i+2,time_difference_CFD_map [i][j]/count_map [i][j]);
          h_time_max_ch->SetBinContent(i+2,time_difference_max_map [i][j]/count_map [i][j]);
        }

        h_corrected_times[i]->SetBinContent(j+2,h_time_CFD_slot->GetBinContent(j+2)-h_time_CFD_slot->GetBinContent(j+2)/2);
        h_time_CFD[i]->SetBinContent(j+2,time_difference_CFD_map [i][j]/count_map [i][j]);
        h_amplitude_slot[i]->SetBinContent(j+2,amplitude_second_peak_map [i][j]/amplitude_first_peak_map[i][j]);
        h_charge_slot[i]->SetBinContent(j+2,charge_second_peak_map [i][j]/count_map [i][j]);
        h_attenuation_amplitude_slot[i]->SetBinContent(j+2,amplitude_first_peak_map[i][j]/amplitude_second_peak_map[i][j]);
        // corrected_times << side << " " << j << " " << i << " " << h_time_CFD_slot->GetBinContent(j+2)/2 << endl;
      }
    }
  }

  corrected_times.close();

  TH2D* h2_time_difference_CFD_moy = from_mean(time_difference_CFD_map,tab_length,count_map,"mean_time_difference_CFD","mean_time_difference_CFD;X;Z");
  TH2D *h2_attenuation_mean = from_mean(attenuation_amplitude_map,tab_length,count_map,"mean_attenuation_amplitude","mean_attenuation_amplitude;X;Z");

  // canvas->Divide(2,2);
  // canvas->cd ();config_histo2D(h2_counts, "Number of pulses sent to each PMT", "Column", "Row","colztext");canvas->SaveAs("plots/counts_maps.pdf");
  // canvas->cd ();config_histo2D(h2_real_length, "Real lengths of cables (with reflectometry)", "Column", "Row","colztext");canvas->SaveAs("plots/real_lengths_maps.pdf");
  // canvas->cd ();config_histo2D(h2_amplitude, "Mean amplitude", "Column", "Row","colztext");canvas->SaveAs("plots/amplitude_maps.pdf");
  // canvas->cd ();config_histo2D(h2_counts, "Mean charge_first_peak", "Column", "Row","colztext");canvas->SaveAs("plots/charge_first_peak_maps.pdf");
  // canvas->cd ();config_histo2D(h2_counts, "Mean charge_second_peak", "Column", "Row","colztext");canvas->SaveAs("plots/charge_second_peak_maps.pdf");
  // canvas->cd ();config_histo2D(h2_amplitude_second_peak, "Mean amplitude of second peak", "Column", "Row","colztext");canvas->SaveAs("plots/amplitude_second_peak_maps.pdf");
  // canvas->cd ();config_histo2D(h2_attenuation_amplitude, "Mean attenuation in amplitude for each cable", "Column", "Row","colztext");canvas->SaveAs("plots/attenuation_amplitude_maps.pdf");
  canvas->cd ();config_histo2D(h2_attenuation_charge, "Mean attenuation charge of each cable", "Column", "Row","colztext");canvas->SaveAs("plots/attenuation_charge_maps.pdf");
  // canvas->cd (1);config_histo2D(h2_time_difference_CFD, "Mean time difference between 2 peaks (CFD)", "Column", "Row","colztext");canvas->SaveAs("plots/time_difference_CFD_maps.pdf");
  // canvas->cd (2);config_histo2D(h2_time_difference_max, "Mean time difference between 2 peaks (max)", "Column", "Row","colztext");canvas->SaveAs("plots/time_difference_max_maps.pdf");
  // canvas->cd ();config_histo2D(h2_attenuation_mean, "Mean attenuation of each cable with respect to the wall mean", "Column", "Row","colztext");canvas->SaveAs("plots/attenuation_mean_maps.pdf");
  // canvas->cd ();config_histo2D(h2_time_difference_CFD_moy, "Mean time_difference_CFD of each cable with respect to the wall mean", "Column", "Row","colztext");canvas->SaveAs("plots/time_moy_maps.pdf");
  // canvas->cd ();config_histo2D(h2_rms_attenuation, "RMS of mean attenuation of each cable", "Column", "Row","colztext");canvas->SaveAs("plots/rms_attenuation_maps.pdf");
  // canvas->cd ();config_histo2D(h2_rms_time_difference_CFD, "RMS of mean time_difference_CFD of each cable", "Column", "Row","colztext");canvas->SaveAs("plots/rms_time_difference_CFD_maps.pdf");
  // canvas->cd ();config_histo2D(h2_expected_length, "Expected cable lengths", "Column",
  // "Row","colztext");canvas->SaveAs("plots/expected_length_maps.pdf");
  // canvas->cd ();config_histo2D(h2_time_difference, "Difference between real and expected times", "Column", "Row","colztext");canvas->SaveAs("plots/time_differences_maps.pdf");
  // canvas->cd (2);config_histo2D(h2_length_difference, "Difference between real and expected lengths", "Column", "Row","colztext");canvas->SaveAs("plots/length_differences_maps.pdf");

  // canvas->cd();config_histo1D(h_time_max_slot,"","Slot","#Delta t (ns)",2,2,1);h_time_max_slot->GetYaxis()->SetRangeUser(90,170);canvas->SaveAs("plots/time_difference_max_slot.pdf");
  // canvas->cd();config_histo1D(h_time_CFD_slot,"","Slot","#Delta t (ns)",2,1);h_time_CFD_slot->GetYaxis()->SetRangeUser(90,170);canvas->SaveAs("plots/time_difference_CFD_slot.pdf");
  // canvas->cd();config_histo1D(h_time_test,"","#Delta t (ns)","#events",2,1,1);/*h_time_CFD_slot->GetYaxis()->SetRangeUser(90,170)*/;canvas->SaveAs("plots/time_test.pdf");



  // TGraph *gr_amplitude_charge = new TGraph (20, v_amplitude, v_charge);
  // canvas->cd();
  // gr_amplitude_charge->SetTitle("Amplitude vs Charge (second peak)");
  // gr_amplitude_charge->GetXaxis()->SetTitle("Amplitude (mV)");
  // gr_amplitude_charge->GetYaxis()->SetTitle("Charge (mV.ns)");
  // gr_amplitude_charge->SetMarkerStyle(8);
  // gr_amplitude_charge->Draw("AC");
  // canvas->SaveAs("plots/amplitude_charge.pdf");



  // temps corrigés par slot/channel

  // canvas->cd();config_histo1D(h_corrected_times[channel_tot_number-1],"","Slot","#Delta t (ns)",2,1,1);
  // auto legend = new TLegend(0.1,0.8,0.65,0.9);
  // legend->AddEntry(h_corrected_times[channel_tot_number-1],"channel12","l");
  // for (int i = channel_tot_number-2; i>=3; i=i-1){
  //   stringstream ss;
  //   ss << i;
  //   string str = ss.str();
  //   str=string("channel")+str;
  //   const char* c_str = str.c_str();

  //   if (i%2!=0) {
  //     if (channel_tot_number-i!=10) {
  //       config_histo1D(h_corrected_times[i],"SAME","Slot","#Delta t (ns)",2,1,channel_tot_number-i);
  //     }
  //     else {

  //       config_histo1D(h_corrected_times[i],"SAME","Slot","#Delta t (ns)",2,1,46);
  //     }
  //     h_corrected_times[i]->SetLineStyle(7);
  //   }
  //   else {
  //     config_histo1D(h_corrected_times[i],"SAME","Slot","#Delta t (ns)",2,1,h_corrected_times[i+1]->GetLineColor());
  //   }

  //   legend->AddEntry(h_corrected_times[i],c_str,"l");

  // }
  // legend->SetNColumns (5);
  // legend->Draw();
  // canvas->SaveAs("plots/corrected_times.pdf");





  // time difference par slot/channel

  // canvas->cd();config_histo1D(h_time_CFD[channel_tot_number-1],"","Slot","#Delta t (ns)",2,1,1);h_time_CFD[channel_tot_number-1]->GetYaxis()->SetRangeUser(90,180);
  // auto legend = new TLegend(0.1,0.8,0.65,0.9);
  // legend->AddEntry(h_time_CFD[channel_tot_number-1],"channel12","l");
  // for (int i = channel_tot_number-2; i>=0; i=i-1){
  //   stringstream ss;
  //   ss << i;
  //   string str = ss.str();
  //   str=string("channel")+str;
  //   const char* c_str = str.c_str();

  //   if (i%2!=0) {
  //     if (channel_tot_number-i!=10) {
  //       config_histo1D(h_time_CFD[i],"SAME","Slot","#Delta t (ns)",2,1,channel_tot_number-i);
  //     }
  //     else {

  //       config_histo1D(h_time_CFD[i],"SAME","Slot","#Delta t (ns)",2,1,46);
  //     }
  //     h_time_CFD[i]->SetLineStyle(7);
  //   }
  //   else {
  //     config_histo1D(h_time_CFD[i],"SAME","Slot","#Delta t (ns)",2,1,h_time_CFD[i+1]->GetLineColor());
  //   }

  //   legend->AddEntry(h_time_CFD[i],c_str,"l");

  // }
  // legend->SetNColumns (5);
  // legend->Draw();
  // canvas->SaveAs("plots/time_difference_CFD.pdf");



  // // amplitude par slot/channel

  // canvas->cd();canvas->SetLogy();config_histo1D(h_amplitude_slot[channel_tot_number-1],"","Slot","Amplitude (mV)",2,1,1);h_amplitude_slot[channel_tot_number-1]->GetYaxis()->SetRangeUser(.3,.55);
  // auto legend = new TLegend(0.35,0.8,0.9,0.9);
  // legend->AddEntry(h_amplitude_slot[channel_tot_number-1],"channel12","l");
  // for (int i = channel_tot_number-2; i>=0; i=i-1){
  //   stringstream ss;
  //   ss << i;
  //   string str = ss.str();
  //   str=string("channel")+str;
  //   const char* c_str = str.c_str();

  //   if (i%2!=0) {
  //     if (channel_tot_number-i!=10) {
  //       config_histo1D(h_amplitude_slot[i],"SAME","Slot","#Delta t (ns)",2,1,channel_tot_number-i);
  //     }
  //     else {

  //       config_histo1D(h_amplitude_slot[i],"SAME","Slot","#Delta t (ns)",2,1,46);
  //     }
  //     h_amplitude_slot[i]->SetLineStyle(7);
  //   }
  //   else {
  //     config_histo1D(h_amplitude_slot[i],"SAME","Slot","#Delta t (ns)",2,1,h_amplitude_slot[i+1]->GetLineColor());
  //   }

  //   legend->AddEntry(h_amplitude_slot[i],c_str,"l");

  // }
  // legend->SetNColumns (5);
  // legend->Draw();
  // canvas->SaveAs("plots/amplitude_slot.pdf");




  // charge par slot/channel
  // canvas->cd();config_histo1D(h_charge_slot[channel_tot_number-1],"","Slot","Charge (mV.ns)",2,1,1);h_charge_slot[channel_tot_number-1]->GetYaxis()->SetRangeUser(9,11);
  // auto legend = new TLegend(0.35,0.8,0.9,0.9);
  // legend->AddEntry(h_charge_slot[channel_tot_number-1],"channel12","l");
  // for (int i = channel_tot_number-2; i>=0; i=i-1){
  //   stringstream ss;
  //   ss << i;
  //   string str = ss.str();
  //   str=string("channel")+str;
  //   const char* c_str = str.c_str();

  //   if (i%2!=0) {
  //     if (channel_tot_number-i!=10) {
  //       config_histo1D(h_charge_slot[i],"SAME","Slot","#Delta t (ns)",2,1,channel_tot_number-i);
  //     }
  //     else {

  //       config_histo1D(h_charge_slot[i],"SAME","Slot","#Delta t (ns)",2,1,46);
  //     }
  //     h_charge_slot[i]->SetLineStyle(7);
  //   }
  //   else {
  //     config_histo1D(h_charge_slot[i],"SAME","Slot","#Delta t (ns)",2,1,h_charge_slot[i+1]->GetLineColor());
  //   }

  //   legend->AddEntry(h_charge_slot[i],c_str,"l");

  // }
  // legend->SetNColumns (5);
  // legend->Draw();
  // canvas->SaveAs("plots/charge_slot.pdf");






  // attenuation amplitude par slot/channel
  // canvas->cd();config_histo1D(h_attenuation_amplitude_slot[channel_tot_number-1],"","Slot","Amplitude (mV.ns)",2,1,1);//h_attenuation_amplitude_slot[channel_tot_number-1]->GetYaxis()->SetRangeUser(9,11);
  // auto legend = new TLegend(0.35,0.8,0.9,0.9);
  // legend->AddEntry(h_attenuation_amplitude_slot[channel_tot_number-1],"channel12","l");
  // for (int i = channel_tot_number-2; i>=0; i=i-1){
  //   stringstream ss;
  //   ss << i;
  //   string str = ss.str();
  //   str=string("channel")+str;
  //   const char* c_str = str.c_str();

  //   if (i%2!=0) {
  //     if (channel_tot_number-i!=10) {
  //       config_histo1D(h_attenuation_amplitude_slot[i],"SAME","Slot","#Delta t (ns)",2,1,channel_tot_number-i);
  //     }
  //     else {

  //       config_histo1D(h_attenuation_amplitude_slot[i],"SAME","Slot","#Delta t (ns)",2,1,46);
  //     }
  //     h_attenuation_amplitude_slot[i]->SetLineStyle(7);
  //   }
  //   else {
  //     config_histo1D(h_attenuation_amplitude_slot[i],"SAME","Slot","#Delta t (ns)",2,1,h_attenuation_amplitude_slot[i+1]->GetLineColor());
  //   }

  //   legend->AddEntry(h_attenuation_amplitude_slot[i],c_str,"l");

  // }
  // legend->SetNColumns (5);
  // legend->Draw();
  // canvas->SaveAs("plots/amplitude_slot.pdf");





  // canvas->cd();config_histo1D(h_time_max_ch,"","Channel","#Delta t (ns)",2,2,1);canvas->SaveAs("plots/time_difference_max_ch.pdf");
  // canvas->cd();config_histo1D(h_time_CFD_ch,"","Channel","#Delta t (ns)",2,1,1);canvas->SaveAs("plots/time_difference_CFD_ch.pdf");


  // canvas->cd();config_histo1D(h_amplitude_first_peak,"Amplitude (mV)","#events",2,1,"");canvas->SaveAs("plots/amplitude.pdf");

  // canvas->cd();config_histo1D(h_charge_first_peak_slot,"","Charge (mV.ns)","#events",2,1,1);config_histo1D(h_charge_second_peak_slot,"SAME","Charge (mV.ns)","#events",2,2,1);
  // h_charge_first_peak_slot->SetTitle("Channel 0 slot 0");
  // auto legend = new TLegend(0.1,0.8,0.3,0.9);
  // legend->AddEntry(h_charge_first_peak_slot,"charge first peak","l");
  // legend->AddEntry(h_charge_second_peak_slot,"charge second peak","l");
  // legend->Draw();
  // canvas->SaveAs("plots/charge_peak.pdf");


  // canvas->cd();config_histo1D(h_amplitude_second_peak,"","Amplitude (mV)","#events",2,2,1);config_histo1D(h_amplitude_first_peak,"SAME","Amplitude (mV)","#events",2,1,1);
  // h_amplitude_second_peak->SetTitle("Channel 0 slot 0");
  // auto legend = new TLegend(0.7,0.8,0.9,0.9);
  // legend->AddEntry(h_amplitude_first_peak,"amplitude first peak","l");
  // legend->AddEntry(h_amplitude_second_peak,"amplitude second peak","l");
  // legend->Draw();
  // canvas->SaveAs("plots/amplitude_peak.pdf");


  // canvas->Modified();

}//end RootPulsesMap


void config_histo1D(TH1F *histo,const char * draw_type,const char * XaxisTitle=0, const char *YaxisTitle=0,int linewidth=0,int linestyle=0,Color_t linecolor=0)
{
  histo->Draw(draw_type);
  histo->GetYaxis()->SetTitleSize(0.048);
  histo->GetXaxis()->SetTitleSize(0.048);
  histo->GetXaxis()->SetTitle(XaxisTitle);
  histo->GetYaxis()->SetTitle(YaxisTitle);
  histo->SetLineWidth(linewidth);
  histo->SetLineStyle(linestyle);
  histo->SetLineColor(linecolor);
}


void config_histo2D(TH2D *histo, const char *histoTitle, const char * XaxisTitle, const char *YaxisTitle,const char * draw_type)
{
  histo->Draw(draw_type);
  histo->GetYaxis()->SetTitleSize(0.048);
  histo->GetXaxis()->SetTitleSize(0.048);
  histo->SetTitle(histoTitle);
  histo->GetXaxis()->SetTitle(XaxisTitle);
  histo->GetYaxis()->SetTitle(YaxisTitle);
}



bool is_labelled_view(string wall,string view){
  bool test;
  if (view=="back") {
    if (wall=="it") {
      test=0;
    }
    else if (wall=="fr") {
      test=1;
    }
    else {
      test=0;
      cout << "ERROR: wrong value for wall!" << endl;
      exit(1);
    }
  }
  else if (view=="front") {
    if (wall=="it") {
      test=1;
    }
    else if (wall=="fr") {
      test=0;
    }
    else {
      test=0;
      cout << "ERROR: wrong value for wall!" << endl;
      exit(1);
    }
  }
  else {
    test=0;
    cout << "ERROR: wrong value for view!" << endl;
    exit(1);
  }

  if (!test) {
    cout << "The wall isn't viewed from electronics" << endl;
  }
  else if (test) {
    cout << "The wall is correctly viewed from electronics" << endl;
  }

  return test;
}

TH2D* from_mean(vector <vector <double> > value,vector <vector <double> > lengths,vector <vector <double> > counts, const char* title, const char* axis){
  double mean = 0;
  int count=0;
  TH2D *histo2D = new TH2D (title,axis, 22, -1, 21, 15, -1, 14);

  for (unsigned i = 0; i<lengths.size(); i++){
    for (unsigned j = 0; j<lengths[i].size(); j++){
      mean += value[i][j]/counts[i][j];
      count++;
    }
  }

  mean=mean/count;

  for (unsigned i = 0; i<lengths.size(); i++){
    for (unsigned j = 0; j<lengths[i].size(); j++){
      histo2D->SetBinContent(j+2,i+2,value[i][j]/counts[i][j]-mean);
    }
  }

  return histo2D;
}
