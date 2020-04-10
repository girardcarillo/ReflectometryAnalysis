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

const double amp_bin = 0.6 ; //mV

bool is_labelled_view(string wall,string view) ;

// Defining constants
double c=2.99e1 ;// cm/ns
double theoretical_celerity = 0.69*c ;
double measured_celerity = 0.7*c ;
int channel_tot_number = 13 ;
int slot_tot_number = 20 ;
int PM_tot_number = channel_tot_number*slot_tot_number ;

int transform_index(int x, int y) ;

void RootTimingAnalysis (string run, string wall, string view, bool enable_drawing = 0) {

  int side=-1 ;
  if (wall=="fr") side=1 ;
  else if (wall=="it") side=0 ;
  else cout << "Wrong side value!" << endl ;

  // Root histos:

  TH2D *h2lengthDifference = new TH2D ("length_difference","Length_Difference ; Y ; Z", 22, -1, 21, 15, -1, 14) ;
  TH2D *h2timediff = new TH2D ("time_difference","time_Difference ; Y ; Z", 22, -1, 21, 15, -1, 14) ;

  TH1F *hlengthDifference = new TH1F("hlengthdiff","", 100, 0, 0) ;

  TProfile *plength = new TProfile ("plength","", 100, 1000, 1850, -100, 250) ;
  TProfile *pattenuation = new TProfile ("pattenuation","", 50, 10, 19, 1.6, 5.1) ;
  TProfile *pattenuation_ch = new TProfile ("pattenuation_ch","", 50, 10, 19, 1.6, 5.1) ;

  // Merging .root files
  TList *list = new TList ;
  for (int slot_counter = 0 ; slot_counter < 21 ; ++slot_counter) {
    if (slot_counter!=10) {
      stringstream ss_slot ;
      ss_slot << slot_counter ;
      string str_slot = ss_slot.str() ;

      string filename_ = "root_files_run"+run+"/slot"+str_slot+".root" ;
      // string filename_ = "slot1.root" ;
      const char* c_filename = filename_.c_str() ;

      TFile *file_adress = TFile::Open(c_filename) ;
      TTree *tree_adress = (TTree*)file_adress->Get("pulses") ;
      list->Add(tree_adress) ;
    }
  }

  TFile* outputfile = TFile::Open("output.root", "recreate") ;
  TTree *TPulses = TTree::MergeTrees(list) ;

  cout << "The TTree has " << TPulses->GetEntries () << " total entries."<< endl ;

  // Declaration of leaf types
  Int_t pulse_slot ;
  Int_t pulse_channel ;
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
  TPulses->SetBranchAddress ("pulse_amplitude_first_peak", &pulse_amplitude_first_peak) ;
  TPulses->SetBranchAddress ("pulse_amplitude_second_peak", &pulse_amplitude_second_peak) ;
  TPulses->SetBranchAddress ("pulse_charge_first_peak", &pulse_charge_first_peak) ;
  TPulses->SetBranchAddress ("pulse_charge_second_peak", &pulse_charge_second_peak) ;
  TPulses->SetBranchAddress ("pulse_baseline", &pulse_baseline) ;
  TPulses->SetBranchAddress ("pulse_attenuation", &pulse_attenuation) ;
  TPulses->SetBranchAddress ("pulse_time_difference_CFD", &pulse_time_difference_CFD) ;
  TPulses->SetBranchAddress ("pulse_time_difference_max", &pulse_time_difference_max) ;
  TPulses->SetBranchAddress ("pulse_time_second_pic", &pulse_time_second_pic) ;

  vector <vector <double> > count_map(channel_tot_number) ;
  vector <vector <vector <double> > > attenuations(channel_tot_number) ;
  vector <vector <vector <double> > > attenuationsStatErr(channel_tot_number) ;
  vector <vector <vector <double> > > attenuations_ch(channel_tot_number) ;
  vector <vector <vector <double> > > time_differences(channel_tot_number) ;
  vector <vector <double> > mean_attenuationStatErr(channel_tot_number) ;
  vector <vector <double> > mean_time_difference(channel_tot_number) ;
  vector <vector <double> > mean_attenuation(channel_tot_number) ;
  vector <vector <double> > mean_attenuation_ch(channel_tot_number) ;
  vector <vector <double> > expected_length(channel_tot_number) ;//[cm]
  vector <vector <double> > real_length(channel_tot_number) ;//[cm]
  vector <vector <double> > tab_tmp(channel_tot_number) ;//[cm]
  vector <double> tab_first_column(channel_tot_number) ;
  int cable_length = 325+700 ;//[cm]
  Double_t x[260], y[260] ;

  for (unsigned i = 0 ; i < expected_length.size() ; ++i) {

    count_map[i]=vector<double>(slot_tot_number);
    expected_length[i]=vector<double>(slot_tot_number) ;
    real_length[i]=vector<double>(slot_tot_number) ;
    tab_tmp[i]=vector<double>(slot_tot_number) ;
    time_differences[i]=vector<vector <double> >(slot_tot_number) ;
    attenuations[i]=vector<vector <double> >(slot_tot_number) ;
    attenuationsStatErr[i]=vector<vector <double> >(slot_tot_number) ;
    attenuations_ch[i]=vector<vector <double> >(slot_tot_number) ;
    mean_attenuationStatErr[i]=vector<double>(slot_tot_number) ;
    mean_time_difference[i]=vector<double>(slot_tot_number) ;
    mean_attenuation[i]=vector<double>(slot_tot_number) ;
    mean_attenuation_ch[i]=vector<double>(slot_tot_number) ;

    tab_first_column[i]=cable_length ;

    if (i%2!=0) {
      cable_length+=50 ;
    }

  }


  // create tab representing theoretical designed length
  for (unsigned i = 0 ; i < channel_tot_number ; i++){
    for (unsigned j = 0 ; j < slot_tot_number ; j++){
      count_map[i][j]=0 ;
      mean_attenuationStatErr[i][j]=0 ;
      mean_time_difference[i][j]=0 ;
      mean_attenuation[i][j]=0 ;
      mean_attenuation_ch[i][j]=0 ;

      if (j==0) {
        expected_length[i][j]=tab_first_column[i] ;
      }
      else {
        expected_length[i][j]=expected_length[i][j-1]+25 ;
      }
      tab_tmp[i][j]=expected_length[i][j] ;

    }
  }

  bool test_label_view = is_labelled_view(wall,view) ;

  if (!test_label_view) {
    for (unsigned i = 0 ; i<expected_length.size() ; i++){
      for (unsigned j = 0 ; j<expected_length[i].size() ; j++){
        expected_length[i][j]=tab_tmp[i][slot_tot_number-j-1] ;
      }
    }
  }

  double StatErr = 0. ;

  for (unsigned int i_pulse = 0 ; i_pulse < TPulses->GetEntries() ; i_pulse++){
    TPulses->GetEntry (i_pulse) ;

    if (pulse_channel<channel_tot_number) {//15 channels on but 13 physical

      attenuations[pulse_channel][pulse_slot].push_back(pulse_attenuation) ;
      attenuations_ch[pulse_channel][pulse_slot].push_back(pulse_charge_first_peak/pulse_charge_second_peak) ;
      time_differences[pulse_channel][pulse_slot].push_back(pulse_time_difference_CFD) ;

      StatErr = sqrt(pow((10*(amp_bin/pulse_amplitude_first_peak))/log(10.),2)+pow((10*(amp_bin/pulse_amplitude_second_peak))/log(10.),2)) ;
      attenuationsStatErr[pulse_channel][pulse_slot].push_back(StatErr) ;

      // if (pulse_channel == 12 && pulse_slot == 19) {
      //   cout << pulse_amplitude_first_peak << " " << pulse_amplitude_second_peak << " " << sqrt(pow((10*(amp_bin/pulse_amplitude_first_peak))/log(10.),2)+pow((10*(amp_bin/pulse_amplitude_second_peak))/log(10.),2)) << endl ;
      // }

    }
  }


  TPulses->Write() ;
  outputfile->Close() ;
  delete outputfile ;

  int counter_OM=0 ;
  double lengthDiff = 0 ;
  double att[260],attRMS[260],length[260], att_ch[260] , attRMS_ch[260], attStatErr[260];

  for (unsigned i = 0 ; i < channel_tot_number ; i++){
    for (unsigned j = 0 ; j < slot_tot_number ; j++){
      // // bug with last value of attenuation vectors. Probably something comming from RootEventBuilder.C (CC)
      attenuations[i][j].pop_back() ;
      attenuations_ch[i][j].pop_back() ;
      time_differences[i][j].pop_back() ;

      count_map[i][j] = attenuations[i][j].size() ;
      mean_attenuation[i][j] = TMath::Mean(attenuations[i][j].begin(),attenuations[i][j].end()) ;
      mean_attenuation_ch[i][j] = TMath::Mean(attenuations_ch[i][j].begin(),attenuations_ch[i][j].end()) ;
      mean_time_difference[i][j] = TMath::Mean(time_differences[i][j].begin(),time_differences[i][j].end()) ;
      mean_attenuationStatErr[i][j] = TMath::Mean(attenuationsStatErr[i][j].begin(),attenuationsStatErr[i][j].end()) ;

      real_length[i][j] = (measured_celerity*mean_time_difference[i][j])/2. ;
      lengthDiff = real_length[i][j]-expected_length[i][j] ;
      h2lengthDifference->SetBinContent(j+2,i+2,lengthDiff) ;
      h2timediff->SetBinContent(j+2,i+2,mean_time_difference[i][j]) ;
      hlengthDifference->Fill(lengthDiff) ;

      x[counter_OM] = expected_length[i][j] ;
      y[counter_OM] = lengthDiff ;
      counter_OM++ ;

      int index = transform_index(i,j) ;
      attRMS[index] = TMath::RMS(attenuations[i][j].begin(),attenuations[i][j].end()) ;
      attRMS_ch[index] = TMath::RMS(attenuations_ch[i][j].begin(),attenuations_ch[i][j].end()) ;
      att[index] = 10*log10(1./mean_attenuation[i][j]) ;
      att_ch[index] = 10*log10(mean_attenuation_ch[i][j]) ;
      length[index] = real_length[i][j]/100. ; // Convert cm in m
      attStatErr[index] = mean_attenuationStatErr[i][j] ;

    }
  }

  gStyle->SetOptStat(0) ;
  gStyle->SetOptTitle(0) ;
  gStyle->SetPaintTextFormat("1.1f") ;
  h2timediff->Draw("COLZTEXT") ;

  // TGraphErrors *gr = new TGraphErrors(260, length, att, 0, attRMS) ;
  TGraphErrors *gr = new TGraphErrors(260, length, att, 0, attStatErr) ;
  TGraphErrors *gr2 = new TGraphErrors(260, length, att_ch, 0, attRMS_ch) ;

  double mediangr =  TMath::Median(260,y) ;
  double inf = mediangr-3*mediangr ;
  double sup = mediangr+3*mediangr ;


  for (int i = 0 ; i < 260 ; ++i) {

    pattenuation->Fill(length[i],att[i]) ;
    pattenuation_ch->Fill(length[i],att_ch[i]) ;

    if (y[i] < sup) {
      if (y[i] > inf) {

        plength->Fill(x[i],y[i]) ;

      }
    }
  }

  if (enable_drawing) {

    TCanvas *canvas = new TCanvas ("canvas", "canvas",10,10,2000,1000) ;
    gStyle->SetOptStat(0) ;
    gStyle->SetOptFit(0) ;

    config_profile(plength, "", "l ^{d} (cm)", "#Delta L (cm)","",2,kCyan+2) ;
    TF1 *fit = new TF1("linfit", "pol1") ;
    plength->Fit("linfit") ;
    plength->GetFunction("linfit")->SetLineColor(kOrange+7) ;

    TLine *hline = new TLine(1000,0,1850,0) ;
    hline->SetLineStyle(7) ;
    hline->SetLineColor(kGray+3) ;
    hline->Draw("SAME") ;

    TLegend* legend1 = new TLegend(0.1,0.677,0.383,0.997) ;

    float chi2lin = plength->GetFunction("linfit")->GetChisquare() ;
    float ndflin = fit->GetNDF() ;

    legend1->AddEntry(plength,"Measured cable","lep") ;
    legend1->AddEntry(plength->GetFunction("linfit"),"Linear fit","l") ;
    legend1->AddEntry((TObject*)0,Form("#chi^{2}/ndf  %1.1f/%1.f",chi2lin,ndflin),"") ;
    legend1->AddEntry((TObject*)0,Form("#alpha          %1.3f #pm %1.3f cm/cm",fit->GetParameter(1),fit->GetParError(1)),"") ;
    legend1->AddEntry((TObject*)0,Form("#beta         %1.1f #pm %1.1f cm",fit->GetParameter(0),fit->GetParError(0)),"") ;
    legend1->Draw() ;


    canvas->SaveAs("plots/profile_length.eps") ;

    // gStyle->SetOptStat(0) ;
    // gStyle->SetPaintTextFormat("1.1f") ;
    // config_histo2D(h2lengthDifference, "Difference between real and expected lengths", "Column", "Row","colztext") ;canvas->SaveAs("plots/length_differences_maps.pdf") ;


    TCanvas *canvas1 = new TCanvas ("canvas1", "canvas1",10,10,2000,1000) ;
    config_histo1D(hlengthDifference,"","#Delta L (cm)","Counts",3,1,kOrange+7) ;
    TF1 *f1 = new TF1("f1","gaus",hlengthDifference->GetMean()-1.5*hlengthDifference->GetStdDev(),hlengthDifference->GetMean()+1.5*hlengthDifference->GetStdDev()) ;
    hlengthDifference->Fit("f1","0") ;

    TLine *hline1 = new TLine(0,0,0,48) ;
    hline1->SetLineStyle(7) ;
    hline1->SetLineColor(kGray+3) ;
    hline1->Draw("SAME") ;

    f1->SetLineColor(kGreen+3) ;
    f1->SetLineStyle(9) ;
    f1->Draw("SAME") ;


    TLegend* legend = new TLegend(0.1,0.698,0.486,0.992) ;

    float chi2gaus = hlengthDifference->GetFunction("f1")->GetChisquare() ;
    float ndfgaus = f1->GetNDF() ;

    legend->AddEntry(hlengthDifference,"#Delta L","l") ;
    legend->AddEntry(f1,"Gaussian fit","l") ;
    legend->AddEntry((TObject*)0,Form("#chi^{2}/ndf  %1.1f/%1.f",chi2gaus,ndfgaus),"") ;
    legend->AddEntry((TObject*)0,Form("Mean  %1.1f #pm %1.1f cm",f1->GetParameter(1),f1->GetParError(1)),"") ;
    legend->AddEntry((TObject*)0,Form("#sigma         %1.1f #pm %1.1f cm",f1->GetParameter(2),f1->GetParError(2)),"") ;
    legend->Draw() ;

    canvas1->SaveAs("plots/length_differences_histo.eps") ;

    gStyle->SetOptFit(0) ;
    // gStyle->SetLegendBorderSize(0) ;
    TCanvas *canvas2 = new TCanvas ("canvas2", "canvas2",10,10,2000,1000) ;

    config_profile(pattenuation, "", "l (m)", "A (dB)","",1,kCyan+2) ;
    config_profile(pattenuation_ch, "", "l (m)", "A (dB)","SAME",1,kMagenta+2) ;
    pattenuation->GetYaxis()->SetRangeUser(1.6,5.1) ;

    TF1 *fit1 = new TF1("linfit1", "pol1") ;
    pattenuation->Fit("linfit1") ;
    pattenuation->GetFunction("linfit1")->SetLineColor(kOrange+7) ;
    fit1->SetLineColor(kOrange+7) ;
    fit1->SetLineWidth(1) ;
    pattenuation->GetFunction("linfit1")->SetLineWidth(1) ;
    TF1 *fit2 = new TF1("linfit2", "pol1") ;
    pattenuation_ch->Fit("linfit2") ;
    pattenuation_ch->GetFunction("linfit2")->SetLineColor(kGreen+2) ;
    fit2->SetLineColor(kGreen+2) ;
    fit2->SetLineWidth(1) ;
    pattenuation_ch->GetFunction("linfit2")->SetLineWidth(1) ;

    TLegend* legend2 = new TLegend(0.112,0.59,0.406,0.997) ;

    legend2->AddEntry(pattenuation,"Amplitude attenuation","lep") ;
    legend2->AddEntry(pattenuation_ch,"Charge attenuation","lep") ;
    legend2->AddEntry(fit1,"Linear fit","l") ;
    legend2->AddEntry((TObject*)0,Form("#alpha^{R, amp}_{att}    %1.3f #pm %1.3f dB/m",fit1->GetParameter(1),fit1->GetParError(1)),"") ;
    legend2->AddEntry((TObject*)0,Form("f_{r}^{amp}         %1.3f #pm %1.3f dB",fit1->GetParameter(0),fit1->GetParError(0)),"") ;
    legend2->AddEntry(fit2,"Linear fit","l") ;
    legend2->AddEntry((TObject*)0,Form("#alpha^{R, ch}_{att}      %1.3f #pm %1.3f dB/m",fit2->GetParameter(1),fit2->GetParError(1)),"") ;
    legend2->AddEntry((TObject*)0,Form("f_{r}^{ch}         %1.3f #pm %1.3f dB",fit2->GetParameter(0),fit2->GetParError(0)),"") ;
    legend2->Draw() ;

    // // // // // multigraph :
    // TMultiGraph *mg = new TMultiGraph();

    // gr->SetMarkerColor(1);
    // gr->SetMarkerStyle(20);
    // gr->SetMarkerSize(0.2);
    // gr->SetLineColor(kCyan+2);
    // gr->SetLineWidth(1);
    // gr->SetLineStyle(1);
    // mg->Add(gr,"AP");


    // gr2->SetMarkerColor(1);
    // gr2->SetMarkerStyle(20);
    // gr2->SetMarkerSize(0.2);
    // gr2->SetLineColor(kMagenta+2);
    // gr2->SetLineWidth(1);
    // gr2->SetLineStyle(1);
    // mg->Add(gr2,"AP");

    // TF1 *fit1 = new TF1("linfit1", "pol1") ;
    // gr->Fit("linfit1") ;
    // gr->GetFunction("linfit1")->SetLineColor(kOrange+7) ;
    // fit1->SetLineColor(kOrange+7) ;
    // fit1->SetLineWidth(1) ;
    // gr->GetFunction("linfit1")->SetLineWidth(1) ;
    // TF1 *fit2 = new TF1("linfit2", "pol1") ;
    // gr2->Fit("linfit2") ;
    // gr2->GetFunction("linfit2")->SetLineColor(kGreen+2) ;
    // fit2->SetLineColor(kGreen+2) ;
    // fit2->SetLineWidth(1) ;
    // gr2->GetFunction("linfit2")->SetLineWidth(1) ;

    // mg->Draw("a");
    // mg->GetYaxis()->SetTitleOffset(0.5) ;
    // mg->GetYaxis()->SetTitleSize(0.048);
    // mg->GetXaxis()->SetTitleSize(0.048);
    // mg->GetXaxis()->SetTitle("l (m)");
    // mg->GetYaxis()->SetTitle("A_{R} (dB)");


    // TLegend* legend2 = new TLegend(0.144,0.631,0.390,0.879) ;

    // legend2->AddEntry(gr,"Amplitude attenuation","lep") ;
    // legend2->AddEntry(gr2,"Charge attenuation","lep") ;
    // legend2->AddEntry(fit1,"Linear fit","l") ;
    // legend2->AddEntry((TObject*)0,Form("#alpha^{R, amp}_{att}    %1.3f #pm %1.3f dB/m",fit1->GetParameter(1),fit1->GetParError(1)),"") ;
    // legend2->AddEntry(fit2,"Linear fit","l") ;
    // legend2->AddEntry((TObject*)0,Form("#alpha^{R, ch}_{att}      %1.3f #pm %1.3f dB/m",fit2->GetParameter(1),fit2->GetParError(1)),"") ;
    // legend2->Draw() ;

    // // // //



    canvas2->SaveAs("plots/attenuation.eps") ;

  }


}//end macro




bool is_labelled_view(string wall,string view){
  bool test ;
  if (view=="back") {
    if (wall=="it") {
      test=0 ;
    }
    else if (wall=="fr") {
      test=1 ;
    }
    else {
      test=0 ;
      cout << "ERROR: wrong value for wall!" << endl ;
      exit(1) ;
    }
  }
  else if (view=="front") {
    if (wall=="it") {
      test=1 ;
    }
    else if (wall=="fr") {
      test=0 ;
    }
    else {
      test=0 ;
      cout << "ERROR: wrong value for wall!" << endl ;
      exit(1) ;
    }
  }
  else {
    test=0 ;
    cout << "ERROR: wrong value for view!" << endl ;
    exit(1) ;
  }

  if (!test) {
    cout << "The wall isn't viewed from electronics" << endl ;
  }
  else if (test) {
    cout << "The wall is correctly viewed from electronics" << endl ;
  }

  return test ;
}

int transform_index(int x, int y){

  int index = -1 ;

  int matrice [channel_tot_number][slot_tot_number] ;
  for (int i = 0; i < channel_tot_number; ++i) {
    for (int j = 0; j < slot_tot_number; ++j) {
      matrice[i][j] = i*slot_tot_number+j ;
    }
  }

  if (x<channel_tot_number && y<slot_tot_number) {
    index = matrice[x][y] ;
  }
  else {
    cout << "Function transform_index(int channel, int slot): Bad channel or slot number!!" << endl ;
  }

  return index ;
}


// double mean(const vector<double> v) {

//   double result ;
//   double sumVal = 0.0 ;            //Calculating sum of all values

//   for (int i = 0 ; i < v.size() ; ++i) {
//     sumVal = sumVal + v.at(i) ;
//   }

//   result = sumVal / v.size() ; //Calculating mean

//   return result ;
// }


// double stddev(const vector<double> v) {
//   double total = 0. ;
//   int length = v.size() ;
//   double mean_value = mean(v) ;
//   for (int i = 0 ; i < v.size() ; ++i) {      //Calcuating standard deviation
//     total = total + (v.at(i) - mean_value)*(v.at(i) - mean_value) ;
//   }

//   return sqrt(total / length) ;
// }
