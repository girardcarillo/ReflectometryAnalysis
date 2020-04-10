// Author: Cloé Girard-Carillo <girardcarillo@lal.in2p3.fr>
#include <iostream>

#include <TH1F.h>

using namespace std;

int trigger_second_peak=75;


const double tdc_to_ns = 0.390625;
const double adc_to_mv = 0.610351563;

TGraph *interpolation(TH1F *histo,int coef,int start_bin,int end_bin,double &maximum,double &minimum){
  int Nbins=end_bin-start_bin;
  TGraph *ginter=new TGraph(coef*Nbins);

  double x[Nbins];
  double y[Nbins];
  for (int i=start_bin; i < end_bin; ++i) {
    x[i-start_bin]=histo->GetBinCenter(i);
    y[i-start_bin]=histo->GetBinContent(i);
  }

  TGraph *graph = new TGraph (Nbins,x,y);

  // Pulse interpolation:
  // spline or math::interpolator

  double x_min = TMath::MinElement (graph->GetN(),graph->GetX());
  double x_max = TMath::MaxElement (graph->GetN(),graph->GetX());

  size_t iinter = 0;

  for (double x = x_min; x < x_max; x += 1/(double)coef) {
    double x_amplitude = graph->Eval (x,0,"S");

    ginter-> SetPoint (iinter, x, x_amplitude);
    iinter++;
  }

  // TGraph *graph_test = new TGraph (Nbins,ginter->GetY(),ginter->GetX());

  // cout << "max = " << TMath::MaxElement(graph_test->GetN(),graph_test->GetX()) << endl;

  // maximum=graph_test->Eval(graph->GetMaximum(),0,"S");
  // cout << maximum << endl;

  return ginter;
}


// double get_pic(int start_bin,TH1F *histo,int &bin_max){
//   if (start_bin>histo->GetNbinsX()) {
//     cout << "get_pic function from Reflectometryanalysis.h: start_bin is greater than the number of bins in the histogram!!" << endl;
//     return 0;
//   }
//   graph = interpolation(histo,10,0,1023);
//   double maximum = TMath::MaxElement(graph->GetN(),graph->GetY());
//   cout << TMath::MaxElement(graph->GetN(),graph->GetY()) << endl;
//   return maximum;
// }



double get_pic(int start_bin,TH1F *histo,int &bin_max){

  if (start_bin>histo->GetNbinsX()) {
    cout << "get_pic() function from Reflectometryanalysis.h: start_bin is greater than the number of bins in the histogram!!" << endl;
    return 0;
  }

  double maximum = histo->GetBinContent(start_bin);

  for (int nbin = start_bin; nbin <= histo->GetNbinsX(); ++nbin) {

    if (histo->GetBinContent(nbin)>maximum) {
      maximum = histo->GetBinContent(nbin);
      bin_max=nbin;
    }

  }

  return maximum;
}


double get_charge(int peak,TH1F* histo,int bin_max_wall=0){
  int bin_max1=0;int bin_max2=0;
  double maximum1 = get_pic(1,histo,bin_max1);
  double maximum2 = get_pic(bin_max1+trigger_second_peak,histo,bin_max2);
  double charge = 0;

  if(peak!=1&&peak!=2) {
    cout << "Bad peak number!" << endl;
    exit(1);
  }
  if (peak==1) {
    charge=histo->Integral(bin_max1-50,bin_max1+100);
  }
  else if (peak==2) {
    charge=histo->Integral(bin_max2-75,bin_max2+bin_max_wall);
  }
  return charge;
}

// double get_charge(int peak,TH1F* histo,int bin_max_wall=0){
//   int bin_max1=0;int bin_max2=0;
//   double maximum1 = get_pic(1,histo,bin_max1);
//   double maximum2 = get_pic(bin_max1+trigger_second_peak,histo,bin_max2);
//   double charge = 0;

//   if(peak!=1&&peak!=2) {
//     cout << "Bad peak number!" << endl;
//     exit(1);
//   }
//   if (peak==1) {
//     charge=histo->Integral(bin_max1-50,bin_max1+100);
//   }
//   else if (peak==2) {
//     charge=histo->Integral(bin_max2+bin_max_wall,1023);
//   }
//   return charge;
// }

double get_time_difference(TH1F *histo, string methode, double fraction=0){

  double time_difference=0;
  Int_t bin_max1=0;
  Int_t bin_max2=0;
  Int_t XCFD = 0;
  Int_t Xmax = 0;
  double maximum1 = get_pic(1,histo,bin_max1);
  double maximum2 = get_pic(bin_max1+trigger_second_peak,histo,bin_max2);
  Double_t CFD = fraction*maximum2;

  // if (maximum2 > 0.1*maximum1) {

    if (methode=="CFD") {
      histo->GetBinWithContent(CFD,XCFD,bin_max1+trigger_second_peak,bin_max2,1000);
      time_difference = fabs((bin_max1-XCFD)*tdc_to_ns);
    }

    else if (methode=="maximum") {
      histo->GetBinWithContent(maximum2,Xmax,bin_max1+trigger_second_peak,bin_max2,1000);
      time_difference=fabs((bin_max1-Xmax)*tdc_to_ns);
    }

  // }

  // else {
  //   time_difference = -1 ;
  // }


  return time_difference;
}

double get_attenuation(TH1F *histo){
  double attenuation=0;
  int bin_max1=0;
  int bin_max2=0;
  double maximum1 = get_pic(1,histo,bin_max1);
  double maximum2 = get_pic(bin_max1+trigger_second_peak,histo,bin_max2);
  attenuation = maximum2/maximum1;
  return attenuation;
}

// baseline computation on the ''nbins'' first bins
double compute_baseline(TH1F *histo, int nbins){

  double baseline = 0. ;

  for (int i = 0; i < histo->GetNbinsX(); ++i) {
    if (i < nbins) {

      baseline += histo->GetBinContent(i) ;

    }
  }

  baseline /= nbins ;

  return baseline ;
}

double get_rising_time(TH1F *histo){

  double rising_time = -1.;

  int bin_max1 = -1 ;
  int bin_max2 = -1 ;
  int bin_max1_10 = -1 ;
  int bin_max1_90 = -1 ;
  double maximum1 = get_pic(1,histo,bin_max1);
  double maximum2 = get_pic(bin_max1+trigger_second_peak,histo,bin_max2);

  bool found_edge_10 = 0 ;
  bool found_edge_90 = 0 ;

  for (int i = bin_max1+trigger_second_peak; i <= bin_max2; ++i) {
    if (found_edge_10 == 0 && histo->GetBinContent(i) > 0.1*maximum2) {

      found_edge_10 = 1 ;
      bin_max1_10 = i ;

    }
    if (found_edge_90 == 0 && histo->GetBinContent(i) > 0.9*maximum2) {

      found_edge_90 = 1 ;
      bin_max1_90 = i ;

    }
  }

  if (bin_max1_10 != -1 && bin_max1_90 != -1) {
    rising_time = (bin_max1_90-bin_max1_10)*tdc_to_ns ;
  }

  return rising_time ;

}
