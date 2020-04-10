// RootEventBuilder.C
// Author: Cloé Girard-Carillo <girardcarillo@lal.in2p3.fr>

// example :
// root 'BuildRootFile.cc("data/RunData_12042018_15h22min33s_Reflecto_Slot0_Ascii.dat")'

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TROOT.h"
#include <TH2D.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TLegend.h>

#include "/home/girardcarillo/Workdir/Analyses/reflectometry_tests/Reflectometryanalysis.h"
#include "BuildRootFile.h"

using namespace std;


void RootBuilderASCII (string filename_){

  std::vector<pulse_t> pulse_list = ReadFile(filename_) ;
  BuildRootFile(pulse_list) ;

}// end BuildRootFile()
