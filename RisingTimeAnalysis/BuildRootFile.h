// BuildRootFile.h is derivated from RootEventBuilder.cc @ CC Lyon, written for the reflectometry analysis

struct pulse_t
{
  int slot ;
  int channel ;
  double rising_time ;
  double amplitude_first_peak ;
  double amplitude_second_peak ;
  double charge_first_peak ;
  double charge_second_peak ;
  double baseline ;
  double signal_attenuation ;
  double time_difference_CFD ;
  double time_difference_max ;
  int time_first_pic ;
  int time_second_pic ;

  void clear ()
  {
    slot = 0 ;
    channel = 0 ;
    amplitude_first_peak = 0.0 ;
    amplitude_second_peak = 0.0 ;
    charge_first_peak = 0.0 ;
    charge_second_peak = 0.0 ;
    baseline = 0.0 ;
    signal_attenuation = 0. ;
    time_difference_CFD = 0. ;
    time_difference_max = 0. ;
    time_first_pic=0. ;
    time_second_pic=0. ;

  }
} ;


int bin_max_wall=70 ;// from longest cables of slot 7 run 105 (to be updated with a new run for FR
                    // main wall)

std::vector<pulse_t> ReadFile(string filename_){

  std::vector<pulse_t> pulse_list ;

  ifstream file ;

  TH1F *hpulse = new TH1F ("hpulse", "hpulse", 1024, 0, 400) ;

  file.open (filename_) ;
  if (!file.is_open ()){
    cerr << "ERROR: file isn't open !" << endl ;
  }
  else {
    cerr << "file " << filename_ << " open!!" << endl ;
  }

  bool read_pulses = false ;

  size_t pulse_read = 0 ;


  while (!file.eof ()){
    if (pulse_read % 10000 == 0) cout << "Read pulse: " << pulse_read << endl ;

    //if (pulse_read > 1000) break ;

    string file_line ;
    getline (file, file_line) ;

    ostringstream oss_slot_channel ;
    oss_slot_channel << "Slot " ;

    // read the event:

    pulse_t new_pulse ;
    new_pulse.clear () ;

    size_t event_found = file_line.find (oss_slot_channel.str ()) ;
    if (event_found != std::string::npos){
      string slot_str, channel_str ;
      int slot, channel ;

      istringstream iss_line (file_line) ;
      iss_line >> slot_str >> slot >> channel_str >> channel ;

      string pulse_line ;
      getline (file, pulse_line) ;

      istringstream iss_pulse (pulse_line) ;

      int isample = 0 ;
      double amplitude=0. ;
      double baseline = 0.0 ;
      double adc = 0.0 ;
      double charge_first_peak = 0.0, charge_second_peak = 0.0 ;
      double signal_attenuation = 0. ;
      double amplitude_first_peak = 0., amplitude_second_peak = 0. ;
      double time_difference_CFD = 0., time_difference_max = 0. ;
      int time_first_pic=0, time_second_pic=0 ;

      double time = 0.0 ;

      vector<double> vec_adc ;

      while (iss_pulse >> adc >> ws){
        vec_adc.push_back(adc) ;
        if (isample < 150) baseline += adc ;
        isample++ ;
      }

      baseline /= 150.0 ;

      size_t ibin = 1 ;

      for (vector<double>::iterator it = vec_adc.begin() ; it != vec_adc.end() ; ++it){
        double iamplitude = (*it - baseline)  * adc_to_mv ;


        if (iamplitude > amplitude){
          amplitude = iamplitude ;
          time = ibin * tdc_to_ns ;
        }

        hpulse->SetBinContent (ibin, iamplitude) ;

        ibin++ ;
      }

      int bin_max1=0 ; int bin_max2=0 ;
      amplitude_first_peak=get_pic(1,hpulse,bin_max1) ;
      amplitude_second_peak = get_pic(bin_max1+trigger_second_peak,hpulse,bin_max2) ;
      signal_attenuation=get_attenuation(hpulse) ;
      time_difference_CFD=get_time_difference(hpulse,"CFD",0.25) ;
      time_difference_max=get_time_difference(hpulse,"maximum") ;
      charge_first_peak=get_charge(1,hpulse) ;
      charge_second_peak=get_charge(2,hpulse,bin_max_wall) ;
      time_first_pic=bin_max1 ;
      time_second_pic=bin_max2 ;

      double bas = compute_baseline(hpulse,150) ;
      double rising_time =  get_rising_time(hpulse) ;

      if (amplitude_second_peak>bas+50.) {

        new_pulse.slot = slot ;
        new_pulse.channel = channel ;
        new_pulse.rising_time = rising_time ;
        new_pulse.amplitude_first_peak = amplitude_first_peak ;
        new_pulse.amplitude_second_peak = amplitude_second_peak ;
        new_pulse.charge_first_peak = charge_first_peak / 1000.0 ;
        new_pulse.charge_second_peak = charge_second_peak / 1000.0 ;
        new_pulse.baseline = baseline ;
        new_pulse.signal_attenuation = signal_attenuation ;
        new_pulse.time_difference_CFD = time_difference_CFD ;
        new_pulse.time_difference_max = time_difference_max ;
        new_pulse.time_first_pic = time_first_pic ;
        new_pulse.time_second_pic = time_second_pic ;

        pulse_list.push_back (new_pulse) ;

      }




    }//end if

    pulse_read++ ;
  }//end while on opened file
  file.close() ;


  return pulse_list ;

}



void BuildRootFile(std::vector<pulse_t> pulse_list){

  // Output root file:

  string root_filename = "OutFile.root" ;
  const char* c_root_filename = root_filename.c_str() ;
  TFile *outfile = new TFile (c_root_filename, "recreate") ;
  TTree *pulses = new TTree ("pulses", "tree generated from ascii file") ;

  // Declaration of leaf types
  Int_t pulse_slot ;
  Int_t pulse_channel ;
  Double_t pulse_rising_time ;
  Double_t pulse_amplitude_first_peak ;
  Double_t pulse_amplitude_second_peak ;
  Double_t pulse_charge_first_peak ;
  Double_t pulse_charge_second_peak ;
  Double_t pulse_baseline ;
  Double_t pulse_attenuation ;
  Double_t pulse_attenuation_amp ;
  Double_t pulse_time_difference_CFD ;
  Double_t pulse_time_difference_max ;
  Int_t pulse_time_first_pic ;
  Int_t pulse_time_second_pic ;

  // Branches:
  pulses->Branch ("pulse_slot", &pulse_slot, "pulse_slot/I") ;
  pulses->Branch ("pulse_channel", &pulse_channel, "pulse_channel/I") ;
  pulses->Branch ("pulse_rising_time", &pulse_rising_time, "pulse_rising_pulse/D") ;
  pulses->Branch ("pulse_amplitude_first_peak", &pulse_amplitude_first_peak, "pulse_amplitude_first_peak/D") ;
  pulses->Branch ("pulse_amplitude_second_peak", &pulse_amplitude_second_peak, "pulse_amplitude_second_peak/D") ;
  pulses->Branch ("pulse_charge_first_peak", &pulse_charge_first_peak, "pulse_charge_first_peak/D") ;
  pulses->Branch ("pulse_charge_second_peak", &pulse_charge_second_peak, "pulse_charge_second_peak/D") ;
  pulses->Branch ("pulse_baseline", &pulse_baseline, "pulse_baseline/D") ;
  pulses->Branch ("pulse_attenuation", &pulse_attenuation, "pulse_attenuation/D") ;
  pulses->Branch ("pulse_time_difference_CFD", &pulse_time_difference_CFD, "pulse_time_difference_CFD/D") ;
  pulses->Branch ("pulse_time_difference_max", &pulse_time_difference_max, "pulse_time_difference_max/D") ;
  pulses->Branch ("pulse_time_first_pic", &pulse_time_first_pic, "pulse_time_first_pic/I") ;
  pulses->Branch ("pulse_time_second_pic", &pulse_time_second_pic, "pulse_time_second_pic/I") ;

  for (vector<pulse_t>::iterator it = pulse_list.begin () ; it != pulse_list.end () ; ++it){
    pulse_t a_pulse = *it ;



    pulse_slot = a_pulse.slot ;
    pulse_channel = a_pulse.channel ;
    pulse_rising_time = a_pulse.rising_time ;
    pulse_amplitude_first_peak = a_pulse.amplitude_first_peak ;
    pulse_amplitude_second_peak = a_pulse.amplitude_second_peak ;
    pulse_charge_first_peak = a_pulse.charge_first_peak ;
    pulse_charge_second_peak = a_pulse.charge_second_peak ;
    pulse_baseline = a_pulse.baseline ;
    pulse_attenuation = a_pulse.signal_attenuation ;
    pulse_time_difference_CFD = a_pulse.time_difference_CFD ;
    pulse_time_difference_max = a_pulse.time_difference_max ;
    pulse_time_first_pic = a_pulse.time_first_pic ;
    pulse_time_second_pic = a_pulse.time_second_pic ;

    pulses->Fill () ;

  }

  pulses->Write () ;
  outfile->Close () ;

}
