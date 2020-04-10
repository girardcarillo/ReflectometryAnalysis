// RTD2RootPulsesFunctions.h
// Author: Mathieu Bongrand <bongrand@lal.in2p3.fr>
// Copyright: 2019 (C) SuperNEMO - LAL (IN2P3/CNRS)

#include "CnRoot.h"

using namespace std;

// Constants:

static const uint16_t MAX_CALO_HITS    = 800;
static const uint16_t MAX_WAVEFORM_SAMPLES = 1024;
static const double   TDC_TO_NS = 0.390625;
static const double   ADC_TO_MV = 0.610351563;

// Functions:

double compute_baseline (int16_t (&waveform)[MAX_WAVEFORM_SAMPLES], size_t limit)
{
  double baseline = 0.0;

  for (size_t i = 0; i < limit; i++)
    {
      baseline += waveform [i];
    }
  baseline /= (double) limit;

  return baseline; // adc
}


double compute_charge (int16_t (&waveform)[MAX_WAVEFORM_SAMPLES], size_t start, size_t stop, double baseline)
{
  double charge = 0.0;

  for (size_t i = start; i < stop; i++)
    {
      charge += (waveform [i] - baseline);
    }
  charge *= ADC_TO_MV * TDC_TO_NS;

  return charge; // mV.ns
}


double compute_amplitude (int16_t (&waveform)[MAX_WAVEFORM_SAMPLES], size_t start, size_t stop, double polarity, double baseline)
{
  double amplitude = 0.0;

  for (size_t i = start; i < stop; i++)
    {
      if (amplitude < (polarity * (waveform [i] - baseline)))
	amplitude = polarity * (waveform [i] - baseline);
    }
  amplitude *= ADC_TO_MV;

  return amplitude; // mV
}


size_t compute_max_position (int16_t (&waveform)[MAX_WAVEFORM_SAMPLES], size_t start, size_t stop, double polarity)
{
  double amplitude = 0.0;
  size_t max_position = 0;

  for (size_t i = start; i < stop; i++)
    {
      if (amplitude > (polarity * waveform [i]))
	{
	  amplitude = polarity * waveform [i];
	  max_position = i;
	}
    }
  amplitude *= ADC_TO_MV;

  return max_position; // bin
}


size_t compute_first_max_position (int16_t (&waveform)[MAX_WAVEFORM_SAMPLES], size_t start_position, double threshold, double baseline)
{
  size_t first_max_position = 0;
  double first_max_amplitude = 0.0;

  for (size_t i = start_position; i < MAX_WAVEFORM_SAMPLES; i++)
    {
      double iampl = (waveform [i] - baseline) * ADC_TO_MV;
      double iampl_next = (waveform [i+1] - baseline) * ADC_TO_MV;

      if (iampl < threshold) continue;

      if (iampl > first_max_amplitude && iampl >= iampl_next)
	{
	  first_max_position = i;
	  first_max_amplitude = iampl;
	  break;
	}
    }

  return first_max_position; // bin
}


double compute_cfd_time (int16_t (&waveform)[MAX_WAVEFORM_SAMPLES], size_t start, size_t stop, double polarity, double baseline, double amplitude)
{
  double t_cfd = 0.0;

  for (size_t i = start; i < stop; i++)
    {
      double ampl = polarity * (waveform [i] - baseline) * ADC_TO_MV;


      if (ampl > (polarity * 0.25 * amplitude))
	{
	  t_cfd = i * TDC_TO_NS;
	  break;
	}
    }

  return t_cfd; // ns
}


double compute_cfd_time (int16_t (&waveform)[MAX_WAVEFORM_SAMPLES], size_t reference_position, double baseline)
{
  double t_cfd = 0.0;

  double reference_amplitude = (waveform [reference_position] - baseline) * ADC_TO_MV;

  for (size_t i = reference_position - 100; i < reference_position; i++)
    {
      double iampl = (waveform [i] - baseline) * ADC_TO_MV;

      if (iampl < 0.25 * reference_amplitude) t_cfd = i * TDC_TO_NS;
    }

  return t_cfd; // ns
}
