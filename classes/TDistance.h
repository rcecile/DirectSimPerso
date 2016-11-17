

#ifndef TDISTANCE_H_
#define TDISTANCE_H_
//_________________________________________
/*
 * TDistance.h
 *
 *  Created on: Feb 14, 2011
 *      Author: gorecki
 *
 *
 *
 *      Compute cosmological distances from cosmological parameters
 *      cf Hogg paper  1999
 */
#include <algorithm>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <vector>
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
using namespace std;

class TDistance {
public:
  //constructor
  TDistance();
  // destructor
  ~TDistance();
  // compute angular diameter distance for a given redshift
  double AngularDiameterDistance(double z);
  // compute the comoving distance
  double ComovingDistance(double z);
  // compute the comoving volume
  double ComovingVolume(double zmin, double zmax, double solidangle);
  // compute the comoving volume element
  double ComovingVolumeElement(double z, double solidangle, double dz);
  // convert m to Parsec
  double ConversionToParsec();
  // return sqrt(fOmega_m0 * pow(1 + z, 3) + fOmega_k0 * pow(1 + z, 2)
  // + fOmega_lambda0 * pow(1 + z, 3 * (fwDE + 1)))
  double FunctionE(double z);
  // return Hubble constante value in unit of km/s/Mpc
  double GetHubbleConstante();
  // comput the Hubble distance
  double HubbleDistance();
  // return H0*FunctionE(z)
  double Hubble_vs_z(double z);
  // return luminosity distance
  double LuminosityDistance(double z);
  // Hubble en km/s/Mpc = 1.e3/(3.086.e16) m/s/m
  double ToHubble();
  // set cosmological parameter
  void SetCosmologicalParameter(double Omega_m0 = 0.3, double Omega_lambda0 =
				0.7, double Omega_k0 = 0, double H0 = 70, double wD = -1);
  // set an interpolation function of the integral of  distane FunctionE(z)
  void SetInterpolationDC();
  double conversionHubbleUnit;
  double fOmega_m0;
  double fOmega_lambda0;
  double fOmega_k0;
  double fwDE;
  double fH0;
 private:
  
  ROOT::Math::Interpolator *interpolator_dc;
  
  const double celerite;
  const double conversionParsectom;
  const double zstep;
  
};

#endif /* TDISTANCE_H_ */
