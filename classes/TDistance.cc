/*
 * TDistance.C
 *
 *  Created on: Feb 14, 2011
 *      Author: gorecki
 */

#include "TDistance.h"
#include "TMath.h"
#include <iostream>

TDistance::TDistance() : /*constructor*/
  celerite(3.e8), conversionParsectom(3.08568025e16), zstep(0.001) {
  //Hubble en km/s/Mpc = 1.e3/(3.086.e16) m/s/m
  conversionHubbleUnit = 1.e3 / (conversionParsectom * 1.e6);
}

TDistance::~TDistance() { //destructor
  delete interpolator_dc;
}

void TDistance::SetCosmologicalParameter(double Omega_m0, double Omega_lambda0,
					 double Omega_k0, double H0, double wDE)
{//set cosmological parameter an fill a	vector D_C
  fOmega_m0 = Omega_m0;
  fOmega_lambda0 = Omega_lambda0;
  fOmega_k0 = Omega_k0;
  fwDE = wDE;
  if (TMath::Abs(fOmega_k0 + fOmega_lambda0 + fOmega_m0 - 1) > 0.01) {
    cout << " Attention, sum of the density Omega is not  = 0 " << endl
	 << " it's equal to" << fOmega_k0 + fOmega_lambda0 + fOmega_m0
	 << endl;
  }
  fH0 = H0 * conversionHubbleUnit; //en m/s/m
  SetInterpolationDC();
}

double TDistance::FunctionE(double z) {
  double E = sqrt(fOmega_m0 * pow(1 + z, 3) + fOmega_k0 * pow(1 + z, 2)
		  + fOmega_lambda0 * pow(1 + z, 3 * (fwDE + 1)));
  //	cout<<" FunctionE = "<<E<<" fwde = "<<fwDE<<endl;
  return E;
}

double TDistance::Hubble_vs_z(double z) {//en m/s/m
  return fH0 * FunctionE(z);
}

double TDistance::ComovingDistance(double z) {//en m/s/m
  double dm = 0;
  double dc = 0;
  double dh = HubbleDistance();
  dc = interpolator_dc -> Eval(z);
  int Case = 0;
  if (fOmega_k0 == 0)
    dm = dc;
  if (fOmega_k0 < 0) {
    dm = dh / sqrt(-fOmega_k0) * TMath::Sin(sqrt(-fOmega_k0) * dc / dh);
    Case = 1;
  }
  if (fOmega_k0 > 0) {
    dm = dh / sqrt(fOmega_k0) * TMath::SinH(sqrt(fOmega_k0) * dc / dh);
    Case = 2;
  }
  return dm;
}

double TDistance::HubbleDistance() {//en m
  return celerite / fH0;
}

double TDistance::LuminosityDistance(double z) {//en m/s/m
  //dl = c/H_0 int_0^z dz/E(z)
  return ComovingDistance(z) * (1 + z);
}

double TDistance::AngularDiameterDistance(double z) {//en m/s/m
  return ComovingDistance(z) / (1 + z);
}

double TDistance::ConversionToParsec() {
  //1 m = 1/conversionParsectom pc
  return 1 / conversionParsectom;
}

double TDistance::ComovingVolumeElement(double z, double solidangle, double dz) {
  double return_value = 0;
  return_value = HubbleDistance() * pow(1 + z, 2) * pow(
							AngularDiameterDistance(z), 2) / FunctionE(z) * solidangle * dz;
  return return_value;
}

double TDistance::ComovingVolume(double zmin, double zmax, double solidangle) {//integrate
  //comoving volume element between zmin and zmax with a step zstep =0.01
  double return_value = 0;
  double z = zmin;
  double deltaz = zstep;
  int nstep = int((zmax - zmin) / deltaz);
  for (int i = 0; i < nstep; i++) {
    z = zmin + i * deltaz;
    return_value += ComovingVolumeElement(z + deltaz / 2, solidangle,
					  deltaz);
    z += deltaz;
  }
  if (z < zmax) {
    
    return_value += ComovingVolumeElement(z + deltaz / 2, solidangle,
					  deltaz);
  }
  return return_value;
}

double TDistance::GetHubbleConstante() { //return fH0
  return fH0;
}

double TDistance::ToHubble() {
  // Before called MeterToHubble()
  //Hubble en km/s/Mpc = 1.e3/(3.086.e16) m/s/m
  return conversionHubbleUnit;
}

void TDistance::SetInterpolationDC() {//fill interpolation vector D_C
  vector<double> vect_Integral_dc;
  vector<double> vect_z;
  double zmin = 0;
  double zmax = 12;
  double deltaz = 0.0001;
  int nz = int((zmax - zmin) / deltaz);
  double dc_int = 0;
  vect_z.push_back(0);
  vect_Integral_dc.push_back(0);
  double z = 0;
  double dh = HubbleDistance();
  for (int i = 0; i < nz; i++) {
    z = zmin + deltaz / 2 + i * deltaz;
    dc_int += 1 / FunctionE(z) * deltaz;
    vect_z.push_back(z);
    vect_Integral_dc.push_back(dc_int * dh);
    
  }
  
  interpolator_dc = new ROOT::Math::Interpolator(vect_z, vect_Integral_dc,
						 ROOT::Math::Interpolation::kLINEAR);
}

