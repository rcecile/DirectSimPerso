#ifndef TAPPARENTMAGNITUDE_H_
#define TAPPARENTMAGNITUDE_H_
//_________________________________
/*
 * TApparentMagnitude.h
 *
 *  Created on: Oct 28, 2010
 *      Author: gorecki
 * Modified on: Sep 30, 2012 by JS Ricol
 *
 *      Compute the apparent magnitude
 */

#include <TROOT.h>
#include <TChain.h>
#include "TFile.h"
#include <TH1F.h>
#include "TF1.h"
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TRandom3.h>
#include "TH2F.h"
#include "TGraph.h"
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
#include "TMath.h"
#include "TKcorrection.h"
#include "TDistance.h"
#include "TFilter.h"
#include "TFlux.h"
#include "TSed.h"
#include "TNtuple.h"
using namespace std;

class TApparentMagnitude {

public:
  //constructor
  TApparentMagnitude(bool p=1);
  //destructor
  ~TApparentMagnitude();
  // return the absolute magnitude for a fiven redshift, a given value of ratio= Sum_eo / Sum_ee,
  // spectral type, E(B-V), and filter
  double AbsoluteMagnitude(double redshift, double ratio, int type,
			   double ebv, int ifilter);
  
  // transform imag into abs mag for P(z) grid purpose
  double ComputeAbsMagFromMag(int ifilter,
			      double mag,
			      int type,
			      double z,
			      double ebv);

  // full mag and error process
  void ComputeApparentMagnitudes(double *mag, double *errmag, double mag_abs,
				 int type_s,double z_s, double ebv_s);
  void EvalApparentMagnitudes(double *mag, double *errmag, double mag_abs,
			      int type_s,double z_s, double ebv_s, int ib=-1);

  // return theoretical apparent magnitude
  void ComputeApparentMagnitudeTheorique(double *mag, double MA, int type,
					 double z, double ebv);
  void EvalApparentMagnitudeTheorique(double *mag, double MA, int type,
				      double z, double ebv, int ib=-1);

  void EvalFluxTheorique(double *flux, double MA, int type,
			 double z, double ebv);
  
  // return the uncertainty on the apparent magnitude given by the LSST science book
  void ComputeApparentMagnitudeErrors(double *mag, double *error, int ib=-1);
  
  void ComputeMagnitudeWithErrors(double *mag, double *err_mag, int ib=-1);
  
  void LoadCFHTLSErrorMag();
  void LoadCFHTLSErrorMag(string filename);
  double CFHTLS_errormag(double mag, int i);

  // return the apparent magnitude affect by an error, which comes from the flux error
  // that is gaussian
  // this need to be corrected one day, in order to consider the flux as "poissonian"
  // when the flux is faint
  void ComputeApparentMagnitudeWithErrors_FluxGaus(double *magnitude, double *error, int ib=-1);
  
  // return in each band bool=1 if observed (as function of error and maglim)
  void ObservedBand_FluxGaus(double *magnitude, double *error, int *band_observed);
  
  // compute the flux and the uncertainty, for a given value of the apparent
  // magnitude, and its uncertainty
  void ComputeFlux(float *mag, float *err_mag, double *flux, double *err_flux, int i);
  void ComputeFlux(float *mag, float *err_mag, double *flux, double *err_flux);

  void ComputeFlux(float *mag, double *flux);
  void ComputeFlux(double *mag, double *flux);

  float ComputeMag(double flux, int i);
  
  // return TFlu::EvalFromTable(type, redshift, ebv, i);
  double ComputeFluxTemplate(double redshift, int type, double ebv, int i);
  
  // set TDistance object
  void SetDistance(TDistance *distance);
  
  // set TKcorrection object
  void SetKcorrectionTable(TKcorrection *table);
  
  // set TFlux object
  void SetFluxTable(TFlux *table);
  
  int GetNBands();
  
  // set number of visit
  void SetNvisit(double *nvis, int n);
  
  bool ParametersInRange(float z, int type, float ebv);
  
  void SetParRange(double *vec);
  int GetNTypes();
  double GetEBV(int type);

  // set distance modulus
  double MD(double z);
  
  void SetSurvey(TDataCard *card);
  
  void LoadDataCard(TDataCard *card);
  void SetF0(TDataCard *card);
  void SetF0(TFilter *filter);
  void SetFlim();

  void SetNewFilter();
  TFilter *GetFilters();
  
  double GetFLUXMAGNULL(int i);

  TSed *GetSed();
  
  void Write();
  
  double *fnvis;
  double *fm5;
  double *F0;
  double *Flim;
  double *maglim;
  double *gamma;
  double err_syst;
  int nband;
  string survey;
  
 private:
  TKcorrection *KcorrectionTable;
  TFlux *FluxTable;
  TDistance *distance;
  const double MAGNULL;
  const double ERRMAGNULL;
  double FLUXMAGNULL[6];
  const double celerite;
  double integral_Bfilter;
  TFile *fcfh;
  TH2F *hsigma_cfh[5];
  double zmin, zmax, ebvmin, ebvmax;
  int typemin, typemax;
  double *ebvlimit;
  
  TSed *sed;
  bool noebv;
  bool print;
};
#endif /* TAPPARENTMAGNITUDE_H_ */

