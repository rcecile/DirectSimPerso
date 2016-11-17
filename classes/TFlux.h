#ifndef TFlux_h
#define TFlux_h

//________________________________________
/*	TKcorrection.h
 *
 *  Created on:
 *      Author: gorecki
 *
 *
 *
 *      Compute flux
 *      ex:
 *      	TFlux *flux = new TFlux();
 *			flux -> SetBFilter(bfilter);
 *			flux -> SetFilter(filter);
 * 			flux -> SetSed(sed);
 */

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
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
#include "TSed.h"
#include "TFilter.h"
#include "TGraph.h"

using namespace std;

class TFlux {

public:
  //constructor
  TFlux();
  //destructor
  ~TFlux();
  // Initialize
  void Initialize();
  // return integral sed(lambda/(1+z))X(lambda)lambda dlambda
  void EvalRedshiftedSedIntegral(double z, double *Flux);
  // return int sed(lambda/(1+z))X(lambda)lambda^2 dlambda
  void EvalRedshiftedSedIntegralLambda2(double z, double *Flux);
  //return int B(lambda)sed(lambda)lambda dlambda
  double EvalSedBIntegral();
  // return the value from a table
  double EvalFromTable(int type, double redshift, double ebv, int filter);
  double linearinterpolation(int type, double redshift, double ebv, int filter);
  double bilinearinterpolation(int type, double redshift, double ebv, int filter);
  double fastinterpolation(int type, double redshift, double ebv, int filter);
  // return the value from a table
  double EvalFromTable(double type, double redshift, double ebv, int filter);
  // return the value from a table
  double EvalIntegralSedBFromTable(int type, double redshift, double ebv);
  // Load table of TFlux
  void LoadTable(string filename);
  // set Tsed object
  void SetSed(TSed *sed);
  // set number of sed
  void SetSedNumber(int ised);
  // set TFilter
  void SetFilter(TFilter *filter);
  // set TFilter  Bfilter object
  void SetBFilter(TFilter *filter);
  // set Integral_Bfilter()
  void SetIntegralBFilter(double intergal);
  // return int B(lambda)/lambda dlambda
  double EvalBIntegral();
  
  bool ParametersNotInRange(double z, int type, double ebv);
  double *GetParRange();

  TSed *GetSed();

  int ised;
  int nband;
  TTree *ttable;
  TTree *tbtable;
  int ntype;
  int nfilter;
  int nz;
  int nebv;
 private:
  
  TFilter *filter_fromfile;
  TFilter *bfilter_fromfile;
  TSed *sed_fromfile;
  
  const double epsilon;
  double lambda_min;
  double lambda_max;
  
  double integral_Bfilter;
  
  TFile *ftable;
  
  double ****az_;
  double ****bz_;
  double ****fa_;
  double ****fb_;
  double ****fc_;
  double ****fd_;
  double ****fflux;
  double ***ffluxb;
  double ebvstep, zstep;
  double zmin, zmax, ebvmin, ebvmax;
  int typemin, typemax;
  
};

#endif /*TFlux_h*/
