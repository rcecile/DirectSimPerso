

#ifndef TKcorrection_h
#define TKcorrection_h
//________________________________________
/*	TKcorrection.h
 *
 *  Created on:
 *      Author: gorecki
 *
 *
 *
 *      Compute Kcorrections
 *      Bfilter if the filter in which the absolute magnitude
 *      is given
 *      cf Hogg 2002
 *      Note that the definitions of the integral depend on the unit of the
 *      sed which is different from Hogg's paper
 *
 *      ex:
 *      TKcorrection *kcorrection = new TKcorrection();
 *		kcorrection -> SetBFilter(bfilter);
 *		kcorrection -> SetFilter(filter);
 *		kcorrection -> SetSed(sed);
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
#include "TKey.h"
#include "TParameter.h"

#include "TSed.h"
#include "TFilter.h"
using namespace std;

class TKcorrection {

public:
  // constructor
  TKcorrection();
  // destructor
  ~TKcorrection();
  // Initialize
  void Initialize();
  // compute integral of the corresponding filter
  double Integral_filter(int filter);
  // compute integral of the sed in the B filter int B(lambda)Sed(lambda)lambda dlambda
  double Integral_SED_Bfilter();
  // compute integral of the sed at redshift z
  // in the  filter int B(lambda)Sed(lambda/(1+z))lambda dlambda
  double Integral_SEDRedshifted_filter(int filter, double z);
  // integral of the B filter
  double Integral_Bfilter();
  // return the K correction for a given filter and redshit
  double Eval(double redshift, int filter);
  // return the K correction for a given filter and redshit , when
  // the values of Integral_Bfilter,  Integral_filter, Integral_SED_Bfilter
  // are known
  double Eval(double z, int filter, double int_B, double* int_X_all,
	      double int_B_SED);
  // return the value from a table
  double EvalFromTable(int type, double redshift, double ebv,
		       int filter, bool debug=false);
  // rerurn Integral_filter()
  void GetIntegralFilter(double *intergal);
  // set Integral_filter()
  void SetIntegralFilter(double *intergal);
  // Load table of Kcorrection
  void LoadTable(string filename);
  // set Integral_Bfilter()
  void SetIntegralBFilter(double intergal);
  // set Integral_SED_Bfilter()
  void SetIntegralBFilterSED(double integral);
  // set Tsed object
  void SetSed(TSed *sed);
  // get Sed
  TSed* GetSed();
  // set TFilter
  void SetFilter(TFilter *filter);
  // set TFilter  Bfilter object
  void SetBFilter(TFilter *filter);
  // set number of sed
  void SetSedNumber(int sed);

  TFilter *GetFilters();
  
  bool ParametersNotInRange(double z, int type, double ebv);
  double *GetParRange();


  int ised;
  double Int_B;
  double Int_X;
  double *Int_X_all;
  double Int_B_SED;
  int nband;

  void Debug();

 private:
  
  TFilter *filter_fromfile;
  TFilter *bfilter_fromfile;
  TSed *sed_fromfile;
  TFile *ftable;
  TTree *ttable;
  
  const double epsilon;
  
  double lambda_min;
  double lambda_max;
  
  
  double ****ktable;
  int ntype;
  int nfilter;
  int nz;
  int nebv;
  double ebvstep, zstep;
  double zmin, zmax, ebvmin, ebvmax;
  int typemin, typemax;
  
};
#endif /* TKcorrection_h*/
