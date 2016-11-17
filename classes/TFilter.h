
#ifndef TFilter_h
#define TFilter_h

//_________________________
///
// TFilter.h
//
//  Created on: Mar 26, 2010
//      Author: gorecki
//      Modified on: Sep 19 2014 by JS Ricol
//

#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <vector>
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include "TRandom3.h"
#include "TDataCard.h"

using namespace std;

class TFilter {

 public:
  //constructor
  TFilter();
  // destructor
  ~TFilter();
  // load input parameters
  void LoadDataCard(TDataCard *card);
  // return the value of the filter
  double Eval(double lambda, int ifilter);
  // Load sed files : number of files and file names
  void LoadFile(int nfile, string *filename);
  // return the effective value of the filter wavelength
  void Lambda_Effective(double *lambda_eff);
  // lambda is in m. If the wavelength of the files are in nm , then must set conversionLambda_filter = 1e9
  void SetConversionUnitLambda(double conversionLambda_filter);
  // specify file type: either txt or ttree or tf1
  //the last tow are root files
  void SetFileType(string extension);
  string GetFileType();
  // return int F(lambda)/lambda dlambda
  double Integral_Eval_FSurlambdaDlambda(int filter);

  int nfilter;
  double fconversionUnit_Lambda;
  string *ffilename;

  double GetLambdaMin(int i);
  double GetLambdaMax(int i);

  void Draw(string opt);

  bool MultiFilter();
  void NewFilters();

  void SetDLamdaNanometer(double dl);
  void SetDTransmission(double dt);

  double GetNorm(int i);

 private:
  double epsilon;
  double *lambdamin;
  double *lambdamax;
  string fextension;
  string bfextension;
  ROOT::Math::Interpolator **ffilter_fromfile;
  TF1 **ffilter_fromfile_tf1;
  TTree **ffilter_fromfile_ttree;
  TFile **f_filter;

  bool multifilter;

  TRandom3 *rand;

  TList **list;
  TKey *key;
  
  double dT;
  double dlambda;
  double *norm;
};

#endif /*TFilter_h*/
