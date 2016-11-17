
#ifndef TSed_h
#define TSed_h

//___________________________________
//
//  Created on:
//       Author: gorecki
//
//   Make operation of spectral energy distribution SED
//   Extinction law Cardelli and Calzetti are implemented with the class TReddening
//
//  Lambda Unit : nm


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

#include "TReddening.h"
#include "TIGM.h"
#include "TRandom3.h"

using namespace std;

class TSed {

 public:
  // constructor
  TSed(bool p=1);
  //destructor
  ~TSed() {
    delete lambdamin_;
    delete lambdamax_;
    delete gtype;
    delete ebvmax;
  }
  void LoadDataCard(TDataCard *card);
  int GetNtype();
  string *GetTypeArray();
  void SetTypes(TDataCard *card);
  void LoadFile(TDataCard *card);
  void LoadType(TDataCard *card);
  void LoadIGM(TDataCard *card);
  double Eval(double lambda, int ised);
  void SetFileType(string extension);
  void SetIGM(bool igm);
  void LoadIGMTable(TIGM *table);
  void SetReddening(TReddening *reddening);
  void SetEBV(double ebv);
  void SetMaxEBVforType(int type, double ebv);
  double GetEBV(int type);
  double FlatEBV(int type);
  void SetEmittedRedshift(double z);
  void SetMix(bool mix);
  
  void BuildTypeArray();
  int TransformTypeFromCWW(int t);
  
  double GetLambdaMin(int i);
  double GetLambdaMax(int i);

  double lambdamin;
  double* lambdamin_;
  double* lambdamax_;
  double lambdamax;
  int nsed, ntype;
  double febv;
  double figmredshift;
  double mix;
  double reddening;
  double igm;
  double lambdaunit;

  string SedType(int ised);
  string SedExt(int ised);
        
 private:
  TReddening *freddening;
  string fextension;
  ROOT::Math::Interpolator **fsed_fromfile;
  TF1 **fsed_fromfile_tf1;
  TTree **fsed_fromfile_ttree;  
  TFile **f_sed;
  TIGM *figmTable;

  string fextinction;
  string finterpolation;
  bool figm;
  bool fmix;
  int *gtype;
  string *stype;
  string *sext;

  double ebvmax_early, ebvmax_late;
  double ebvmax_sb,  ebvmax_default;
  double ebv_max_out;

  double *ebvmax;
  TRandom3 *rand;

  vector<int> list_early;
  vector<int> list_late;
  vector<int> list_sb;
  vector<int> *vec;
  bool print;
};

#endif /*TSed_h*/
