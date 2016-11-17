
#ifndef TIGM_h
#define TIGM_h
//___________________________
/*
 * TIGM.h
 *
 *  Created on: Sep 27, 2010
 *      Author: gorecki
 *
 *      cf Madau  1995/2000
 *
 *      Modified on: Sep 19 2014 by JS Ricol
 *

 */

#include <TROOT.h>
#include <TChain.h>
#include "TFile.h"
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

#include "TDataCard.h"

using namespace std;

class TIGM {

public:
  // constructor
  TIGM(){;}
  // destructor
  ~TIGM();
  void Load(TDataCard *card);
  // Load IGM table that can be made with the prog: $PhotoZ/progs/IGM/MakeIGMTable.exe
  void LoadIGMTable();
  // Set igm table name
  void SetIGMMadau(string filename);
  // Compute IGM extinsion law as a function of lambda and redshift
  void GetIGMTable(double *lambda, double *z, double **absorption);
  // get  number of set in lambda and in redshift for which the
  // table is computed
  void GetParameters(int &n_lambda, int &n_z);
  // eval the igm extinsion law from the table
  double IGMTableEval(double lambda, double z);
  // return igm extinsion at a given emitted wavelength for source at a given redshift
  double IGMExtinction(double lambda /*emis*/, double zsource);
  // return term of equation 16 in Madau's paper
  double sigmaEq16(double lambda_obs, double z);
  // return term from Madau's paper
  double dNH1(double nh1, double z);
  // compute double integral
  double IntegralDouble(double ze, double zs);
 private:
  int fmadau;
  double **fabsorption;
  double *vect_z;
  double *vect_lambda;
  double tabsorption;
  double tlambda;
  double tz;
  int nlambda;
  int nredshift;
  double tlambda_min;
  double tz_min;
  TTree *TreeIGM;
  TFile *fileIGM;
  double lambdamax;
};
#endif /* TIGM_h*/
