/*
 * TIGM.C
 *
 *  Created on: Sep 27, 2010
 *      Author: gorecki
 */

#include "TIGM.h"
#include "TKey.h"
#include "TList.h"
#include <iostream>

TIGM::~TIGM() {

}


void TIGM::Load(TDataCard *datacard){
  SetIGMMadau(datacard->igmtable);
  LoadIGMTable();
}

void TIGM::SetIGMMadau(string filename) {
  // Nead to specify the file name of the IGM table
  // Must be a root file with a tree with branch
  // lambda IGMabsorption and z
  fileIGM = new TFile(filename.c_str());
  TList* list = fileIGM -> GetListOfKeys();
  string name = "aaa";
  for (int j = 0; j < list->GetEntries(); ++j) {
    TKey* key = (TKey*) list->At(j);
    if (TClass::GetClass(key->GetClassName())->InheritsFrom("TTree")) {
      name = key->GetName();
    } else {
      cout << " could  not open IGM file" << endl;
    }
  }
  TreeIGM = (TTree*) fileIGM -> Get(name.c_str());
  TreeIGM -> SetBranchAddress("z", &tz);
  TreeIGM -> SetBranchAddress("lambda", &tlambda);
  TreeIGM -> SetBranchAddress("IGMabsorption", &tabsorption);
  nlambda = 1580;
  nredshift = 1201;
  tlambda_min = 20.e-9;
  tz_min = 0.;
}

void TIGM::LoadIGMTable() {
  // Load the table
  fabsorption = new double*[nredshift];
  for (int j = 0; j < nredshift; j++) {
    fabsorption[j] = new double[nlambda];
  }
  
  vect_z = new double[nredshift];
  vect_lambda = new double[nlambda];
  
  /***********************
   * recupere les valeurs
   ***********************/
  
  int indx_lambda = 0;
  int indx_z = 0;
  
  lambdamax=0;
  
  for (int jentries = 0; jentries < TreeIGM -> GetEntries(); jentries++) {
    TreeIGM -> GetEntry(jentries);
    fabsorption[indx_z][indx_lambda] = tabsorption;
    if (indx_z == 0) {
      if (tlambda>lambdamax)
	lambdamax=tlambda;
      vect_lambda[indx_lambda] = tlambda;
    }
    indx_lambda++;
    if (indx_lambda == nlambda) {
      vect_z[indx_z] = tz;
      indx_z++;
      indx_lambda = 0;
    }
    
  }
}

void TIGM::GetIGMTable(double *lambda, double *z, double **absorption) {
  for (int i = 0; i < nredshift; i++) {
    for (int j = 0; j < nlambda; j++) {
      absorption[i][j] = fabsorption[i][j];
      if (i == 0) {
	lambda[j] = vect_lambda[j];
      }
    }
    z[i] = vect_z[i];
  }
}

void TIGM::GetParameters(int &n_lambda, int &n_z) {
  n_lambda = nlambda;
  n_z = nredshift;
}

double TIGM::IGMTableEval(double lambda, double z) {
  // return the value of Exp(-tau) from the table
  // for a given value of lambda_emis and z
  double return_value = 0;
  int i_lambda = 0, i_z = 0;
  
  if (lambda>lambdamax)
    return 1;

  if (lambda == tlambda_min)
    i_lambda = 0;
  else {
    while (vect_lambda[i_lambda] < lambda) {
      i_lambda++;
    }
    i_lambda--;
  }

  if (i_lambda>=nlambda-2)
    i_lambda=nlambda-2;
  
  if (z == tz_min)
    i_z = 0;
  else {
    while (vect_z[i_z] < z)
      i_z++;
    i_z--;
  }

  // JS : wtf ?
  if (i_lambda < tlambda_min)
    i_lambda = 0;
  
  double a_z_inf, b_z_inf;
  double a_z_sup, b_z_sup;
  double a_lambda, b_lambda;
  double f_lambda_inf, f_lambda_sup;
  /****************************
   * interpolation sur redshift
   ****************************/
  a_z_inf = (fabsorption[i_z][i_lambda] - fabsorption[i_z + 1][i_lambda])
    / (vect_z[i_z] - vect_z[i_z + 1]);
  b_z_inf = fabsorption[i_z][i_lambda] - a_z_inf * vect_z[i_z];
  f_lambda_inf = a_z_inf * z + b_z_inf;
  
  a_z_sup = (fabsorption[i_z][i_lambda + 1]
	     - fabsorption[i_z + 1][i_lambda + 1]) / (vect_z[i_z] - vect_z[i_z + 1]);
  b_z_sup = fabsorption[i_z][i_lambda + 1] - a_z_sup * vect_z[i_z];
  f_lambda_sup = a_z_sup * z + b_z_sup;
  
  /****************
   * lambda
   *******************************/
  a_lambda = (f_lambda_inf - f_lambda_sup) / (vect_lambda[i_lambda]
					      - vect_lambda[i_lambda + 1]);
  b_lambda = f_lambda_inf - a_lambda * vect_lambda[i_lambda];
  return_value = a_lambda * lambda + b_lambda;
  
  return return_value;
}

double TIGM::dNH1(double nh1, double z) {
  double return_value = 0;
  double NH1inter = 1.59e17;
  if (nh1 <= NH1inter)
    return_value = 2.4e7 * pow(nh1, -1.5) * pow(1 + z, 2.46); // eq 10 of [85]
  if (nh1 >= NH1inter)
    return_value = 1.9e8 * pow(nh1, -1.5) * pow(1 + z, 0.68);
  return return_value;
}

double TIGM::sigmaEq16(double lambda_obs, double z) // sigma in eq 16 of [85]
{
  double return_value = 0;
  double lambda_L = 91.2;
  double m2nm = 1.e-9;
  return_value = 6.3e-18 * pow(lambda_obs / (lambda_L * m2nm), 3) * pow(1 + z, -3);
  //	cout << " lambda_obs = " << lambda_obs << "  lambda_L = " << lambda_L
  //			* m2nm << " rapport lambda = " << pow(lambda_obs
  //			/ (lambda_L * m2nm), 3) << " A+Z = " << pow(1 + z, -3)
  //			<< " sigma = " << return_value << endl;
  return return_value;
}

double TIGM::IntegralDouble(double ze, double zs) {
  double return_value = 0;
  double xe = 1 + ze;
  double xc = 1 + zs;
  double term1 = 0.25 * pow(xc, 3) * (pow(xe, 0.46) - pow(xc, 0.46));
  double term2 = 9.4 * pow(xc, 1.5) * (pow(xe, 0.18) - pow(xc, 0.18));
  double term3 = -0.7 * pow(xc, 3) * (pow(xc, -1.32) - pow(xe, -1.32));
  double term4 = -0.023 * (pow(xe, 1.68) - pow(xc, 1.68));
  return_value = term1 + term2 + term3 + term4;
  return return_value;
}

double TIGM::IGMExtinction(double lambda /*emis*/, double zsource) {
  /*
   * IGM extinction for a source at redshift zsource, and wavelength lambda
   * */
  double lambda_obs = lambda * (1 + zsource);
  double m2nm = 1.e-9;
  double return_value;
  double tau = 0;
  int cas = 0;
  double zmin;
  double lambda_forest[] = {121.6, 102.6, 97.3, 95.0};
  double lambda_L = 91.2;
  double A[] = { 3.6e-3, 1.7e-3, 1.2e-3, 9.3e-4 };
  if (lambda_obs > lambda_forest[0] * m2nm)
    return_value = 1;
  if (lambda_L * m2nm <= lambda && lambda <= lambda_forest[0] * m2nm) {// eq 15 of [85]
    if (lambda < lambda_forest[3] * m2nm) {
      for (int i = 0; i < 4; i++) {
	tau += A[i] * pow(lambda_obs / (lambda_forest[i] * m2nm), 3.46);
      }
    }
    else {
      for (int i = 0; i < 3; i++) {
	if (lambda_forest[i + 1] * m2nm <= lambda
	    && lambda <= lambda_forest[i] * m2nm)
	  for (int j = 0; j < i + 1; j++) {
	    tau += A[j] * pow(lambda_obs / (lambda_forest[j] * m2nm), 3.46);
	  }
      }
    }
    
    return_value = TMath::Exp(-tau);
    cas = 1;
  }
  double templow, temphigh, taulow = tau, tauhigh = tau, zlow, zhigh;
  double lambdalow, lambdahigh;
  lambdalow = (lambda - 1.e-11) * (1 + zsource);
  lambdahigh = (lambda - 1.e-11) * (1 + zsource);
  
  if (lambda <= lambda_L * m2nm) { //eq 16
    for (int i = 0; i < 4; i++) {
      taulow += A[i] * pow(lambdalow / (lambda_forest[i] * m2nm), 3.46);
      tauhigh += A[i] * pow(lambdahigh / (lambda_forest[i] * m2nm), 3.46);
      tau += A[i] * pow(lambda_obs / (lambda_forest[i] * m2nm), 3.46);
    }
    double pasz = 0.01;
    zmin = lambda_obs / (lambda_L * m2nm) - 1;
    zlow = lambdalow / (lambda_L * m2nm) - 1;
    zhigh = lambdahigh / (lambda_L * m2nm) - 1;
    if (zmin < 0)
      zmin = 0.;
    
    if (zlow < 0)
      zlow = 0.;
    if (zhigh < 0)
      zhigh = 0.;
    /*
     * Integral sur z
     * */
    int nz = 100;
    pasz = (zsource - zmin) / nz;
    /*Integral double
     */
    
    tau += IntegralDouble(zsource, zmin);
    taulow += IntegralDouble(zsource, zlow);
    tauhigh += IntegralDouble(zsource, zhigh);
    if ((taulow - tau) <= 0)
      return_value = 1.e-11;
    else
      return_value = TMath::Exp(-tau);
    cas = 2;
    
  }
  return return_value;
}
