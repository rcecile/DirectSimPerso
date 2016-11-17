/*
 * TFlux.C
 *
 *  Created on: Mar 29, 2010
 *      Author: gorecki
 */

#include "TFlux.h"

#include "TKey.h"
#include "TParameter.h"
using namespace std;

TFlux::TFlux() :
  epsilon(1.e-9){
  lambda_min=1.e-7;
  lambda_max=3.e-5; 
}

void TFlux::SetSed(TSed *sed) {
  sed_fromfile = sed;
}

TSed *TFlux::GetSed() {
  return sed_fromfile;
}

void TFlux::SetSedNumber(int sed) {
  ised = sed;
}

void TFlux::SetFilter(TFilter *filter) {
  nband = filter -> nfilter;
  filter_fromfile = filter;

}

void TFlux::SetBFilter(TFilter *filter) {
  bfilter_fromfile = filter;

}

TFlux::~TFlux() {
  if (bfilter_fromfile != NULL)
    delete bfilter_fromfile;
  delete filter_fromfile;
  delete sed_fromfile;
}

void TFlux::EvalRedshiftedSedIntegral(double z, double *Flux) {
  
  for (int i = 0; i < nband; i++) {
    double integral = 0.;
    
    lambda_min = filter_fromfile -> GetLambdaMin(i);
    lambda_max = filter_fromfile -> GetLambdaMax(i);

    for (double lambda = lambda_min; lambda < lambda_max; lambda += epsilon) {
      integral += filter_fromfile -> Eval(lambda, i)
	* sed_fromfile -> Eval(lambda / (1. + z), ised) * lambda
	* epsilon;
    }
    Flux[i] = integral;

  }
}

void TFlux::EvalRedshiftedSedIntegralLambda2(double z, double *Flux) {
  
  for (int i = 0; i < nband; i++) {
    double integral = 0.;
    
    lambda_min = filter_fromfile -> GetLambdaMin(i);
    lambda_max = filter_fromfile -> GetLambdaMax(i);
    
    for (double lambda = lambda_min; lambda < lambda_max; lambda += epsilon) {
      integral += filter_fromfile -> Eval(lambda, i)
	* sed_fromfile -> Eval(lambda / (1. + z), ised) * lambda
	* lambda * epsilon;
    }
    Flux[i] = integral;
  }
}

double TFlux::EvalSedBIntegral() {
  double integral = 0.;

  lambda_min = bfilter_fromfile -> GetLambdaMin(0);
  lambda_max = bfilter_fromfile -> GetLambdaMax(0);
  
  for (double lambda = lambda_min; lambda < lambda_max; lambda += epsilon) {
    integral += bfilter_fromfile -> Eval(lambda, 0) * sed_fromfile -> Eval(lambda, ised) * lambda * epsilon;
  }
  return integral;
}

double TFlux::EvalBIntegral() {
  double integral = 0.;

  lambda_min = bfilter_fromfile -> GetLambdaMin(0);
  lambda_max = bfilter_fromfile -> GetLambdaMax(0);
  
  for (double lambda = lambda_min; lambda < lambda_max; lambda += epsilon) {
    integral += bfilter_fromfile -> Eval(lambda, 0) / lambda * epsilon;
  }
  return integral;
  
}

void TFlux::SetIntegralBFilter(double intergal) {
  integral_Bfilter = intergal;
}

void TFlux::LoadTable(string filename) {

  // Load table of Kcorrection
  // the root file must contain a ttree produces
  // by the program $PhotoZ/progs/Kcorrection/MakeTableKcorrection.cxx
  
  ftable = new TFile(filename.c_str());
  TList* list = ftable -> GetListOfKeys();
  
  TList* UserInfo;
  //	string name = "aaa";
  //	ftable -> ls();
  
  for (int j = 0; j < list->GetEntries(); ++j) {
    TKey* key = (TKey*) list->At(j);
    //cout<< key->GetName()<<endl;
    if (TClass::GetClass(key->GetClassName())->InheritsFrom("TTree")) {
      
      ttable = (TTree*) ftable -> Get("tree");
      
      UserInfo = (TList*) ttable -> GetUserInfo();
      
      zmin = ((TParameter<double>*) UserInfo -> At(0)) -> GetVal();
      zmax = ((TParameter<double>*) UserInfo -> At(1)) -> GetVal();
      zstep = ((TParameter<double>*) UserInfo -> At(2)) -> GetVal();
      nz = ((TParameter<int>*) UserInfo -> At(3)) -> GetVal();
      
      ebvmin = ((TParameter<double>*) UserInfo -> At(4)) -> GetVal();
      ebvmax = ((TParameter<double>*) UserInfo -> At(5)) -> GetVal();
      ebvstep = ((TParameter<double>*) UserInfo -> At(6)) -> GetVal();
      nebv = ((TParameter<int>*) UserInfo -> At(7)) -> GetVal();
      
      ntype = ((TParameter<int>*) UserInfo -> At(8)) -> GetVal();
      typemin=0;
      typemax=ntype;

      nfilter = ((TParameter<int>*) UserInfo -> At(9)) -> GetVal();
      nband = nfilter;
      tbtable = (TTree*) ftable -> Get("treeB");
      //			cout << " sed = " << ntype << endl;
      break;
      
    }
    if (!TClass::GetClass(key->GetClassName())->InheritsFrom("TTree")) {
      cout << " could not open Flux table" << endl;
    }
  }
  
  //cout << "Loading Table Flux ... redshift from " << zmin << " to " << zmax << " step=" << zstep << " n=" << nz << endl;
  //cout << "Loading Table Flux ... E(B-V) from " << ebvmin << " to " << ebvmax << " step=" << ebvstep << " n=" << nebv << endl;
  //cout << "Loading Table Flux ... Types = " << ntype << endl;

  double k[nfilter], fb;
  double ebv, z;
  int type;
  ttable -> SetBranchAddress("f", &k);
  ttable -> SetBranchAddress("z", &z);
  ttable -> SetBranchAddress("ebv", &ebv);
  ttable -> SetBranchAddress("type", &type);
  tbtable -> SetBranchAddress("f", &fb);
  
  fflux = new double***[ntype + 1];
  ffluxb = new double**[ntype + 1];

  az_ = new double***[ntype + 1];
  bz_ = new double***[ntype + 1];
  fa_ = new double***[ntype + 1];
  fb_ = new double***[ntype + 1];
  fc_ = new double***[ntype + 1];
  fd_ = new double***[ntype + 1];
  

  for (int i = 0; i < ntype + 1; i++) {
    fflux[i] = new double**[nz + 1];
    az_[i] = new double**[nz + 1];
    bz_[i] = new double**[nz + 1];
    fa_[i] = new double**[nz + 1];
    fb_[i] = new double**[nz + 1];
    fc_[i] = new double**[nz + 1];
    fd_[i] = new double**[nz + 1];
    ffluxb[i] = new double*[nz + 1];
    for (int j = 0; j < nz + 1; j++) {
      fflux[i][j] = new double*[nebv + 1];
      az_[i][j] = new double*[nebv + 1];
      bz_[i][j] = new double*[nebv + 1];
      fa_[i][j] = new double*[nebv + 1];
      fb_[i][j] = new double*[nebv + 1];
      fc_[i][j] = new double*[nebv + 1];
      fd_[i][j] = new double*[nebv + 1];
      ffluxb[i][j] = new double[nebv + 1];
      for (int k = 0; k < nebv + 1; k++) {
	fflux[i][j][k] = new double[nfilter];
	az_[i][j][k] = new double[nfilter];
	bz_[i][j][k] = new double[nfilter];
	fa_[i][j][k] = new double[nfilter];
	fb_[i][j][k] = new double[nfilter];
	fc_[i][j][k] = new double[nfilter];
	fd_[i][j][k] = new double[nfilter];
      }
    }
  }
  
 
  /***********************
   * recupere les valeurs
   ***********************/
  
  int iebv, itype, iz;

  //cout << ttable->GetEntries() << endl;

  for (int ientry=0; ientry<ttable->GetEntries(); ientry++){
      ttable -> GetEntry(ientry);
      iebv = int((ebv-ebvmin)/ebvstep+ebvstep/10.);
      iz = int((z-zmin)/zstep+zstep/10.);
      itype = type;
      
      tbtable -> GetEntry(ientry);
      
      for (int i_band = 0; i_band < nfilter; i_band++) {
	if (iebv < nebv && iz < nz && itype < ntype){
	  fflux[itype][iz][iebv][i_band] = k[i_band];
	  ffluxb[itype][iz][iebv] = fb;
	}
      }
  }
  
  //cout << "done" << endl;
  
  for (int i_band = 0; i_band < nfilter; i_band++){
    iebv = nebv;
    for (itype = 0; itype < ntype; itype++) {
      for (iz = 0; iz < nz; iz++) {
	//cout << iz << " " << itype << " " << iebv << endl;
	fflux[itype][iz][iebv][i_band] = fflux[itype][iz][iebv-1][i_band];
	ffluxb[itype][iz][iebv] = ffluxb[itype][iz][iebv-1];
      }
    }
    iz = nz;
    for (itype = 0; itype < ntype; itype++) {
      for (iebv = 0; iebv < nebv; iebv++) {
	fflux[itype][iz][iebv][i_band] = fflux[itype][iz-1][iebv][i_band];
	ffluxb[itype][iz][iebv] = ffluxb[itype][iz-1][iebv];
      }
    }
    itype = ntype;
    for (iz = 0; iz < nz; iz++) {
      for (iebv = 0; iebv < nebv; iebv++) {
	fflux[itype][iz][iebv][i_band] = fflux[itype-1][iz][iebv][i_band];
	ffluxb[itype][iz][iebv] = ffluxb[itype-1][iz][iebv];
      }
    }
  }

  for (int iband = 0; iband < nfilter; iband++)
    for (itype = 0; itype < ntype; itype++)
      for (iz = 0; iz < nz; iz++){
	double x1 = double(iz)*zstep;
	double x2 = double(iz+1)*zstep;
	for (iebv = 0; iebv < nebv; iebv++) {
	  double y1 = double(iebv)*ebvstep;
	  double y2 = double(iebv+1)*ebvstep;
	  double Q11 = fflux[itype][iz][iebv][iband];
	  double Q21=fflux[itype][iz+1][iebv][iband];
	  double Q12 = fflux[itype][iz][iebv+1][iband];
	  double Q22=fflux[itype][iz+1][iebv+1][iband];
	  fa_[itype][iz][iebv][iband] = (-Q11*y2+Q21*y2+Q12*y1-Q22*y1)/(zstep*ebvstep);
	  fb_[itype][iz][iebv][iband] = (-Q11*x2+Q21*x1+Q12*x2-Q22*x1)/(zstep*ebvstep);
	  fc_[itype][iz][iebv][iband] = (Q11-Q21-Q12+Q22)/(zstep*ebvstep);
	  fd_[itype][iz][iebv][iband] = (Q11*x2*y2-Q21*x1*y2-Q12*x2*y1+Q22*x1*y1)/(zstep*ebvstep);
	  az_[itype][iz][iebv][iband] = (Q21-Q11)/zstep;
	  bz_[itype][iz][iebv][iband] = Q11-x1*(Q21-Q11)/zstep;
	}
      }
}

double *TFlux::GetParRange(){
  double *vec = new double[6];
  int i=0;
  vec[i++] = zmin;
  vec[i++] = zmax;
  vec[i++] = typemin;
  vec[i++] = typemax;
  vec[i++] = ebvmin;
  vec[i++] = ebvmax;
  return vec;
}

bool TFlux::ParametersNotInRange(double z, int type, double ebv){
  bool ok=1;
  if (z<zmin || z>zmax){
    cout << "WARNING z= " << z << " not in range of FluxTable [" << zmin << ", " << zmax << "]" << endl;
    ok=0;
  }
  if (type<typemin || type>typemax){
    cout << "WARNING type= " << type << " not in range of FluxTable [" << typemin << ", " << typemax << "]" << endl;
    ok=0;
  }
  if (ebv<ebvmin || ebv>ebvmax){
    cout << "WARNING ebv= " << ebv << " not in range of FluxTable [" << ebvmin << ", " << ebvmax << "]" << endl;
    ok=0;
  }
  return ok;
}



double TFlux::EvalFromTable(double type, double redshift, double ebv,
			    int filter) {
  int itype=int(type+0.5);
  return EvalFromTable(itype, redshift, ebv, filter);
}


double TFlux::EvalFromTable(int type, double redshift, double ebv,
			    int filter) {
  return linearinterpolation(type, redshift, ebv, filter);
  //return bilinearinterpolation(type, redshift, ebv, filter);
  //return fastinterpolation(type, redshift, ebv, filter);
}

double TFlux::linearinterpolation(int type, double redshift, double ebv,
				  int filter) {
  int iz = int(redshift/zstep+0.0001);
  int iebv = int(ebv/ebvstep+0.0001);
  
  return az_[type][iz][iebv][filter]*redshift+bz_[type][iz][iebv][filter];
  
}
double TFlux::fastinterpolation(int type, double redshift, double ebv,
				    int filter) {
  int iz = int(redshift/zstep+0.000001);
  int iebv = int(ebv/ebvstep+0.000001);

  return fa_[type][iz][iebv][filter]*redshift+fb_[type][iz][iebv][filter]*ebv+fc_[type][iz][iebv][filter]*redshift*ebv+fd_[type][iz][iebv][filter];
}

double TFlux::bilinearinterpolation(int type, double redshift, double ebv,
				      int filter) {
  /*If the min lambda of "B filter" is lower than
   * the lyman absorption, the integral of the SED in the "B filter"
   * is not redshift dependent
   *
   * */
  
  // cout << "flux eval : ebv = " << ebv << endl;

  double return_value = 0.;
  int i_filter = filter;
  
  int i_temp = type;
  int i_z = int(redshift/zstep+0.000001);
  int i_ebv = int(ebv/ebvstep+0.000001);
  
  if (i_z>=nz){
    cout << "WARNING z not in range of FluxTable : " << redshift << endl;
    i_z = nz-1;
  }
  
  if (i_ebv>=nebv){
    cout << "WARNING EBV not in range of FluxTable : " << ebv << endl;
    i_ebv = nebv-1;
  }

  if (i_temp>=ntype){
    cout << "WARNING Type not in range of FluxTable : " << i_temp << endl;
    i_temp = ntype-1;
  }
  
  if (ebv>0){
    double x1 = double(i_z)*zstep;
    double x2 = double(i_z+1)*zstep;
    double y1 = double(i_ebv)*ebvstep;
    double y2 = double(i_ebv+1)*ebvstep;
    double x = redshift;
    double y = ebv;
    int i_zp1 = i_z+1;
    int i_ebvp1 = i_ebv+1;

    if (i_ebv==nebv-1)
      i_ebvp1=i_ebv;
    if (i_z==nz-1)
      i_zp1=i_z;
    
    double fQ11=fflux[i_temp][i_z][i_ebv][i_filter];
    double fQ21=fflux[i_temp][i_zp1][i_ebv][i_filter];

    double fR1=(x2-x)/(x2-x1)*fQ11+(x-x1)/(x2-x1)*fQ21;
    double fQ12=fflux[i_temp][i_z][i_ebvp1][i_filter];
    double fQ22=fflux[i_temp][i_zp1][i_ebvp1][i_filter];
    double fR2 = (x2-x)/(x2-x1)*fQ12+(x-x1)/(x2-x1)*fQ22;
    
    return_value = (y2-y)/(y2-y1)*fR1+(y-y1)/(y2-y1)*fR2;
    
    return return_value;
  }
  else{
    if (i_z==nz-1)
      return fflux[i_temp][i_z][i_ebv][i_filter];
    int i_zp1 = i_z+1;
    double x1 = double(i_z)*zstep;
    double x2 = double(i_zp1)*zstep;
    double y1 = fflux[i_temp][i_z][i_ebv][i_filter];
    double y2 = fflux[i_temp][i_zp1][i_ebv][i_filter];
    double x = redshift;
    double y = y1 + (y2-y1)/(x2-x1)*(x-x1);
    return y;
  }
}



double TFlux::EvalIntegralSedBFromTable(int type, double redshift, double ebv) {
  /*If the min lambda of "B filter" is lower than
   * the lyman absorption, the integral of the SED in the "B filter"
   * is not redshift dependent
   *
   * */
  
  double return_value = 0.;
  int i_temp = type;
  int i_z = int(redshift/zstep+0.000001);
  int i_ebv = int(ebv/ebvstep+0.000001);
  
  if (i_z>=nz){
    cout << "WARNING z not in range of FluxTable : " << redshift << endl;
    i_z = nz-1;
  }
  
  if (i_ebv>=nebv){
    cout << "WARNING EBV not in range of FluxTable : " << ebv << endl;
    i_ebv = nebv-1;
  }

  if (i_temp>=ntype){
    cout << "WARNING Type not in range of FluxTable : " << i_temp << endl;
    i_temp = ntype-1;
  }
  
  if (ebv>0){
    double x1 = double(i_z)*zstep;
    double x2 = double(i_z+1)*zstep;
    double y1 = double(i_ebv)*ebvstep;
    double y2 = double(i_ebv+1)*ebvstep;
    double x = redshift;
    double y = ebv;
    int i_zp1 = i_z+1;
    int i_ebvp1 = i_ebv+1;
    if (i_ebv==nebv-1)
      i_ebvp1=i_ebv;
    if (i_z==nz-1)
      i_zp1=i_z;
    
    double fQ11=ffluxb[i_temp][i_z][i_ebv];
    double fQ21=ffluxb[i_temp][i_zp1][i_ebv];
    double fR1=(x2-x)/(x2-x1)*fQ11+(x-x1)/(x2-x1)*fQ21;
    double fQ12=ffluxb[i_temp][i_z][i_ebvp1];
    double fQ22=ffluxb[i_temp][i_zp1][i_ebvp1];
    double fR2 = (x2-x)/(x2-x1)*fQ12+(x-x1)/(x2-x1)*fQ22;
    
    return_value = (y2-y)/(y2-y1)*fR1+(y-y1)/(y2-y1)*fR2;
    
    return return_value;
  }
  else{
    if (i_z==nz-1)
      return ffluxb[i_temp][i_z][i_ebv];
    int i_zp1 = i_z+1;
    double x1 = double(i_z)*zstep;
    double x2 = double(i_zp1)*zstep;
    double y1 = ffluxb[i_temp][i_z][i_ebv];
    double y2 = ffluxb[i_temp][i_zp1][i_ebv];
    double x = redshift;
    double y = y1 + (y2-y1)/(x2-x1)*(x-x1);
    return y;
  }
}
