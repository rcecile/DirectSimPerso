//#define TKcorrection_cxx
#include "TKcorrection.h"

using namespace std;

TKcorrection::TKcorrection() :
  epsilon(1.e-9) {
  lambda_min = 1.e-7;
  lambda_max = 3.e-5;
}

void TKcorrection::SetSed(TSed *sed) {
  sed_fromfile = sed;
}

TSed *TKcorrection::GetSed() {
  return sed_fromfile;
}

void TKcorrection::SetSedNumber(int sed) {
  ised = sed;
}

void TKcorrection::SetFilter(TFilter *filter) {
  nband = filter -> nfilter;
  filter_fromfile = filter;
  Int_X_all = new double[nband];
}

TFilter* TKcorrection::GetFilters() {
  return filter_fromfile;
}

void TKcorrection::SetBFilter(TFilter *filter) {
  bfilter_fromfile = filter;
  Int_B = Integral_Bfilter();

}

TKcorrection::~TKcorrection() {
  delete[] Int_X_all;
  delete filter_fromfile;
  delete sed_fromfile;
}

void TKcorrection::GetIntegralFilter(double *integral) {
  for (int i = 0; i < nband; i++) {
    Int_X_all[i] = Integral_filter(i);
    integral[i] = Int_X_all[i];
  }
}

void TKcorrection::SetIntegralFilter(double *integral) {
  for (int i = 0; i < nband; i++) {
    Int_X_all[i] = integral[i];
  }
}

void TKcorrection::SetIntegralBFilter(double integral) {
  Int_B = integral;
}

void TKcorrection::SetIntegralBFilterSED(double integral) {
  Int_B_SED = integral;
}

double TKcorrection::Integral_SEDRedshifted_filter(int filter, double z) {
  double return_value = 0.;
  sed_fromfile -> SetEmittedRedshift(z);

  lambda_min = filter_fromfile -> GetLambdaMin(filter);
  lambda_max = filter_fromfile -> GetLambdaMax(filter);

  for (double lambda = lambda_min; lambda < lambda_max; lambda += epsilon) {
    double val = filter_fromfile -> Eval(lambda, filter)
      * sed_fromfile -> Eval(lambda / (1 + z), ised) * epsilon
      * lambda;
    return_value += val;
    //cout << "Band " << filter << " l=" << lambda << "; filter=" << filter_fromfile -> Eval(lambda, filter)
    //<< "; sed=" << sed_fromfile -> Eval(lambda / (1 + z), ised) << "; val = " << val << endl;
  }
  //cout << "TOTAL = " << return_value << endl;
  
  return return_value;
}

double TKcorrection::Integral_SED_Bfilter() {
  double return_value = 0.;

  lambda_min = bfilter_fromfile -> GetLambdaMin(0);
  lambda_max = bfilter_fromfile -> GetLambdaMax(0);
  
  for (double lambda = lambda_min; lambda < lambda_max; lambda += epsilon) {
    return_value += bfilter_fromfile -> Eval(lambda, 0)
      * sed_fromfile -> Eval(lambda, ised) * epsilon * lambda;
  }
  
  return return_value;
}

double TKcorrection::Integral_filter(int filter) {
  double return_value = 0.;

  lambda_min = filter_fromfile -> GetLambdaMin(filter);
  lambda_max = filter_fromfile -> GetLambdaMax(filter);
  
  for (double lambda = lambda_min; lambda < lambda_max; lambda += epsilon) {
    return_value += filter_fromfile -> Eval(lambda, filter) * epsilon
      / lambda;
  }
  return return_value;
}

double TKcorrection::Integral_Bfilter() {
  double return_value = 0.;
  
  lambda_min = bfilter_fromfile -> GetLambdaMin(0);
  lambda_max = bfilter_fromfile -> GetLambdaMax(0);

  // cout << "bfilter : " << lambda_min << "; " << lambda_max << endl;
  
  for (double lambda = lambda_min; lambda < lambda_max; lambda += epsilon) {
    return_value += bfilter_fromfile -> Eval(lambda, 0) * epsilon / lambda;
  }
  return return_value;
}

double TKcorrection::Eval(double z, int filter, double int_B,
			  double* int_X_all, double int_B_SED) {
  double K = 999;
  double integral_X_z = Integral_SEDRedshifted_filter(filter, z);
  
  if (integral_X_z != 0) {
    K = -2.5 * log(integral_X_z * int_B / int_X_all[filter] / int_B_SED
		   / (1. + z)) / log(10.);
  }
  return K;
}

double TKcorrection::Eval(double z, int filter) {
  double K = 999.;
  double integral_X_z = Integral_SEDRedshifted_filter(filter, z);
  double int_B = Integral_Bfilter();
  double int_X_all = Integral_filter(filter);
  double int_B_SED = Integral_SED_Bfilter();
  
  if (integral_X_z != 0) {
    K = -2.5 * log(integral_X_z * int_B / int_X_all / int_B_SED / (1. + z)) / log(10.);
  }

  //cout << z << " " << filter << " " << K << " " << integral_X_z << " " << int_B
  //   << " " << int_X_all << " " << int_B_SED << endl;
  
  return K;
}

void TKcorrection::LoadTable(string filename) {
  ftable = new TFile(filename.c_str());
  TList* list = ftable -> GetListOfKeys();
  TList* UserInfo;
  string name = "aaa";

  for (int j = 0; j < list->GetEntries(); ++j) {
    TKey* key = (TKey*) list->At(j);
    if (TClass::GetClass(key->GetClassName())->InheritsFrom("TTree")) {
      name = key->GetName();
      ttable = (TTree*) ftable -> Get(name.c_str());
      
      UserInfo = (TList*) ttable -> GetUserInfo();
      if (UserInfo -> At(0)!=NULL){
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
      }


      nband = nfilter;
    } else {
      cout << " could not open Kcorrection table" << endl;
    }
  }
  
  /*
  cout << "Loading Table Kcorrection ... redshift from "
       << zmin << " to " << zmax << " step=" << zstep << " n=" << nz << endl;
  cout << "Loading Table Kcorrection ... E(B-V) from "
       << ebvmin << " to " << ebvmax << " step=" << ebvstep << " n=" << nebv << endl;
  cout << "Loading Table Kcorrection ... Types = "
       << ntype << endl;
  */
  
  double k[nfilter];
  double z;
  double ebv;
  int type;
  ttable -> SetBranchAddress("k", &k);
  ttable -> SetBranchAddress("z", &z);
  ttable -> SetBranchAddress("ebv", &ebv);
  ttable -> SetBranchAddress("type", &type);
  
  //cout << "table initialization" << endl;
  
  ktable = new double***[ntype + 1];
  for (int i = 0; i < ntype+1; i++) {
    ktable[i] = new double**[nz + 1];
    for (int j = 0; j < nz + 1; j++) {
      ktable[i][j] = new double*[nebv + 1];
      for (int k = 0; k < nebv + 1; k++)
	ktable[i][j][k] = new double[nfilter];
    }
  }
  
  /***********************
   * recupere les valeurs
   ***********************/
  
  int iebv, itype, iz;
  
  //cout << "filling table" << endl;
  
  for (int ientry=0; ientry<ttable->GetEntries(); ientry++)
    {
      ttable -> GetEntry(ientry);
      iebv = int((ebv-ebvmin)/ebvstep+ebvstep/10.);
      iz = int((z-zmin)/zstep+zstep/10.);
      itype = type;
      
      for (int i_band = 0; i_band < nfilter; i_band++) {
	if (iebv < nebv && iz < nz && itype < ntype)
	  ktable[itype][iz][iebv][i_band] = k[i_band];
	if (iebv == nebv)
	  ktable[itype][iz][iebv][i_band] = ktable[itype][iz][iebv - 1][i_band];
	if (itype == ntype)
	  ktable[itype][iz][iebv][i_band] = ktable[itype - 1][iz][iebv][i_band];
	if (iz == nz)
	  ktable[itype][iz][iebv][i_band] = ktable[itype][iz - 1][iebv][i_band];
      }
    }
  
  for (int iebv=0; iebv<=nebv; iebv++)
    for (int iz=0; iz<=nz; iz++)
      for (int itype=0; itype<=ntype; itype++)
	{
	  if(iebv == nebv)
	  ktable[itype][iz][iebv] = ktable[itype][iz][iebv - 1];
	  if (itype == ntype)
	    ktable[itype][iz][iebv] = ktable[itype - 1][iz][iebv];
	  if (iz == nz)
	    ktable[itype][iz][iebv] = ktable[itype][iz - 1][iebv];
	}
  
  //Debug();

}

void TKcorrection::Debug(){
  for (int i_ebv=0; i_ebv<nebv; i_ebv++)
    for (int i_z=0; i_z<nz; i_z++)
      for (int i_temp=0; i_temp<ntype; i_temp++)
	cout << i_ebv << " " << i_z << " " << i_temp << " " << ktable[i_temp][i_z][i_ebv][3] << endl;;
}

double *TKcorrection::GetParRange(){
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

bool TKcorrection::ParametersNotInRange(double z, int type, double ebv){
  bool ok=1;
  if (z<zmin || z>zmax){
    cout << "WARNING z= " << z << " not in range of KCorrectionTable [" << zmin << ", " << zmax << "]" << endl;
    ok=0;
  }
  if (type<typemin || type>typemax){
    cout << "WARNING type= " << type << " not in range of KCorrectionTable [" << typemin << ", " << typemax << "]" << endl;
    ok=0;
  }
  if (ebv<ebvmin || ebv>ebvmax){
    cout << "WARNING ebv= " << ebv << " not in range of KCorrectionTable [" << ebvmin << ", " << ebvmax << "]" << endl;
    ok=0;
  }
  return ok;
}


double TKcorrection::EvalFromTable(int type, double redshift, double ebv,
				   int filter, bool debug) {

  // Eval Kcorrection from table

  double return_value = 0.;
  int i_filter =  filter;
  
  int i_temp = type;
  int i_z = int(redshift/zstep+0.000001);
  int i_ebv = int(ebv/ebvstep+0.000001);
  

  // cout << i_z << " " << i_temp << " " << i_ebv << endl;


  if (i_z>=nz){
    
    i_z = nz-1;
  }
  
  if (i_ebv>=nebv){
    cout << "WARNING EBV not in range of KCorrectionTable : " << ebv << endl;
    i_ebv = nebv-1;
  }

  if (i_temp>=ntype){
    cout << "WARNING Type not in range of KCorrectionTable : " << i_temp << endl;
    i_temp = ntype-1;
  }
  
  
  if (ebv>0){
    double x1 = double(i_z)*zstep;
    double x2 = double(i_z+1)*zstep;
    double y1 = double(i_ebv)*ebvstep;
    double y2 = double(i_ebv+1)*ebvstep;
    double x = redshift;
    double y = ebv;
    double fQ11, fQ21, fQ12, fQ22;
    double fR1, fR2;
    
    int i_zp1 = i_z+1;
    int i_ebvp1 = i_ebv+1;
    
    if (i_ebv==nebv-1)
      i_ebvp1=i_ebv;
    if (i_z==nz-1)
      i_zp1=i_z;
    
    fQ11=ktable[i_temp][i_z][i_ebv][i_filter];
    fQ21=ktable[i_temp][i_zp1][i_ebv][i_filter];
    fR1=(x2-x)/(x2-x1)*fQ11+(x-x1)/(x2-x1)*fQ21;
    
    fQ12=ktable[i_temp][i_z][i_ebvp1][i_filter];
    fQ22=ktable[i_temp][i_zp1][i_ebvp1][i_filter];
    fR2 = (x2-x)/(x2-x1)*fQ12+(x-x1)/(x2-x1)*fQ22;
    
    return_value = (y2-y)/(y2-y1)*fR1+(y-y1)/(y2-y1)*fR2;
    return return_value;
  }
  else{
    if (i_z==nz-1)
      return ktable[i_temp][i_z][i_ebv][i_filter];
    int i_zp1 = i_z+1;
    double z1 = double(i_z)*zstep;
    double z2 = double(i_zp1)*zstep;
    double y1 = ktable[i_temp][i_z][i_ebv][i_filter];
    double y2 = ktable[i_temp][i_zp1][i_ebv][i_filter];
    double y = y1+(y2-y1)/(z2-z1)*(redshift-z1);
    return y;
  }

}
