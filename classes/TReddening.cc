/*
 * TReddening.C
 *
 *  Created on: 3 aožt 2011
 *      Author: alexiagorecki
 */

#include "TReddening.h"

TReddening::TReddening() {
  // TODO Auto-generated constructor stub
  SetRVCardelli();
  SetRVCalzetti();

}

TReddening::~TReddening() {
  // TODO Auto-generated destructor stub
}

void TReddening::ReadTabCalz1(){
  double lambdaf, extf=0;
  
  string FICHIER =  "/sps/lsst/data/LePhare/lephare_dev/ext/SB_calzetti.dat";

  FILE *file = fopen(FICHIER.c_str(),"r");
  if(file==NULL) {
    printf("impossible d'ouvrir %s en mode écriture\n",FICHIER.c_str());
    exit(1);
  } 
  while(fscanf(file,"%lf %lf",&lambdaf, &extf)!=EOF) {
    tab_lambda_calz1.push_back(lambdaf);
    tab_ext_calz1.push_back(extf);
  }
  fclose(file);
  // cout << "calzetti tab : " << tab_ext_calz1.size() << endl;
}

void TReddening::ReadTabCalz2(){
  double lambdaf, extf=0;
  
  string FICHIER =  "/sps/lsst/data/LePhare/lephare_dev/ext/SB_calzetti_bump1.dat";

  FILE *file = fopen(FICHIER.c_str(),"r");
  if(file==NULL) {
    printf("impossible d'ouvrir %s en mode écriture\n",FICHIER.c_str());
    exit(1);
  } 
  while(fscanf(file,"%lf %lf",&lambdaf, &extf)!=EOF) {
    tab_lambda_calz2.push_back(lambdaf);
    tab_ext_calz2.push_back(extf);
  }
  fclose(file);
  //cout << "calzetti tab : " << tab_ext_calz2.size() << endl;
}

void TReddening::ReadTabCalz3(){
  double lambdaf, extf=0;
  
  string FICHIER =  "/sps/lsst/data/LePhare/lephare_dev/ext/SB_calzetti_bump2.dat";

  FILE *file = fopen(FICHIER.c_str(),"r");
  if(file==NULL) {
    printf("impossible d'ouvrir %s en mode écriture\n",FICHIER.c_str());
    exit(1);
  } 
  while(fscanf(file,"%lf %lf",&lambdaf, &extf)!=EOF) {
    tab_lambda_calz3.push_back(lambdaf);
    tab_ext_calz3.push_back(extf);
  }
  fclose(file);
  //cout << "calzetti tab : " << tab_ext_calz3.size() << endl;
}


void TReddening::ReadTabPrev(){
  double lambdaf, extf=0;
  string FICHIER =  "/sps/lsst/data/LePhare/lephare_dev/ext/SMC_prevot.dat";
  FILE *file = fopen(FICHIER.c_str(),"r");
  if(file==NULL) {
    printf("impossible d'ouvrir %s en mode écriture\n",FICHIER.c_str());
    exit(1);
  } 
  while(fscanf(file,"%lf %lf",&lambdaf, &extf)!=EOF) {
    tab_lambda_prev.push_back(lambdaf);
    tab_ext_prev.push_back(extf);
  }
  fclose(file);
  //cout << "prevot tab : " << tab_ext_prev.size() << endl;
}


double TReddening::a_FarUV(double *xx, double *par) {
  double x = 1. / xx[0];
  double y = x + par[4];
  double return_value = par[0] + par[1]*y + par[2]*y*y + par[3]*y*y*y;
  return return_value;
}

double TReddening::b_FarUV(double *xx, double *par) {
  double x = 1. / xx[0];
  double y = x + par[4];
  double return_value = par[0] + par[1]*y + par[2]*y*y + par[3]*y*y*y;
  return return_value;
}

double TReddening::cardelli_highl(double *lambda, double *par) {
  double x = 1. / lambda[0];
  //double y = x-1.82;
  double return_value = 0.;
  double a_card = par[0] * pow(x, par[1]);
  double b_card = par[2] * pow(x, par[1]);
  return_value = a_card + b_card / par[3];
  return return_value;
}

double TReddening::a_lowl(double *lambda, double *par) {
  double x = 1. / lambda[0];
  double y = x - 1.82;
  double return_value = 0.;
  return_value = par[0] + par[1] * y + par[2] * pow(y, 2) + par[3]
    * pow(y, 3) + par[4] * pow(y, 4) + par[5] * pow(y, 5) + par[6]
    * pow(y, 6) + par[7] * pow(y, 7);
  return return_value;
}

double TReddening::b_lowl(double *lambda, double *par) {
  double x = 1. / lambda[0];
  double y = x - 1.82;
  double return_value = 0.;
  return_value = par[0] * y + par[1] * pow(y, 2) + par[2] * pow(y, 3)
    + par[3] * pow(y, 4) + par[4] * pow(y, 5) + par[5] * pow(y, 6)
    + par[6] * pow(y, 7);
  return return_value;
}

double TReddening::a_lowVL(double *lambda, double *par) {
  
  double x = 1. / lambda[0];
  double y = x + par[6];
  double FA = 0.;
  if (x >= 5.9 && x <= 8)
    FA = par[5] * pow(y, 2) + par[7];
  else
    FA = 0.;
  double return_value = 0.;
  return_value = par[0] + par[1] * x + par[2] / (pow(x + par[3], 2) + par[4])
    + FA;
  return return_value;
}

double TReddening::b_lowVL(double *lambda, double *par) {
  double x = 1. / lambda[0];
  double y = x + par[6];
  double FB = 0.;
  if (x >= 5.9 && x <= 8)
    FB = par[5] * pow(y, 2) + par[7];
  else
    FB = 0.;
  double return_value = 0.;
  return_value = par[0] + par[1] * x + par[2] / (pow(x + par[3], 2) + par[4])
    + FB;
  return return_value;
}

double TReddening::Prevot(double lambda) {
  return Prevot_law(lambda);
}

double TReddening::Prevot_law(double lambda) {
  double ext = 0;
  //if (lambda>=6.10e-8 && lambda<1.46e-7) 
  if (lambda>=0 && lambda<1.46e-7) 	
    { 
      ext = 39.8158 + lambda*(-1.84437e8);
    }
  if (lambda>=1.46e-7 && lambda<3.1e-7)
    { 
      ext  = 1.92234e14*lambda*lambda + (-1.32706e8)*lambda + 27.903;
    }
  
  if (lambda>=3.10e-7 && lambda<5.30e-7) 
    { 
      ext = 8.542 + lambda*(-1.1122e7);
    }
  
  if (lambda>=5.30e-7 && lambda<1.25e-6) 
    { 
      ext = 4.20257 + lambda*(-2.80205e6);
    }
  return ext*(4.05/2.72); //rescale on calz Av
}

double TReddening::Prevot_file(double lambda) {
  lambda *= 1e10;
 
 if (tab_ext_prev.size()==0)
    ReadTabPrev();
  
  int j = 0;
  int n = tab_ext_prev.size();
  while(tab_lambda_prev[j]<lambda && j<n) {
    j++;
  }
  
  if(j>=n)
    return tab_ext_prev[n-1];
  
  if(j==0){
    return tab_ext_prev[0];}
  
  double weight=(lambda - tab_lambda_prev[j-1])/(tab_lambda_prev[j]-tab_lambda_prev[j-1]);
  double ext=tab_ext_prev[j]*weight + tab_ext_prev[j-1]*(1-weight);
  
  return ext;
}

double TReddening::Late(double lambda) {
  //return Cardelli(lambda);
  return Prevot(lambda);
}

double TReddening::Cardelli(double lambda) {// cardelli law
  double return_value = 0.;
  lambda *= 1.e6; //on passe en microns
  
  
  double cardelli_highl_param[] = { 0.574, 1.61, -0.527, fRV_card };
  double a_law_param[] = { 1, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979,
			   -0.77530, 0.32999 };
  double b_law_param[] = { 1.41338, 2.28305, 1.07233, -5.38434, -0.62251,
			   5.30260, -2.09002 };
  double a_lowVL_param[] = { 1.752, -0.316, 0.104, -4.67, 0.341, -0.04473,
			     -5.9, -0.009779 };
  double b_lowVL_param[] = { -3.090, 1.825, 1.206, -4.62, 0.263, 0.2130,
			     -5.9, 0.1207 };
  double a_FARUV_param[] = { -1.073, -0.628, 0.137, -0.070, -8. };
  double b_FarUV_param[] = { 13.670, 4.257, -0.420, 0.374, -8. };
  
  if (lambda >= 0.1 && lambda < 0.125)
    return_value = a_FarUV(&lambda, a_FARUV_param) + b_FarUV(&lambda,
							     b_FarUV_param) / fRV_card;
  if (lambda >= 0.125 && lambda < 0.3)
    return_value = a_lowVL(&lambda, a_lowVL_param) + b_lowVL(&lambda,
							     b_lowVL_param) / fRV_card;
  
  if (lambda >= 0.3 && lambda < 0.9)
    return_value = a_lowl(&lambda, a_law_param) + b_lowl(&lambda,
							 b_law_param) / fRV_card;
  
  if (lambda >= 0.9 && lambda <= 3.)
    return_value = cardelli_highl(&lambda, cardelli_highl_param);
  //	cout<<" Rv = "<<fRV_card<<" value = "<<return_value<<endl;
  return_value *= fRV_card;
  return return_value;
}

double TReddening::Calzetti(double lambda) {
  return Calzetti_old(lambda);
}


double TReddening::Calzetti_newSB1(double lambda) { 
  lambda *= 1e10;

  if (tab_ext_calz1.size()==0){
    ReadTabCalz1();}

  int j = (lambda-600)/10;
  int n = tab_ext_calz1.size();
  
  if(j>=n-1){
    return tab_ext_calz1[n-1];}
  
  if(j<1){
    return tab_ext_calz1[0];}
  
  double weight=(lambda - tab_lambda_calz1[j])/10;
  double ext=tab_ext_calz1[j]*(1-weight) + tab_ext_calz1[j+1]*weight;

  return ext;
}

double TReddening::Calzetti_newSB2(double lambda) { 
  lambda *= 1e10;

  if (tab_ext_calz2.size()==0){
    ReadTabCalz2();}

  int j = (lambda-600)/10;
  int n = tab_ext_calz2.size();
  
  if(j>=n-1){
    return tab_ext_calz2[n-1];}
  
  if(j<1){
    return tab_ext_calz2[0];}
  
  double weight=(lambda - tab_lambda_calz2[j])/10;
  double ext=tab_ext_calz2[j]*(1-weight) + tab_ext_calz2[j+1]*weight;

  return ext;
}

double TReddening::Calzetti_newSB3(double lambda) { 
  lambda *= 1e10;

  if (tab_ext_calz3.size()==0){
    ReadTabCalz3();}

  int j = (lambda-600)/10;
  int n = tab_ext_calz3.size();
  
  if(j>=n-1){
    return tab_ext_calz3[n-1];}
  
  if(j<1){
    return tab_ext_calz3[0];}
  
  double weight=(lambda - tab_lambda_calz3[j])/10;
  double ext=tab_ext_calz3[j]*(1-weight) + tab_ext_calz3[j+1]*weight;

  return ext;
}

double TReddening::Calzetti_old(double lambda) { // calzetti law
  lambda *= 1.e6;
  double return_value = 0.;
  double derivee_120 = 2.659 * (-1.509 / pow(0.120, 2) + 2 * 0.198 / pow(0.120, 3) - 3. * 0.011 / pow(0.120, 4));
  double k_120 = 2.659 * (-2.156 + 1.509 / 0.120 - 0.198 / pow(0.120, 2)
			  + 0.011 / pow(0.120, 3)) + fRV_calz;
  double derivee_2200 = -2.659 * 1.040 / 2.200 / 2.200;
  double k_2200 = 2.659 * (-1.857 + 1.040 / 2.200) + fRV_calz;
  
  if (lambda >= 0.120 && lambda < 0.630)
    return_value = 2.659 * (-2.156 + 1.509 / lambda - 0.198
			    / pow(lambda, 2) + 0.011 / pow(lambda, 3)) + fRV_calz;
  if (lambda >= 0.630 && lambda <= 2.200)
    return_value = 2.659 * (-1.857 + 1.040 / lambda) + fRV_calz;
  if (lambda < 0.120)
    return_value = derivee_120 * (lambda - 0.120) + k_120;
  if (lambda > 2.200)
    return_value = derivee_2200 * (lambda - 2.200) + k_2200;
  
  return return_value;
}

void TReddening::SetRVCardelli(double RV) {
  fRV_card = RV;
}

void TReddening::SetRVCalzetti(double RV) {
  fRV_calz = RV;
}

double TReddening::AttenuationFactor(double ebv, double lambda, string sedext) {
  
  double law = 0;

  if (ebv<=0)
    return 1;
  
  if (sedext == "Null") {
    return 1;
  }
  
  if (sedext == "Cardelli" || sedext == "Card"|| sedext == "card"|| sedext == "cardelli") {
    law = Cardelli(lambda);
  }
  else if (sedext == "Prevot" || sedext == "Prev" || sedext == "prevot" || sedext == "prev") {
    law = Prevot_law(lambda);
  }
  else if (sedext == "CalzFile" || sedext == "Calzetti_File" || sedext == "Calzetti_file" || sedext == "calzetti_file" ||  sedext == "calz_file" || sedext == "Calz_file" || sedext == "Calz_File" || sedext == "Calz") {
    law = Calzetti_newSB1(lambda);  
  }
  else if (sedext == "CalzLaw" || sedext == "Calzetti_Law" || sedext == "Calzetti_law" || sedext == "calzetti_law" ||  sedext == "calz_law" || sedext == "Calz_law" || sedext == "Calz_Law" ) {
    law = Calzetti_old(lambda);  
  }
  else if (sedext == "Calz_bump1" || sedext == "calz_bump1" || sedext == "Calzetti_Bump1" || sedext == "calzetti_bump1" || sedext == "Calzetti_bump1" || sedext == "Calz_Bump1" || sedext == "CalzBump1" ) {
    law = Calzetti_newSB2(lambda);
  }
  else if (sedext == "Calz_bump2" || sedext == "calz_bump2" || sedext == "Calzetti_Bump2" || sedext == "calzetti_bump2" || sedext == "Calzetti_bump2" || sedext == "Calz_Bump2" || sedext == "CalzBump2" ) {
    law = Calzetti_newSB3(lambda);
  }
  else
    cout << "law " << sedext << " not known" << endl;
  
  //cout << " ext = " << sedext << endl;

  tempeval = pow(10., -0.4 * ebv * law);
  if (TMath::IsNaN(tempeval))
    cout << " ebv = " << ebv << " sedext = " << sedext << " lambda = "
	 << lambda << " law = " << law << endl;
  return pow(10., -0.4 * ebv * law);
}
