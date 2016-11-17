/*
 * TFilter.C
 *
 *  Created on: Mar 26, 2010
 *      Author: gorecki
 *  Modified on: Sep 19 2014 by JS Ricol
 */

#include "TFilter.h"
#include "TList.h"
#include "TKey.h"

using namespace std;

TFilter::TFilter() :
  epsilon(1.e-10) {
  multifilter = 0;
  rand = new TRandom3();
  rand -> SetSeed(0);
  dT = 0;
  dlambda = 0;
}

TFilter::~TFilter() {
  delete lambdamin;
  delete lambdamax;
  delete norm;
  ffilename = NULL;
  if (fextension == "tf1"){
    for (int ifile = 0; ifile < nfilter; ifile++) {
      f_filter[ifile]->Close();
    }
  }
}

void TFilter::LoadDataCard(TDataCard *card){
  //cout << card->filterlambdaunit << " " << card->filterformat << " " << card->nfilter << endl;
  SetConversionUnitLambda(card->filterlambdaunit);
  SetFileType(card->filterformat);
  LoadFile(card->nfilter, card->filterlist);
}

void TFilter::Draw(string opt){
  int col[6]= {619, 601, 419, 798, 801, 632};
  int fcol;
  for (int ifilter = 0; ifilter < nfilter; ifilter++){
    if (ifilter<6)
      fcol = col[ifilter];
    else
      fcol=1;
    ffilter_fromfile_tf1[ifilter] -> SetLineColor(fcol);
    ffilter_fromfile_tf1[ifilter] -> SetLineStyle(1);
    ffilter_fromfile_tf1[ifilter] -> SetLineWidth(3);
    ffilter_fromfile_tf1[ifilter] -> SetMaximum(1.);
    ffilter_fromfile_tf1[ifilter]->Draw(opt.c_str());
    opt="same";
  }
}


void TFilter::SetFileType(string extension) {
  fextension = extension;
  
  //cout << "extension : " << fextension << " " << multifilter << endl;
  
  if (fextension == "tf1list")
    multifilter = 1;
}

string TFilter::GetFileType(){
  return fextension;
}

void TFilter::LoadFile(int nfile, string *filename) { // Load files

  //cout << "tfilter : " << nfile << " " << fextension << endl;
  
  nfilter = nfile;

  lambdamin = new double[nfilter];
  lambdamax = new double[nfilter];
  norm = new double[nfilter];
  
  //cout << "Loading filters ... " << endl;
  ffilename = new string[nfile];
  
  for (int i = 0; i < nfile; i++) {
    ffilename[i] = filename[i];
    //cout << i << " " << ffilename[i] << endl;
  }
  
  cout << nfile << " " << fextension << endl;
  
  if (fextension == "txt") {// extension file : txt
    vector<double> lambda[nfilter];
    vector<double> filter[nfilter];
    ffilter_fromfile = new ROOT::Math::Interpolator*[nfilter];
    for (int ifile = 0; ifile < nfilter; ifile++) {
      //cout << filename[ifile] << endl;
      ifstream file(filename[ifile].c_str());
      if (!file.is_open()) {
	cout << "File " << file << " " << filename[ifile] << " could not be open"
	     << endl;
      } else {
	double entry[2];
	while (true) {
	  file >> entry[0];
	  if (file.eof())
	    break;
	  file >> entry[1];

	  lambda[ifile].push_back(entry[0]);
	  filter[ifile].push_back(entry[1]);
	  lambdamax[ifile] = entry[0]*fconversionUnit_Lambda;
	  
	}
      }
      
      lambdamin[ifile] = lambda[ifile][0]*fconversionUnit_Lambda;
      
      /*
      /cout << ifile << " " << filename[ifile]
	    << "; lambda range = [" << lambdamin[ifile] << ", " << lambdamax[ifile] << "]"
	    << " npoints = " << filter[ifile].size() << endl;
      */
      
      ffilter_fromfile[ifile] = new ROOT::Math::Interpolator(lambda[ifile], filter[ifile],
							     ROOT::Math::Interpolation::kLINEAR);
      
    }
    f_filter = NULL;
    ffilter_fromfile_tf1 = NULL;
    ffilter_fromfile_ttree = NULL;
  }
  if (fextension == "tf1" || fextension == "tf1list" ) {// extension file : root contains TF1
    //		cout << " fextension = " << fextension << endl;
    f_filter = new TFile*[nfilter];
    ffilter_fromfile_tf1 = new TF1*[nfilter];
    list = new TList*[nfilter];
    
    for (int ifile = 0; ifile < nfilter; ifile++) {
      f_filter[ifile] = new TFile(filename[ifile].c_str());
      list[ifile] = f_filter[ifile] -> GetListOfKeys();
      string name = "aaa";
      for (int j = 0; j < list[ifile]->GetEntries(); ++j) {
	key = (TKey*) list[ifile]->At(j);
	if (TClass::GetClass(key->GetClassName())->InheritsFrom("TF1")) {
	  name = key->GetName();
	}
	if (name == "aaa") {
	  cout << " Wrong extension  " << filename
	       << " file does not contain TF1 object" << endl;
	  break;
	} else {
	  // cout << " name = " << name << endl;
	  ffilter_fromfile_tf1[ifile] = (TF1*) f_filter[ifile] -> Get(name.c_str());
	  //cout << "LOADING " << name << endl;
	  break;
	}
      }
      
      double dl=1e-7;
      double l=0;
      while (dl>1e-10){
	l+=dl;
	if (ffilter_fromfile_tf1[ifile]->Eval(l/fconversionUnit_Lambda)>0){
	  l-=dl;
	  dl/=10;
	}
      }
      lambdamin[ifile]=l;
      dl=1e-7;
      l=3e-5;
      while (dl>1e-10){
	l-=dl;
	// cout << l << " " << ffilter_fromfile_tf1[ifile]->Eval(l/fconversionUnit_Lambda) << endl;
	if (ffilter_fromfile_tf1[ifile]->Eval(l/fconversionUnit_Lambda)>0){
	  l+=dl;
	  dl/=10;
	}
      }
      lambdamax[ifile]=l;
      //cout << ifile << " lambdarange  = [" << lambdamin[ifile] << ", " << lambdamax[ifile] << "]" << endl;
    }
    ffilter_fromfile = NULL;
    ffilter_fromfile_ttree = NULL;
  }
  //	cout << " lambdamin = " << lambdamin << " lambdamax = " << lambdamax << endl;
}

void TFilter::NewFilters(){
  if (fextension == "tf1list") {
    for (int ifile = 0; ifile < nfilter; ifile++) {
      int n = list[ifile]->GetEntries();
      int j = rand->Integer(n);
      
      string name = list[ifile]->At(j)->GetName();
      delete ffilter_fromfile_tf1[ifile];
      ffilter_fromfile_tf1[ifile] = (TF1*) f_filter[ifile] -> Get(name.c_str());
      // cout << "new filter : " << name << endl;
    }
  }
}

double TFilter::GetLambdaMin(int ifilter){
  return lambdamin[ifilter];
}
double TFilter::GetLambdaMax(int ifilter){
  return lambdamax[ifilter];
}

void TFilter::SetDLamdaNanometer(double dl){
  dlambda = dl/1e9;
}

void TFilter::SetDTransmission(double dt){
  dT = dt;
}

double TFilter::Eval(double lambda, int ifilter) {
  // gives flux vs lambda
  
  //cout << " lambda = " << lambda << "; convUnit = " << fconversionUnit_Lambda << endl;

  if (lambda<lambdamin[ifilter] || lambda>lambdamax[ifilter])
    return 0;
  
  lambda /= fconversionUnit_Lambda;

  double return_value = 0;
  if (fextension == "txt") {
    return_value = ffilter_fromfile[ifilter] -> Eval(lambda-dlambda);
  }
  if (fextension == "tf1" || fextension == "tf1list" || fextension == "spline") {
    return_value = ffilter_fromfile_tf1[ifilter] -> Eval(lambda-dlambda);
  }
  
  return_value*=(1+dT);
  return return_value;
}

void TFilter::SetConversionUnitLambda(double conversionLambda_filter) {
  fconversionUnit_Lambda = conversionLambda_filter;
}

void TFilter::Lambda_Effective(double *lambda_eff) {
  
  for (int filter = 0; filter < nfilter; filter++) {
    lambda_eff[filter] = 0.;
    double normalization = 0;
    for (double lambda = lambdamin[filter]; lambda <= lambdamax[filter]; lambda += epsilon) {
      lambda_eff[filter] += Eval(lambda, filter) * epsilon * lambda;
      normalization += Eval(lambda, filter) * epsilon;
    }
    lambda_eff[filter] /= normalization;/*en m*/
    norm[filter] = normalization;
  }
}

double TFilter::GetNorm(int i){
  return norm[i];
}

double TFilter::Integral_Eval_FSurlambdaDlambda(int filter) {
  double return_value = 0;
  //cout << "lambda_min = " << lambdamin << " lambda_max = " << lambdamax << "; epsilon = " << epsilon << endl;
  for (double lambda = lambdamin[filter]; lambda <= lambdamax[filter]; lambda += epsilon) {
    return_value += Eval(lambda, filter) * epsilon / lambda;
  }
  return return_value;
}

bool TFilter::MultiFilter(){
  return multifilter;
}
