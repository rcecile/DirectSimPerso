/*
 * TSed.C
 *
 *  Created on: Mar 26, 2010
 *      Author: gorecki
 *  Modified on: Sep 19, 2014 by JS Ricol
 */
#include "TSed.h"
#include "TList.h"
#include "TF1.h"
#include "TKey.h"
#include "TDataCard.h"
#include <string>

using namespace std;

////////////////////////////////////////////////////////////////
//                                                            //
// This is the description block.                             //
//                                                            //
////////////////////////////////////////////////////////////////

TSed::TSed(bool p){
  fmix=0;
  ebvmax_early = 0.1;
  ebvmax_late = 0.3;
  ebvmax_sb = 0.3;
  ebvmax_default = 0.3;
  rand = new TRandom3();
  print = p;
}

void TSed::LoadDataCard(TDataCard *card){
  SetTypes(card);
  SetFileType(card->sedformat);
  LoadFile(card);
  LoadType(card);
  LoadIGM(card);
}

void  TSed::LoadIGM(TDataCard *card){
  SetIGM(1);
  figmTable = new TIGM();
  figmTable->SetIGMMadau(card->igmtable);
  figmTable->LoadIGMTable();
}

void TSed::SetTypes(TDataCard *card){
  SetMix(card->fmix);
  nsed = card->nsed;
  
  // cout << "ntype : " << nsed << endl;

  ntype = nsed;
  
  if (fmix)
    ntype=(nsed-1)*10+1;
}

int TSed::GetNtype(){
  return ntype;
}

string *TSed::GetTypeArray(){
  return stype;;
}

void TSed::SetFileType(string extension) {
  fextension = extension;
  if (fextension=="dat" || fextension=="sed" || fextension=="ascii")
    fextension="txt";
}

void TSed::LoadType(TDataCard *card){
  string *sedtype = card->sedtype;
  string *sedext = card->sedext;

  if (print){
    cout << "ntype = " << ntype << endl;
    cout << sedtype << endl;
  }
  
  gtype = new int[ntype];
  stype = new string[ntype];
  sext = new string[ntype];
  ebvmax = new double[ntype];

  if (card->ebvmax_early_defined)
    ebvmax_early = card->ebvmax_early;
  if (card->ebvmax_late_defined)
    ebvmax_late = card->ebvmax_late;
  if (card->ebvmax_sb_defined)
    ebvmax_sb = card->ebvmax_sb;
  ebvmax_default = card->ebvmax_default;
  
  if (print)
    cout << "ntype = " << ntype << endl;

  for (int i=0; i<ntype; i++){
    if (!fmix){
      stype[i] = sedtype[i];
      sext[i] = sedext[i];
    }
    else{
      int k = i/10;
      int r = i%10;
      if (r>=5)
	k++;
      stype[i] = sedtype[k];
      sext[i] = sedext[k];
    }

    //cout << "Sed " << i << " type = " << stype[i] << endl;

    if (stype[i] == "Early")
      SetMaxEBVforType(i, ebvmax_early);
    else if (stype[i] == "Late")
      SetMaxEBVforType(i, ebvmax_late);
    else if (stype[i] == "SB")
      SetMaxEBVforType(i, ebvmax_sb);
    else
      SetMaxEBVforType(i, ebvmax_default);
  }
  BuildTypeArray();
}

void TSed::BuildTypeArray(){
  list_early.reserve(ntype);
  list_late.reserve(ntype);
  list_sb.reserve(ntype);
  for (int i=0; i<ntype; i++){
    if (stype[i] == "Early"){
      list_early.push_back(i);
    }
    else if (stype[i] == "Late"){
      list_late.push_back(i);
    }
    else if (stype[i] == "SB"){
      list_sb.push_back(i);
    }
  }
}

int TSed::TransformTypeFromCWW(int t){
  vec = &list_early;
  if (t>=5)
    vec = &list_late;
  if (t>=25)
    vec = &list_sb;
  
  float r = rand->Uniform(0, vec->size());
  if (r<0)
    r=0;
  if (r>=vec->size())
    r-=0.5;
  int ivec = int(r);
  
  t = (*vec)[ivec];
  if (t<0 || t>=ntype){
    cout << "WRONG TYPE" << endl;
    cout << t << " " << vec->size() << " " << r << " " << ivec << " " << (*vec)[ivec] << endl;
    cout << "type set to average" << endl;
    t = ntype/2;
  }

  return t;
}

double TSed::GetEBV(int type){
  return FlatEBV(type);
}

double TSed::FlatEBV(int type){
  return rand->Uniform(0, ebvmax[type]);
}

void TSed::SetMaxEBVforType(int type, double ebv){
  ebvmax[type] = ebv;
}


void TSed::LoadFile(TDataCard *card) { // Load files
  string *filename = card->sedlist;
  lambdaunit = card->sedlambdaunit;
  lambdamin = 100;
  lambdamax = 0;
  
  if (fextension == "txt") {// extension file : txt
    //		cout << " fextension = " << fextension << endl;
    vector<double> lambda[nsed];
    vector<double> sed[nsed];
    
    lambdamin_ = new double[nsed];
    lambdamax_ = new double[nsed];
    
    fsed_fromfile = new ROOT::Math::Interpolator*[nsed];

    int lambdavec_size = 3000;

    for (int ifile = 0; ifile < nsed; ifile++) {
      ifstream file(filename[ifile].c_str());
      
      lambdamin_[ifile] = 100;
      lambdamax_[ifile] = 0;

      lambda[ifile].reserve(lambdavec_size);
      sed[ifile].reserve(lambdavec_size);

      if (!file.is_open()) {
	cout << "File " << filename[ifile] << " could not be open"
	     << endl;
      } else {
	double entry[2];
	string line;
	int nl=0;
	while ( getline (file,line) && nl<lambdavec_size){
	  //cout << line << endl;
	  int found = line.find("#");
	  if (found>=0)
	    continue;
	  std::istringstream iss(line);
	  for (int i=0; i<2; i++)
	    iss >> entry[i];
	  if (file.eof())
	    break;
	  entry[0]*=lambdaunit;
	  
	  //if (ifile==0 && print)
	  //cout << "Reading " << nl << " : " << entry[0] << " " << entry[1] << endl;
	  
	  if ((lambda[ifile])[nl-1]<entry[0]){
	    lambda[ifile].push_back(entry[0]);
	    sed[ifile].push_back(entry[1]);
	    
	    nl++;
	    
	    if (lambdamin_[ifile] > entry[0])
	      lambdamin_[ifile] = entry[0];
	    if (lambdamax_[ifile] < entry[0])
	      lambdamax_[ifile] = entry[0];
	  }
	}
      }
      
      if (ifile==0 && print)
	cout << "done : " << (lambda[ifile])[0] << " " << (lambda[ifile])[(lambda[ifile]).size()-1] << endl;
      
      if (print)
	cout << filename[ifile] << " lambda range = [" << lambdamin_[ifile] << ", " <<  lambdamax_[ifile] << "] "
	     << (lambda[ifile]).size() << " points" << endl;
      
      fsed_fromfile[ifile] = new ROOT::Math::Interpolator(lambda[ifile],
							  sed[ifile],
							  ROOT::Math::Interpolation::kLINEAR);
      
      
      if (print)
	cout << "Eval @4000A� --> " <<  fsed_fromfile[ifile]->Eval(4e-7) << endl;
      
      
    }
    f_sed = NULL;
    fsed_fromfile_tf1 = NULL;
  }
  if (fextension == "tf1") {// extension file : root contains TF1
    //cout << " fextension = " << fextension << endl;
    f_sed = new TFile*[nsed];
    fsed_fromfile_tf1 = new TF1*[nsed];
    for (int ifile = 0; ifile < nsed; ifile++) {
      f_sed[ifile] = new TFile(filename[ifile].c_str());
      TList* list = f_sed[ifile] -> GetListOfKeys();
      string name = "aaa";
      for (int j = 0; j < list->GetEntries(); ++j) {
	TKey* key = (TKey*) list->At(j);
	if (TClass::GetClass(key->GetClassName())->InheritsFrom("TF1")) {
	  name = key->GetName();
	}
	if (name == "aaa") {
	  cout << " Wrong extension  " << filename
	       << " file does not contain TF1 object" << endl;
	  break;
	} else {
	  fsed_fromfile_tf1[ifile] = (TF1*) f_sed[ifile] -> Get(name.c_str());
	  if (print)
	    cout << "sed eval @ 4000 A : " << fsed_fromfile_tf1[ifile]->Eval(4e-7) << endl;
	}
      }
    }
    fsed_fromfile = NULL;
  }
  if (fextension == "spline") {// extension file : root contains TF1
    fsed_fromfile = NULL;
  }
}

double TSed::GetLambdaMin(int ised){
  return lambdamin_[ised];
}

double TSed::GetLambdaMax(int ised){
  return lambdamax_[ised];
}

void TSed::SetReddening(TReddening *red) {
  freddening = red;
}

void TSed::SetIGM(bool igm) {
  figm = igm;
}
void TSed::LoadIGMTable(TIGM *table) {
  if (figm) {
    figmTable = table;
    
  } else
    figmTable = NULL;
}
void TSed::SetEBV(double ebv) {
  febv = ebv;
}

void TSed::SetEmittedRedshift(double z) {
  figmredshift = z;
}

void TSed::SetMix(bool mix){
  fmix=mix;
}


string TSed::SedType(int type){
  return stype[type];
}

string TSed::SedExt(int type){
  return sext[type];
}

double TSed::Eval(double lambda, int itype) {
  int ised = itype;
  int sedlow = itype;
  float mix=1;
  
  if (fmix){
    // cout << "fmix " << itype%10 << endl;
    if (itype % 10 == 0) {
      ised = int(itype / 10.);
    }
    else {
      sedlow = int(itype / 10);
      mix = 1 - (itype) / 10. + sedlow;
    }
  }
  
  double return_value = 0;
  
  if (fextension == "txt") {
    if (lambda >= lambdamin_[ised] && lambda <= lambdamax_[ised]) {
      if (mix == 1){
	return_value = fsed_fromfile[ised] -> Eval(lambda)* freddening -> AttenuationFactor(febv, lambda, SedExt(itype));
      }
      else
	return_value = fsed_fromfile[sedlow] -> Eval(lambda) * mix * freddening -> AttenuationFactor(febv, lambda, SedExt(sedlow*10))
	  + (1. - mix) * fsed_fromfile[sedlow + 1] -> Eval(lambda) * freddening -> AttenuationFactor(febv, lambda, SedExt((sedlow+1)*10));
    }
  }
  
  if (fextension == "tf1" || fextension == "spline") {
    if (mix == 1){
      return_value = fsed_fromfile_tf1[ised] -> Eval(lambda) * freddening -> AttenuationFactor(febv, lambda, SedExt(itype));
    }
    else{
      return_value = fsed_fromfile_tf1[sedlow] -> Eval(lambda) * mix *  freddening -> AttenuationFactor(febv, lambda, SedExt(sedlow*10))
	+ (1. - mix) * fsed_fromfile_tf1[sedlow + 1] -> Eval(lambda) * freddening -> AttenuationFactor(febv, lambda, SedExt((sedlow+1)*10));
    }
  }
  
  igm = 1;
  if (figm) {
    igm = figmTable -> IGMTableEval(lambda, figmredshift);
  }
  
  return_value *= igm;
  
  return return_value;
}

