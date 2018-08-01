/*
 * TDataCard.C
 *
 *  Created on: Sept. 18, 2014
 *      Author: JS Ricol
 */

#include "TDataCard.h"
#include "fstream"

#include "iostream"
#include <math.h>
#include <sstream>
#include <string.h>
#include <unistd.h>

#include <dirent.h>

using namespace std;

TDataCard::TDataCard(bool p) {
  zgrid[0] = 0;
  zgrid[1] = 4.5;
  zgrid[2] = 0.01;
  ebvgrid[0] = 0;
  ebvgrid[1] = 0;
  ebvgrid[2] = 1;
  deltatypegrid = 1;
  ebvmax_early_defined = 0;
  ebvmax_late_defined = 0;
  ebvmax_sb_defined = 0;
  fitonlytype = 0;
  ebvmax_default = 0.3;
  
  CFHTLS_ErrMagFile = "/sps/lsst/PhotozBAO/CFHTLS/ErrMag_CFHTLS.root";

  inputdir_am = "franzona";
  print = p;
}

TDataCard::~TDataCard() {
  delete filterlist;
  delete sedlist;
}

void TDataCard::Print(){
  cout << "##################################################################" << endl;
  cout << "###      Input Parameters in " << filename << endl;
  cout << "###      Survey = " << survey << endl;
  cout << "###      PhotoZ/AppMag dir = " << outputdir_pz << " " << outputdir_am << endl;
  cout << "###      Filters : " << nfilter << " " << filterformat << " " << endl;
  for (int i=0; i<nfilter; i++)
    cout << "###               " << filterlist[i] << endl;
  cout << "###      BFilter : " << bfiltername << endl;
  cout << "###      Z = [" << zmin << ", " << zmax << "] zbin=" << zbin << endl;
  cout << "###      EBV = [" << ebvmin << ", " << ebvmax << "] ebvbin=" << ebvbin << endl;
  cout << "###      SEDs : " << nsed << " " << sedformat << endl;
  for (int i=0; i<nsed; i++)
    cout << "###               " << sedlist[i] << endl;
  cout << "###      MIX = " << fmix << endl;
  cout << "###      TablesDir = " << tabledir << endl;
  cout << "###      FluxTable = " << fluxfile << endl;
  cout << "###      KcorrTable = " << kcorrfile << endl;
  cout << "###      Prior = " << prior << " file = " << priorfile << endl;
  cout << "#################################################################" << endl;
}

double *TDataCard::GetParRange(){
  double *vec = new double[6];
  int i=0;
  vec[i++] = zmin;
  vec[i++] = zmax;
  vec[i++] = 0;
  vec[i++] = ntype-1;
  vec[i++] = ebvmin;
  vec[i++] = ebvmax;
  return vec;
}

string TDataCard::Getbfiltername(){
  return bfiltername;
}

int TDataCard::Read(string name) {
  if (access(name.c_str(), 0 ) != 0){
    cout << "can not read cardfile " << name << endl;
    return -1;
  }
  filename = name;
  ifstream file(name.c_str());
  string par;
  
  smix = "0";
  fmix=0;
  savePDF = 0;

  while (!file.eof()){
    file >> par;
    
    //cout << "par : " << par << endl;

    if (file.eof())
      break;
    if (par=="survey:")
      file >> survey;
    if (par=="outputdir_pz:")
      file >> outputdir_pz;
    if (par=="outputdir_am:")
      file >> outputdir_am;
    if (par=="inputdir_am:")
      file >> inputdir_am;
    if (par=="calib0point:")
      file >> calib0point;
    if (par=="filterformat:")
      file >> filterformat;
    if (par=="filterlambdaunit:")
      file >> filterlambdaunit;
    if (par=="nfilter:")
      file >> nfilter;
    if (par=="filterlist:"){
      filterlist = new string[nfilter];
      for (int i=0; i<nfilter; i++){
	file >> filterlist[i];
	if (print)
	  cout << i << " " << filterlist[i] << endl;
      }
    }
    if (par=="filterdir:"){
      string filterdir;
      file >> filterdir;
      filterlist = new string[nfilter];
      string postfix = ".txt";
      if (filterformat=="tf1")
	postfix=".root";
      DIR *dpdf;
      struct dirent *epdf;
      dpdf = opendir(filterdir.c_str());
      if (dpdf != NULL){
	int i=0;
	while (epdf = readdir(dpdf)){
	  string name = epdf->d_name;
	  if (name.find(postfix)!=std::string::npos){
	    if (print)
	      cout << name << endl;
	    filterlist[i]=filterdir+"/"+name;
	    i++;
	    //cout << i << " " << nfilter << endl;
	    if (i==nfilter)
	      break;
	  }
	}
      }
    }
    if (par=="bfilter:")
      file >> bfiltername;
    if (par=="bfilterlambdaunit:")
      file >> bfilterlambdaunit;
    if (par=="Zrange:")
      file >> zmin >> zmax >> zbin;
    if (par=="EBVrange:")
      file >> ebvmin >> ebvmax >> ebvbin;
    if (par=="nsed:")
      file >> nsed;
    if (par=="sedlambdaunit:")
      file >> sedlambdaunit;
    if (par=="sedformat:")
      file >> sedformat;
    if (par=="sedlist:"){
      sedlist = new string[nsed];
      sedtype = new string[nsed];
      sedext = new string[nsed];
      string line;
      getline(file, line);
      for (int i=0; i<nsed; i++){
	if (!getline(file, line))
	  break;
	std::istringstream iss(line);
	iss >> sedlist[i] >> sedtype[i] >> sedext[i];
	if (sedext[i]==""){
	  if (sedtype[i]=="Early" || sedtype[i] =="Late")
	    sedext[i] = "Cardelli";
	  if (sedtype[i]=="SB")
	    sedext[i]="CalzLaw";
	}
	cout << "Ext : " << i << " " << sedtype[i] << " " << sedext[i] << endl;
      }
    }
    if (par=="sedfile:"){
      string sedfilename;
      file >> sedfilename;
      ifstream sedfile(sedfilename.c_str());
      sedlist = new string[nsed];
      sedtype = new string[nsed];
      sedext = new string[nsed];
      string line;
      for (int i=0; i<nsed; i++){
	if (!getline(sedfile, line))
	  break;
	std::istringstream iss(line);
	iss >> sedlist[i] >> sedtype[i] >> sedext[i];
	if (sedext[i]==""){
	  if (sedtype[i]=="Early" || sedtype[i] =="Late")
	    sedext[i] = "Cardelli";
	  if (sedtype[i]=="SB")
	    sedext[i]="CalzLaw";
	}
      }
      sedfile.close();
    } 
    if (par=="ebvmax_early:"){
      file >> ebvmax_early;
      ebvmax_early_defined = 1;
    }
    if (par=="ebvmax_late:"){
      file >> ebvmax_late;
      ebvmax_late_defined = 1;
    }
    if (par=="ebvmax_sb:"){
      file >> ebvmax_sb;
      ebvmax_sb_defined = 1;
    }
    if (par=="ebvmax_default:"){
      file >> ebvmax_default;
    }
    if (par=="ebvtreatment:")
      file >> ebvtreatment;
    if (par=="sortedtypes:")
      file >> sortedtypes;
    if (par=="sedmix:")
      file >> smix;
    if (par=="IGMtable:")
      file >> igmtable;
    if (par=="tabledir:")
      file >> tabledir;
    if (par=="FluxTable:")
      file >> fluxfile;
    if (par=="KcorrTable:")
      file >> kcorrfile;
    if (par=="nVisit:"){
      nVisit = new double[nfilter];
      for (int i=0; i<nfilter; i++){
	file >> nVisit[i];
      }
    }
    if (par=="CFHTLS_ErrMagFile:")
      file >> CFHTLS_ErrMagFile;
    if (par=="libcat:")
      file >> libcat;
    if (par=="prior:")
      file >> prior;
    if (par=="priorfile:")
      file >> priorfile;
    
    if (par=="zgrid:")
      file >> zgrid[0] >> zgrid[1] >> zgrid[2];
    if (par=="ebvgrid:")
      file >> ebvgrid[0] >> ebvgrid[1] >> ebvgrid[2];
    if (par=="deltatypegrid:")
      file >> deltatypegrid;
    
    if (par=="fitonlytype:")
      file >> fitonlytype;

    if (par=="selcut:")
      file >> selcut;

    if (par=="BDTdir:")
      file >> BDTdir;
    if (par=="savePDF:")
      file >> savePDF;
 
  }
   
  ntype = nsed;
  if (smix=="mix"){
    fmix=1;
    ntype = (nsed-1)*10+1;
    cout << "nnnntype = " << ntype << " " << smix << endl;
  }

  if (access(outputdir_am.c_str(), 0 ) != 0 ){
    cout << "WARNING output dir specified in DataCard  : " << outputdir_am << ", does not exist" << endl;
    return -1;
  }
  
  if (selcut.find("gold") != std::string::npos || selcut.find("GOLD") != std::string::npos || selcut.find("Gold") != std::string::npos)
    selcut = "golden";
  
  if (print){
    cout << "prior = " << prior << endl;
    cout << "priorfile = " << priorfile << endl;
  }


  cout << "termine" << endl;

  file.close();
  return 1;
  
}

