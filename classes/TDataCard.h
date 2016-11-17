#ifndef TDataCard_H_
#define TDataCard_H_

//___________________________________
// TDataCard.h
//
//
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include <dirent.h>

using namespace std;

class TDataCard {
 public:
  //constructor
  TDataCard(bool print=1);
  //destructor
  ~TDataCard();
  // read
  int Read(string file);
  double* GetParRange();
  void Print();
  
  // DataCard filename
  string filename;

  // survey
  string survey;
  
  // output
  string inputdir_am;
  string outputdir_am;
  string outputdir_pz;
  
  // calib
  string calib0point;

  // Z range
  double zmin;
  double zmax;
  double zbin;
  
  // ebv range
  double ebvmin;
  double ebvmax;
  double ebvbin;
    
  // filters
  string filterformat;
  double filterlambdaunit;
  double bfilterlambdaunit;
  int nfilter;
  string *filterlist;
  

  // SEDs
  int nsed;
  int ntype;
  string sedformat;
  double sedlambdaunit;
  string *sedlist;
  string *sedtype;
  string *sedext;
  string smix;
  bool fmix;
  string ebvtreatment;
  string sortedtypes;
  bool ebvmax_early_defined, ebvmax_late_defined, ebvmax_sb_defined;
  double ebvmax_early, ebvmax_late, ebvmax_sb, ebvmax_default;
  
  // IGM
  string igmtable;

  // Tables dir
  string tabledir;
  
  // Flux Table
  string fluxfile;
  
  // Kcorr Table
  string kcorrfile;

  // nVisit
  double *nVisit;

  string CFHTLS_ErrMagFile;
  
  // lib
  string libcat;

  // Prior
  string prior;
  string priorfile;

  // PhotoZ
  double zgrid[3];
  double ebvgrid[3];
  int deltatypegrid;

  string selcut;
  bool fitonlytype;
  
  // PDF
  int savePDF;

  // BDT
  string BDTdir;
  
  string Getbfiltername();
 private:
  string bfiltername;
  bool print;

};
#endif /*TDataCard_H_*/
