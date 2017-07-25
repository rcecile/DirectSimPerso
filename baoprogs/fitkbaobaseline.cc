
/**
  * @file  fitkbaobaseline.cc
  * @brief Given an input power spectrum + errors fit the BAO scale to "wiggles only"
  *        power spectrum
  *
  *
  * @author Alex Abate
  * Contact: abate@email.arizona.edu
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <functional>
#include <numeric>
#include <algorithm>

// sophya
#include "sopnamsp.h"
#include "histinit.h"
#include "hisprof.h"
#include "histerr.h"
#include "histos.h"
#include "datatable.h"
#include "fitshdtable.h"
#include "swfitsdtable.h"
#include "fitsarrhand.h"
#include "fiosinit.h"
#include "tarray.h"
#include "datacards.h"

// DirectSim
#include "geneutils.h"
#include "cosmocalcs.h"
#include "pkspectrum.h"
#include "fitkbaobaselinescale.h"




void usage(void);
void usage(void) {

	cout << endl<<" Usage: fitkbaobaselinescale [...options...]              " << endl<<endl;
	
	cout << "  Given an input power spectrum + errors fit the BAO scale."<<endl;
	cout << endl;
 
	cout << "  The input power spectrum has already been corrected for  "<<endl;
	cout << "  shot noise, selection, photo-z etc and is supplied to the"<<endl;
	cout << "  program with option -P                                   "<<endl;
	cout << endl;

	cout << "  Method: divide observed power spectrum by a reference    "<<endl;
	cout << "  power spectrum and fit a sine wave described by some     "<<endl;
	cout << "  amplitude, and a characteristic scale. To compute the    "<<endl;
	cout << "  reference power spectrum the redshift of the observed "<<endl;
	cout << "  power spectrum must be supplied with option -z, and the  "<<endl;
	cout << "  values of sigma_8 and the spectral index parameters must "<<endl;
	cout << "  supplied with option -c "<<endl;
	cout << endl;
	
	cout << "  Results are written to files starting with the root name "<<endl;
	cout << "  supplied with the -O option. The chi-square values,      "<<endl;
	cout << "  reference power spectrum, best-fit sinusoid and fit      "<<endl;
	cout << "  results are written to files "<<endl;
	cout << endl;
	
	cout << " -P : PSFile: power spectrum file to read in               "<<endl;
	cout << "              (4 columns: k (Mpc^-1), P(k) Mpc^3, SN, err)     "<<endl;
	cout << " -O : outfile_root: file root name to write results to     "<<endl; 
	cout << " -s : if input file is simulation (no shotnoise, no sigma) "<<endl; 
	cout << endl;
}



int main(int narg, char *arg[]) {
  
  SophyaInit();
  FitsIOServerInit();
  
  // input power spectrum
  string ps_file, cosmo_file, ref_file;
  // output file
  string outfile_root;
  double maxk  = 0.2;
  double mink = 0.02;
  double Sigma8, n_s;
  bool simu_mode = false;
  
  //--- decoding command line arguments 
  char c;
  while ((c = getopt(narg,arg,"hsP:O:k:d")) != -1) {
    switch (c) {
    case 'P' :
      ps_file = optarg;
      break;
    case 'O' :
      outfile_root = optarg;
      break;
    case 's' :
      simu_mode = true;
      break;
    case 'k' :
      sscanf(optarg,"%lf%lf",&mink,&maxk);
      break;
    case 'h' :
    default :
      usage(); return -1;
    }
  }

  /*cout << "     Printing command line arguments ... "<<endl<<endl;
  cout << "     Reading in observed power spectrum from: "<< ps_file <<endl;
  if (simu_mode) cout << " Will read shot noise and sigmaP from the input power spectrum file" << endl;
  cout << "     Saving results to files beginning "<< outfile_root <<endl;
  cout <<endl;*/
  
    try {
	
      if (maxk <= 0.02) {
	cout << "Value of maxk must be larger than the hardcoded minimum value of k, which is 0.02 !"<< endl;
	exit(-1);
      }
      
      double h;
      TArray<r_8> power_spectrum;
      ifstream ifs(ps_file.c_str());
      sa_size_t nr, nc;

      if (simu_mode){
	//Read power spectrum file
	h=0.679;
	cout << "my h "<< h << endl;

	power_spectrum.ReadASCII(ifs,nr,nc);

      }
      else{
	//Read power spectrum file
	DataCards dc(ps_file);
	h = dc.DParam("HUBBLE");
	cout << "my h "<< h << endl;

	power_spectrum.ReadASCII(ifs,nr,nc,'@');
      }


      // Initialise FitBAOScale
      int degmax=7;  //max degree of the polynomial used to create the baseline
      
      FitBAObaselineScale fitbao( power_spectrum , simu_mode, maxk, mink, h);
      cout<<endl<<endl<<"-------------------------------------------------------------------------"<<endl;
      cout<<"    mink = "<<mink<<"                maxk = "<<maxk<<endl;
      cout<<"-------------------------------------------------------------------------"<<endl;
      fitbao.ComputeChisq(degmax);
      
      
      // print info to a file
      cout << "     Print chisq and results to files"<<endl;
      string outfile;
      outfile = outfile_root + "_chisq.txt";
      cout << "     Write chi^2 to file "<< outfile <<endl;
      //fitbao.WriteChisq(outfile);
      
      outfile = outfile_root + "_result.txt";
      cout << "     Write results to file "<< outfile <<endl;
      fitbao.WriteResults(outfile);
      cout << endl;
      
      outfile = outfile_root + "_baseline.txt";
      cout << "     Write  --k----P(k)----Psmooth(k)--  to file :"<< outfile <<endl;
      fitbao.WriteTable(outfile);
      cout << endl;

      outfile = outfile_root + "_avrg.txt";
      cout << "     Write  --k----AvrgPoints(*k^2)--  to file :"<< outfile <<endl;
      fitbao.WriteTableAvrg(outfile);
      cout << endl;
      
    }// end of try
    
    
    catch(PThrowable exc ) {
    cerr << "fitkbaobaseline.cc , Catched exception: \n" << exc.what() << endl;
    }
    catch(std::exception ex) {
      cerr << "fitkbaobaseline.cc , Catched exception ! " << (string)(ex.what()) << endl;
    }
    catch(...) {
      cerr << "fitkbaobaseline.cc , Catched ... ! " << endl;
    }
    
    cout << "--------------- fitkbao.cc / END --------------------- " << endl;
}// end of main
