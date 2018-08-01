/**
  * @file  computepsfromgrids.cc
  * @brief compute power spectrum from gridded data
  *
  * herited from Alex Abate computepsfromarray
  * but very simplified
  *
  */

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

// sophya
#include "machdefs.h"
#include "sopnamsp.h"
#include "timing.h"
#include "array.h"
#include "hisprof.h"
#include "histerr.h"
#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "swfitsdtable.h"
#include "resusage.h"

// DirectSim
#include "mydefrg.h"
#include "geneutils.h"
#include "cat2grid.h"
#include "mass2gal.h"
#include "pkspectrum.h"
#include "fitkbaoscale.h"
#include "chisqstats.h"
#include "powerspecfromgrids.h"


void usage(void);
void usage(void) {

	cout << endl<<" Usage: computepsfromgrids [...options...]" << endl<<endl;
	
	cout << "  Compute power spectrum from gridded galaxy data. The     "<<endl;
	cout << "  output power spectrum is correctly normalized and the    "<<endl;
	cout << "  distortion in the simulated density distribution (from   "<<endl;
	cout << "  setting over-densities with delta<-1 equal to -1) is     "<<endl;
	cout << "  properly taken account of. This had to be done because   "<<endl;
	cout << "  delta<-1 corresponds to a negative (unphysical) density. "<<endl;
	cout << "  This can be interpreted to arise from structure formation"<<endl;   
	cout << "  on nonlinear scales not included in the simulation method."<<endl;
	cout << endl;
	
	cout << "  The file containing the gridded data for power spectrum  "<<endl;
	cout << "  analysis is supplied with option -C. This file is        "<<endl;
	cout << "  probably output from the subfromfull program.            "<<endl;
	cout << endl;
		
	cout << "  The shot noise power spectrum is computing using    "<<endl;
	cout << "  gridded data made from a random catalog read in from the "<<endl;
	cout << "  same file as the gridded galaxy data. "<<endl;
	cout << endl;
				
	cout << "  This code reads the cosmology in the header of the input catalog."<<endl;
	cout <<endl;
	
	cout << "  EXAMPLE: "<<endl;
	cout << endl;
	
	cout << "  $ computepsfromgrids -C subgrids.fits  -O powerspectra -m 0.5 "<<endl;
	cout << endl;
	
	cout << " -C : infile : file containing gridded data                        "<<endl;
	cout << " -O : outfile : root filename of text file the galaxy power spectra"<<endl;
	cout << "                are written to                                     "<<endl;
	cout << " -w : w0,wa : if dark energy differs from grid header one          "<<endl;
	cout << " -N : NormNgalMean: to take into account negative cells            "<<endl;
	cout << " -P : additionnal tests and prints "<<endl;
	cout << endl;
	}
	
int main(int narg, char* arg[]) {
	cout << " ==== computepsfromgrids_main.cc program , compute power spectrum";
	cout << " from grided data  ==== " <<endl;
	
	
	// Make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	

	// Set defaults etc ....
	// FILES TO READ IN/OUT
	string infile, outfileroot, subinfo;
	bool do_change_w = false;
	double w0_other =-99.;
	double wa_other =-99.;
	double ratio_AngDiam = 1.;
	double NormNgalMean = 1.;
	
	// decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hC:O:w:N:")) != -1) {
	  switch (c) {
	  case 'C' :
	    infile = optarg;
	    break;
	  case 'O' :
	    outfileroot	= optarg;
	    break;
	  case 'w' :
	    sscanf(optarg,"%lf,%lf",&w0_other,&wa_other);
	    do_change_w = true;
	    break;
	  case 'N' :
	    sscanf(optarg,"%lf",&NormNgalMean); // parameter to correct the mean galaxies nb due to negative cells set to 0
	    break;
	  case 'h' :
	  default :
	    usage(); return -1;
	  }
	}
 
	// Command line argument printing
	cout << "     Printing command line arguments ... "<<endl;
	cout << "     Galaxy sub-array file is "<< infile <<endl;
	cout << "     Galaxy power spectrum will be output to "<< outfileroot <<endl;
	cout << endl;

  
	int rc = 1;  
	try {  // exception handling try bloc at top level
	  
	  // monitor memory usage
	  ResourceUsage res;
	  
	  // Read in gridded galaxy data
	  cout <<"     Read in file "<< infile <<endl;
	  FitsInOutFile fin(infile, FitsInOutFile::Fits_RO);  


	  int NbHDU = fin.NbHDUs();
	  fin.MoveAbsToHDU (NbHDU);  // numero de HDU commence a 1  (convention FITS !) 
	  
	  // modif Adeline : read header of grid-data output (header of subgride output is alos modified
	  double z_center = atof(fin.KeyValue("ZREF").c_str()); 
	  double grid_res = atof(fin.KeyValue("DX").c_str());  // RQ : DX = DY = DZ = R
	  double nGalGrid = atof(fin.KeyValue("NWGRID").c_str()); 
	  int   GridThick = atof(fin.KeyValue("NZ").c_str()); 
	  if (atof(fin.KeyValue("NX").c_str()) < GridThick) GridThick = atof(fin.KeyValue("NX").c_str());
	  if (atof(fin.KeyValue("NY").c_str()) < GridThick) GridThick = atof(fin.KeyValue("NY").c_str());
	  // Set cosmology 
	  cout << "     Initialise cosmology:"<<endl;
	  
	  //Modif Adeline : read cosmo parameters in file header
	  string H0_s, OmegaM_s, OmegaL_s, OmegaB_s, OmegaR_s, wDE_s, wDA_s, Sigma8_s, Ns_s, Ngrids_s;
	  double h, OmegaM, OmegaL, OmegaB, OmegaR, wDE, wDA, Sigma8, n_s;
	  int Ngrids;
	  H0_s = fin.KeyValue("H0");
	  OmegaM_s = fin.KeyValue("OMEGAM0");
	  OmegaL_s = fin.KeyValue("OMEGADE0");
	  OmegaB_s = fin.KeyValue("OMEGAB0");
	  OmegaR_s = fin.KeyValue("OMEGAR0");
	  wDE_s = fin.KeyValue("DE_W0");
	  wDA_s = fin.KeyValue("DE_WA");
	  Sigma8_s = fin.KeyValue("SIGMA8");
	  Ns_s = fin.KeyValue("N_S");

	  h = atof(H0_s.c_str()) / 100;
	  OmegaM = atof(OmegaM_s.c_str());
	  OmegaL = atof(OmegaL_s.c_str());
	  OmegaB = atof(OmegaB_s.c_str());
	  OmegaR = atof(OmegaR_s.c_str());
	  wDE = atof(wDE_s.c_str());
	  wDA = atof(wDA_s.c_str());
	  Sigma8 = atof(Sigma8_s.c_str());
	  n_s = atof(Ns_s.c_str());

	  SimpleUniverse su(h, OmegaM, OmegaL);
	  su.SetOmegaBaryon(OmegaB);
	  su.SetOmegaRadiation(OmegaR);
	  su.SetSigma8(Sigma8);
	  su.SetSpectralIndex(n_s);
	  su.SetFlatUniverse_OmegaLambda(); // Cecile modif - to be sure that it is flat by adjusting OmegaLambda
	  cout << "     OmegaK="<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
	  cout << ", OmegaL="<< su.OmegaLambda() <<", OmegaB="<< su.OmegaBaryon();
	  cout << ", Omega_rad=" << su.OmegaRadiation() << ", Omega_cdm=" << su.OmegaCDM() <<", H0="<< su.H0() << endl;
	  cout << "check flatness: OmegaTot=" << su.OmegaTotal() << endl;
	  if (wDE != -1 or wDA !=0)  
	    su.SetDarkEnergy(su.OmegaLambda(),wDE,wDA);
	  
	  cout << " and w0=" << su.wDE() << ", wA=" << su.waDE() << ", sigma8=" << su.Sigma8() << endl;
	  cout << "Spectral index=" << su.Ns() << endl;
	  cout << "_____________________________________________________________________________________"<< endl;
	  cout << endl;

	  ///////////////////////////////////////////////////////////////////////////////////////////////////////
	  Ngrids_s = fin.KeyValue("NumberOfGrids");
	  Ngrids = atof(Ngrids_s.c_str());
	  
	  string outfile;
	  double meangw, siggw, meangr, siggr;
	  r_4 volcat;
	  double nGalGrid_tot;

	  double sum_FourierCoeffs; // sum of Fourier coefficients (for checking)
	  // POWER SPECTRUM COMPUTATION PARAMETERS
	  double kmax = PI/grid_res;
	  double kmin = PI/grid_res/GridThick; 
	  int nbin=GridThick/2.; // so = kmax / kmin as kmin = deltak
	  cout << "Compute PS with in [" << kmin << " - " << kmax << "] with nbin = "<< nbin << endl;

	  HProf histogram_weighted(kmin, kmax, nbin);
	  HProf histogram_weighted_tot(kmin, kmax, nbin);
	  HProf histogram_random(kmin, kmax, nbin);
	  HProf histogram_random_tot(kmin, kmax, nbin);

	  if (do_change_w) {
	    su.SetEmissionRedShift(z_center);
	    double add_ref = su.AngularDiameterDistance();
	    su.SetDarkEnergy(su.OmegaLambda(),w0_other,wa_other);
	    su.SetEmissionRedShift(z_center);
	    double add_new =  su.AngularDiameterDistance();
	    ratio_AngDiam = add_new / add_ref;
	    cout << "Ratio of angular diameters is " << ratio_AngDiam << " for w0,wa = "<< w0_other <<  ",  " << wa_other << endl;
	    su.SetDarkEnergy(su.OmegaLambda(),su.wDE(),su.waDE());
	  } 
	  else cout << "Dark energy parametrization of the grid FITS header kept"<< endl;
	  
	  cout << "_____________________________________________________________________________________"<< endl;
	  
	  if (NbHDU > 1) {	    

	    for (int igd=0; igd< Ngrids;igd ++) {
	    
	      cout << "____ GRID # " << igd << "_______________________________________________________________________"<< endl;
	    
	    
	      // read wngals ///////////////////////////////////////////////////////////////////////////
	  
	      fin.MoveAbsToHDU (igd+1);
	      TArray<r_8> wngals;
	      fin >> wngals; // in the case where there is no selection effects on the galaxy catalog: ngals=wngals	  
	      
	      // Read data from file header
	      sa_size_t nx = wngals.SizeX(); 
	      sa_size_t ny = wngals.SizeY(); 
	      sa_size_t nz = wngals.SizeZ(); 
	      
	      cout << "    Size of sub-array Nx,Ny,Nz = "<< nx <<","<< ny <<","<< nz;
	      cout <<", resolution = "<< grid_res;
	      cout <<" number of galaxies in grid = "<< nGalGrid << endl;
	      
	      RandomGenerator rg; // need this for cat2grid
	      
	      // grid of galaxies nG is replaced by (nG - <nG>)/<nG> whith <nG> corrected for cells set to 0 by NormNgalMean
	      MeanSigma(wngals, meangw, siggw);
	      cout << "    Input weighted galaxy grid: Mean="<< meangw <<", Sigma="<< siggw <<endl;
	      double MeanNgBar;
	      MeanNgBar = meangw / NormNgalMean ;
	      wngals -= MeanNgBar ;
	      wngals /= MeanNgBar ;
	      
	      MeanSigma(wngals, meangw, siggw);
	      cout << "   Normalized weighted galaxy grid: Mean="<< meangw <<", Sigma="<< siggw <<endl;
	      cout << endl;
	      
	      volcat = wngals.SizeX()*wngals.SizeY()*wngals.SizeZ()*pow(grid_res,3); 
	      cout <<"    Grid volume = "<< volcat <<endl<<endl; 
	      
	      // Compute power spectrum
	      cout <<"     Compute power spectrum of gridded galaxy data "<<endl;
	      
	      PowerSpecFromGrids powerSpectrum_weighted(wngals,grid_res,ratio_AngDiam); // does FT in constructor
	      powerSpectrum_weighted.Setzc(z_center);
	      
	      cout << "AccumulatePowerSpectra for weighted gal-grid"<< endl;
	      sum_FourierCoeffs = 
		powerSpectrum_weighted.AccumulatePowerSpectra(histogram_weighted);
	      
	      cout <<"     Check: sum of Fourier coefficients = "<< sum_FourierCoeffs <<endl;
	      cout <<"            variance of real space field / 2 = "<< siggw*siggw/2 <<endl;
	      
	      
	      histogram_weighted_tot += histogram_weighted;
	      nGalGrid_tot += nGalGrid;

	      // read wrngals ///////////////////////////////////////////////////////////////////////////
	      // Compute shot noise power spectrum
	      cout <<"     Compute shot noise power spectrum from random catalog grid"<<endl;
	      fin.MoveAbsToHDU (Ngrids+igd+1);  // numero de HDU commence a 1  (convention FITS !) 
	      TArray<r_8> wrgals;
	      
	      fin >> wrgals;
	      
	      // grid of Poisson nP is replaced by (nP - <nP>)/<nP> whith <nP> corrected for cells set to 0 by NormNgalMean
	      MeanSigma(wrgals, meangr, siggr);
	      cout << "    Input random galaxy grid: Mean="<< meangr <<", Sigma="<< siggr <<endl;
	      MeanNgBar = meangr / NormNgalMean ;
	      wrgals -= MeanNgBar ;
	      wrgals /= MeanNgBar ;
	      
	      MeanSigma(wrgals, meangr, siggr);
	      cout << "    Normalized random galaxy grid: Mean="<< meangr <<", Sigma="<< siggr <<endl;
	      
	      PowerSpecFromGrids powerSpectrum_random(wrgals,grid_res,ratio_AngDiam);
	      powerSpectrum_random.Setzc(z_center);
	      
	      cout << "AccumulatePowerSpectra for random grid"<<endl;
	      sum_FourierCoeffs = powerSpectrum_random.AccumulatePowerSpectra(histogram_random);
	      histogram_random_tot += histogram_random;
	      
	      cout <<"     Check: sum of Fourier coefficients = "<< sum_FourierCoeffs <<endl;
	      cout <<"            variance of real space field / 2 = "<< siggr*siggr/2 <<endl;
	    
	      cout << endl;
	      
	      // Write out power spectrum
	      string sigd = static_cast<ostringstream*>( &(ostringstream() << igd) )->str();
	      outfile = outfileroot + "_G"+ sigd + "_wngal.txt";
	      
	      powerSpectrum_weighted.WritePS(outfile,histogram_weighted,volcat,histogram_random,su.H0()/100.,nGalGrid);
	      cout << "_____________________________________________________________________________________"<< endl;
	      cout << "_____________________________________________________________________________________"<< endl;
	      
	      // at the end, save the mean PS with extrapolated shot-noise value
	      if (igd == Ngrids-1) { 
		histogram_random_tot /= Ngrids;
		histogram_weighted_tot /= Ngrids;
		outfile = outfileroot + "_wngal.txt";
		powerSpectrum_weighted.WritePS(outfile,histogram_weighted_tot,volcat,histogram_random_tot,su.H0()/100.,nGalGrid_tot,Ngrids,true);
	      }
	    }    
	  }
	   
	} // End of try bloc 
	
	
	catch (PThrowable & exc) {  // catching SOPHYA exceptions
	  cerr << " computepsfromgrids.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
	       << "\n...exc.Msg= " << exc.Msg() << endl;
	  rc = 96;
	}
	catch (std::exception & e) {  // catching standard C++ exceptions
	  cerr << " computepsfromgrids.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
	  rc = 98;
	}
	catch (...) {  // catching other exceptions
	  cerr << " computepsfromgrids.cc: some other exception (...) was caught ! " << endl;
	  rc = 97;
	}
	cout << " ==== End of computepsfromgrids.cc program  Rc= " << rc << endl;
	return rc;	
}
