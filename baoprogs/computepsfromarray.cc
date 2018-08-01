/**
  * @file  computepsfromarray.cc
  * @brief compute power spectrum from gridded data
  *
  *
  * @author Alex Abate
  * Contact: abate@email.arizona.edu
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
#include "powerspec.h"
#include "mass2gal.h"
#include "pkspectrum.h"
#include "fitkbaoscale.h"
#include "chisqstats.h"


void usage(void);
void usage(void) {

	cout << endl<<" Usage: computepsfromarray [...options...]" << endl<<endl;
	
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
	
	cout << "  In order to do the correction described above, either the"<<endl;
	cout << "  power spectra of the undistorted and distorted over-     "<<endl;
	cout << "  density distribution must be supplied to the program, or "<<endl;
	cout << "  the original over-density distribution itself (so those  "<<endl;
	cout << "  power spectra can be computed here). The file containing "<<endl;
	cout << "  the density distribution or its power spectra are        "<<endl;
	cout << "  supplied with option -S. The density distribution file   "<<endl;
	cout << "  is probably output from the simdensity program.          "<<endl;
	cout << endl;
	
	cout << "  The shot noise power spectrum is also computing using    "<<endl;
	cout << "  gridded data made from a random catalog read in from the "<<endl;
	cout << "  same file as the gridded galaxy data. "<<endl;
	cout << endl;
	
	cout << "  The mean density of the over-density distribution is     "<<endl;
	cout << "  needed to properly normalized the power spectrum. It is  "<<endl;
	cout << "  either read from the file header, or using option -a it  "<<endl;
	cout << "  is passed to the program as an argument (overriding any  "<<endl;
	cout << "  value in the file header).                               "<<endl;
	cout << endl;
	
	cout << "  If the galaxies have approximately Gaussian photo-z errors"<<endl;
	cout << "  the magnitude of this error sigma_z (in format sigma_z*1+z)"<<endl;
	cout << "  should be supplied to the program with the -e option. Then"<<endl;
	cout << "  the power spectrum can be undampled accordingly. To turn  "<<endl;
	cout << "  off the undamping even if the photo-z error is non-zero  "<<endl;
	cout << "  use option -d.                                           "<<endl;    
	cout << endl;
	
	cout << "  The maximum k to use in the power spectrum analysis is set"<<endl;
	cout << "  with option -m.                                          "<<endl;
	cout << endl;
	
				
	cout << "  This code reads the cosmology in the header or the input catalog."<<endl;
	cout <<endl;
	
	cout << "  EXAMPLE: "<<endl;
	cout << endl;
	
	cout << "  $ computepsfromarray -C subgrids.fits -S overdensity.fits "<<endl;
	cout << "                       -O powerspectra -d -m 0.5            "<<endl;
	cout << endl;
	
	cout << " -C : infile : file containing gridded data                        "<<endl;
	cout << " -S : overdensityfile : file containing over-density distribution  "<<endl;
	cout << "                        or over-density power spectra              "<<endl;
	cout << " -O : outfile : root filename of text file the galaxy power spectra"<<endl;
	cout << "                are written to                                     "<<endl;
	cout << " -a : meandens : specify mean density of overdensity distribution  "<<endl;
	cout << " -k : compute power spectra as functions of k_z and k_xy NOT YES OPERATIONAL          "<<endl;
	cout << " -d : Don't undamp photo-z error damping of Fourier coefficients   "<<endl;
	cout << " -e : photoZerror : size of photometric redshift error             "<<endl;
	cout << " -m : maxk_in_calc : maximum kradial used in power spectrum comp   "<<endl;
	cout << " -t : tol_corr; min damoing value in k// -> max_corr=1/tol_corr    "<<endl;
	cout << " -x : doPixCorr : turn off pixel shape correction                  "<<endl;
	cout << " -w : w0,wa : if dark energy differs from grid header one          "<<endl;
	cout << " -N : NormNgalMean: to take into account negative cells            "<<endl;
	cout << " -P : additionnal tests and prints "<<endl;
	cout << endl;
	}
	
int main(int narg, char* arg[]) {
	cout << " ==== computepsfromarray_main.cc program , compute power spectrum";
	cout << " from grided data  ==== " <<endl;
	
	
	// Make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	

	// Set defaults etc ....
	// FILES TO READ IN/OUT
	string infile, overdensityfile, outfileroot, subinfo;
	// CATALOG-TYPE AND REDSHIFT-TYPE PARAMETERS
	double photoZerror = 0;		// Photo-z error size is 0
	// IF HAVE PHOTO-Z ERR>0 PARAMETERS
	//double coeff = 1;		// keep all kradial <= coeff / sig_r
	bool doUnDamp = true;	// undamp Fourier components
	// POWER SPECTRUM COMPUTATION PARAMETERS
	//	int nbin=175;			// Number of k bins in power spectrum
	int nbin=80;			// Number of k bins in power spectrum
	bool doPixCorr = true;	// Correct for pixel size smoothing
	bool doSpectra_z_xy = false;  // True if you want to also compute k_perp= k_xy and k_parall=k_z
	double maxk_in_calc = 1000; // Set maximum radial k in ps calc
	// bool setMaxK=false;		// Set max k separately from z error // Useless, Cecile
	double tol_corr=0.15;          // minimum value of damping in k// for correction 
	double meandens;		// mean density of simlss grid after delta<-1 =-1
	bool do_change_w = false;
	double w0_other =-99.;
	double wa_other =-99.;
	double ratio_AngDiam = 1.;
	bool isMeanDensitySpec = false;
	double NormNgalMean = 1.;
	
	// decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hdxkC:S:O:E:e:w:N:m:a:t:")) != -1) {
	   switch (c) {
	  	    case 'C' :
			    infile = optarg;
		        break;
	        case 'S' :
		        overdensityfile = optarg;
		        break;
	        case 'O' :
		        outfileroot	= optarg;
		        break;
	        case 'd' :
		        doUnDamp = false; // don't undamp Fourier components
		        break;
	        case 'a' :
		        sscanf(optarg,"%lf",&meandens);
		        isMeanDensitySpec=true;
		        break;
	        case 'E' :
		        sscanf(optarg,"%lf",&photoZerror);
		        break;
	        case 'e' :
		        sscanf(optarg,"%lf",&photoZerror);
			cout << "Undamping error on redshift works properly only if zaxis is radial"<< endl;
		        return -1;
	        case 'w' :
		        sscanf(optarg,"%lf,%lf",&w0_other,&wa_other);
			do_change_w = true;
		        break;
	        case 'N' :
		        sscanf(optarg,"%lf",&NormNgalMean); // parameter to correct the mean galaxies nb due to negative cells set to 0
		        break;
	        case 't' :
       		        tol_corr=atof(optarg);
		        break;
	        case 'x' :
		        doPixCorr = false; 
		        break;
	        case 'k' :
		        doSpectra_z_xy = true;
		        break;
	        case 'm' :
		        sscanf(optarg,"%lf",&maxk_in_calc);
		        break;
	        case 'h' :
	            default :
		        usage(); return -1;
	        }
	}


 
	// Command line argument printing
	cout << "     Printing command line arguments ... "<<endl;
	cout << "     Galaxy sub-array file is "<< infile <<endl;

	cout << "     Maximum k_radial in analysis given by "<< maxk_in_calc <<endl;
	if (photoZerror>0) {
		if(doUnDamp)
			cout << "     Undamping Fourier coefficients"<<endl;
		else
			cout << "     NOT undamping Fourier coefficients"<<endl;
		}
	if (photoZerror)
		cout << "     Photo-z error = "<< photoZerror <<endl;
	if(doPixCorr)
		cout << "     Pixel correction is ON"<<endl;
	else
		cout << "     Pixel correction is OFF"<<endl;
	cout << "     Reading simulated over-density grid from file "<< overdensityfile <<endl;
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
	  r_4 volcat, volsim;
	  double photoZerrorMpc = 0;

	  double sum_FourierCoeffs; // sum of Fourier coefficients (for checking)
	  double kmax = PI/grid_res;
	  double kmin = 0.; 
	  
	  HProf histogram_weighted(kmin, kmax, nbin);
	  HProf histogram_random(kmin, kmax, nbin);
	  HProf histogram_weighted_z(kmin, kmax, nbin);
	  HProf histogram_random_z(kmin, kmax, nbin);
	  HProf histogram_weighted_xy(kmin, kmax, nbin);
	  HProf histogram_random_xy(kmin, kmax, nbin);
	  HProf hp(kmin, kmax, nbin);
	  HProf hpf(kmin, kmax, nbin);
	  HProf hp_z(kmin, kmax, nbin);
	  HProf hpf_z(kmin, kmax, nbin);
	  HProf hp_xy(kmin, kmax, nbin);
	  HProf hpf_xy(kmin, kmax, nbin);
	  
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
	  
	  // read ngals ///////////////////////////////////////////////////////////////////////////
	  if (NbHDU > 1) {	    
	    
	    cout << "_____________________________________________________________________________________"<< endl;
	    
	    // read wrngals ///////////////////////////////////////////////////////////////////////////
	    fin.MoveAbsToHDU (Ngrids+1);  // numero de HDU commence a 1  (convention FITS !) 
	    // Compute shot noise power spectrum
	    cout <<"     Compute shot noise power spectrum from random catalog grid"<<endl;
	    {
	      TArray<r_8> wrgals;
	      
	      fin >> wrgals;
	      
	      // grid of Poisson nP is replaced by (nP - <nP>)/<nP> whith <nP> corrected for cells set to 0 by NormNgalMean
	      MeanSigma(wrgals, meangr, siggr);
	      cout << "    Input random galaxy grid: Mean="<< meangr <<", Sigma="<< siggr <<endl;
	      double MeanNgBar;
	      MeanNgBar = meangr / NormNgalMean ;
	      wrgals -= MeanNgBar ;
	      wrgals /= MeanNgBar ;

	      MeanSigma(wrgals, meangr, siggr);
	      cout << "    Normalized random galaxy grid: Mean="<< meangr <<", Sigma="<< siggr <<endl;

	      PowerSpec powerSpectrum_random(wrgals,grid_res,ratio_AngDiam);
	      powerSpectrum_random.Setzc(z_center);
	      
	      cout << "AccumulatePowerSpectra for random grid with maxk_in_calc="<<maxk_in_calc<<" tol_corr="<<tol_corr<<endl;
	      sum_FourierCoeffs = powerSpectrum_random.AccumulatePowerSpectra(histogram_random, 
									      doPixCorr, maxk_in_calc, photoZerrorMpc, doUnDamp, tol_corr);
	      
	      cout <<"     Check: sum of Fourier coefficients = "<< sum_FourierCoeffs <<endl;
	      cout <<"            variance of real space field / 2 = "<< siggr*siggr/2 <<endl;
	    }
	    cout << endl;
	  	  
	    cout << "_____________________________________________________________________________________"<< endl;
	  }
	  // read wngals ///////////////////////////////////////////////////////////////////////////
	  
	  for (int igd=0; igd< Ngrids;igd ++) {
	    
	    fin.MoveAbsToHDU (igd+1);
	    TArray<r_8> wngals;
	    fin >> wngals; // in the case where there is no selection effects on the galaxy catalog: ngals=wngals	  
	    
	    // Read data from file header
	    sa_size_t nx = wngals.SizeX(); 
	    sa_size_t ny = wngals.SizeY(); 
	    sa_size_t nz = wngals.SizeZ(); 
	    
	    
	    if(!isMeanDensitySpec)
	      meandens=atof(fin.KeyValue("MeanOverDensity").c_str());
	    cout << "    Size of sub-array Nx,Ny,Nz = "<< nx <<","<< ny <<","<< nz;
	    cout <<", resolution = "<< grid_res;
	    cout <<" Mpc, mean density of distorted overdensity grid = "<< meandens <<endl;
	    cout <<" number of galaxies in grid = "<< nGalGrid << endl;
	    
	    RandomGenerator rg; // need this for cat2grid
	    
	    // Calculate photo-z error if needed
	    if (photoZerror>0) {
	      cout << "     Calculate comoving photo-z error:"<<endl;
	      photoZerrorMpc = ZErr2CoDistErr(su,photoZerror,z_center); 
	      cout <<"    Redshift of array center = "<< z_center <<", photo-z error (redshift or Zaxis) = ";
	      cout << photoZerror <<", photo-z co-distance error = "<< photoZerrorMpc;
	      cout <<" Mpc"<<endl;
	      cout <<endl;
	    }

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
	    
	    // over-density grid read in - compute its power spectrum
	    
	    cout << "     Compute over-density power spectra "<<endl;
	    cout << "     i) with unmodified density distribution "<<endl;
	    cout << "    ii) when grid cells of delta<-1 were set to delta=-1"<<endl;
	    
	    cout << "     Reading over-density file " << overdensityfile << endl<<endl;  
	    FitsInOutFile fsin(overdensityfile, FitsInOutFile::Fits_RO);
	    TArray<r_8> drho; 
	    fsin >> drho; // delta distribution
	      
	    cout << " test drho.size() "<< drho.SizeX()  <<endl;
	    double min, max;
	    drho.MinMax(min, max);
	    cout <<"     min = "<< min <<", max = "<< max <<endl;	    
	    
	    int xplanes=1;
	    Mass2Gal m2g(drho, su, rg, xplanes);
	    cout << endl;
	    m2g.ReadHeader(fsin);
	    // double zc = m2g.ReturnZref();
	    cout << endl;
	    TArray<r_8> dens, densf;
	    
	    cout <<"     Read out density distribution ... "<<endl;
	    m2g.MassArray(dens);
	    cout << endl;
	    
	    cout <<"     Zero size of original over-density array read in"<<endl;
	    drho.ZeroSize();
	    cout << endl;
	      
	    cout <<"     Get new array with grid cells of delta<-1 were set to delta=-1"<<endl;
	    cout <<"     Mean and Sigma of density fluctuations BEFORE ..."<<endl;
	    double meanm, sigm;
	    MeanSigma(dens, meanm, sigm);
	    cout <<"    Mean="<< meanm <<", Sigma="<< sigm <<", Variance="<< sigm*sigm <<endl;
	    m2g.CleanNegativeMassCells();	      
	    cout <<"     Read out distorted density distribution ... "<<endl;
	    m2g.MassArray(densf);    
	    volsim = m2g.ReturnCubeVol();
	    
	    cout << endl;
	    
	    cout <<"     Mean and Sigma of density fluctuations AFTER ..."<<endl;
	    double meanmf, sigmf;
	    MeanSigma(densf, meanmf, sigmf);
	    cout <<"     Mean="<< meanmf <<", Sigma="<< sigmf <<", Variance="<< sigmf*sigmf <<endl;
	    
	    
	    if ( (meanmf-1)!=meandens ) {
	      cout <<"     mean of fudged density cube (="<< meanmf-1 <<")";
	      cout <<" not equal to stored mean of fudge density cube";
	      cout <<" (="<< meandens <<")"<<endl;
	      cout <<"     probably because mean was never stored in the first";
	      cout <<" place (likely if 2nd value above is 0)"<<endl;
	      if (meandens==0)
		meandens=meanmf-1;
	      cout <<"     Mean density using to correct power spectrum = ";
	      cout << meandens <<endl;
	    }
	    cout <<endl;
	    
	    res.Update();
	    cout << " Memory size (KB):" << res.getMemorySize() << endl;
	    cout << " Resource usage info : \n" << res << endl;
	    
	    cout <<"     Zero size of arrays"<<endl;
	    m2g.ZeroSizeMassArrays();
	    
	    cout <<"     Computing overdensity power spectra"<<endl<<endl;
	    double Dx = m2g.ReturnDX(); double Dy=m2g.ReturnDY(); double Dz=m2g.ReturnDZ();
	    if(Dx-Dy>0||Dz-Dy>0) cout <<"   CAREFUL! Pixels are not cuboid"<<endl;
	    
	    
	    cout <<"     Power spectrum defined using:"<<endl;
	    cout <<"     kmin="<< kmin <<", kmax="<< kmax <<", nbin="<< nbin <<endl;
	    
	    {
	      PowerSpec powerSpectrum_overdensity(dens,Dx,Dy,Dz,ratio_AngDiam);
	      sum_FourierCoeffs = powerSpectrum_overdensity.AccumulatePowerSpectra(hp,doPixCorr);
	      cout <<"     Check: sum of Fourier coefficients = "<< sum_FourierCoeffs <<endl;
	      cout <<"            variance of real space field / 2 = "<< sigm*sigm/2 <<endl;
	    }
	    {
	      PowerSpec powerSpectrum_overdensityf(densf,Dx,Dy,Dz,ratio_AngDiam);
	      sum_FourierCoeffs = powerSpectrum_overdensityf.AccumulatePowerSpectra(hpf, doPixCorr);
	      cout <<"     Check: sum of Fourier coefficients = "<< sum_FourierCoeffs <<endl;
	      cout <<"            variance of real space field / 2 = "<< sigmf*sigmf/2 <<endl;
	    }
	    
	    
	    cout <<"     Overdensity cube volume="<< volsim <<" Mpc";
	    cout <<", survey volume="<< volcat <<" Mpc"<<endl<<endl;
	    
	    // Compute power spectrum
	    cout <<"     Compute power spectrum of gridded galaxy data "<<endl;
	    
	    PowerSpec powerSpectrum_weighted(wngals,grid_res,ratio_AngDiam); // does FT in constructor
	    powerSpectrum_weighted.Setzc(z_center);
	    
	    cout << "AccumulatePowerSpectra for weighted gal-grid with maxk_in_calc="<<maxk_in_calc <<" tol_corr="<<tol_corr<<endl;
	    sum_FourierCoeffs = 
	      powerSpectrum_weighted.AccumulatePowerSpectra(histogram_weighted, 
							    doPixCorr, maxk_in_calc, photoZerrorMpc, doUnDamp, tol_corr);
	    
	    cout <<"     Check: sum of Fourier coefficients = "<< sum_FourierCoeffs <<endl;
	    cout <<"            variance of real space field / 2 = "<< siggw*siggw/2 <<endl;
	      
	    // Write out power spectrum
	    // Write out shot noise power spectrum not done anymore as main information in coloumn 6 of _wngal file
	    string sigd = static_cast<ostringstream*>( &(ostringstream() << igd) )->str();
	    outfile = outfileroot + "_G"+ sigd + "_wngal.txt";
	    
	    powerSpectrum_weighted.WritePS(outfile,histogram_weighted, 
					   volcat,hp,hpf,volsim, histogram_random, meandens, nGalGrid);
	    cout << "_____________________________________________________________________________________"<< endl;
	  }
	   
	} // End of try bloc 
	
	
	catch (PThrowable & exc) {  // catching SOPHYA exceptions
	  cerr << " computepsfromarray.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
	       << "\n...exc.Msg= " << exc.Msg() << endl;
	  rc = 96;
	}
	catch (std::exception & e) {  // catching standard C++ exceptions
	  cerr << " computepsfromarray.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
	  rc = 98;
	}
	catch (...) {  // catching other exceptions
	  cerr << " computepsfromarray.cc: some other exception (...) was caught ! " << endl;
	  rc = 97;
	}
	cout << " ==== End of computepsfromarray.cc program  Rc= " << rc << endl;
	return rc;	
}
