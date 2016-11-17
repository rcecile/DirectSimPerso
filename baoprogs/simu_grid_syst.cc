
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

int main(int narg, char* arg[]) {

  cout << " ==== simu_grid_syst.cc program , produce grids and spectra to correct for cell size and sigma effects"<<endl;
  
  // Make sure SOPHYA modules are initialized 
  SophyaInit();  
  InitTim();

  int  cell_test[2];
  int  cell_ref = 1;
  cell_test[0]= 8;
  cell_test[1]= 16;

  double meandens = 1.;

  int rc = 1;  
  try {  // exception handling try bloc at top level
    TArray<r_8> gbase,gcell0,gcell1; 
    int ndim=3;
    int Nxyz = 256;
    sa_size_t mydim[3] = {Nxyz,Nxyz,Nxyz}; // a l'envers
    gbase.SetSize(ndim, mydim);  
    gcell0.SetSize(ndim, mydim);  
    gcell1.SetSize(ndim, mydim);  

    RandomGenerator rg; 
    long seed=1;
    uint_8 npoiss;
    rg.SetSeed(seed);
    for(int i=0; i<Nxyz; i++) { 
      for(int j=0; j<Nxyz; j++) { 
	for(int k=0; k<Nxyz; k++)  {
	  npoiss = rg.PoissonAhrens(meandens); // Poisson fluctuate
	  gbase(i,j,k) = (double)npoiss;

	  double cell=0;
	  for (int id=0; id<cell_test[0];  id++) {
	    npoiss = rg.PoissonAhrens(meandens); 
	    cell += (double)npoiss;
	  }
	  cell /= (double)cell_test[0];
	  gcell0(i,j,k) = cell;

	  cell=0;
	  for (int id=0; id<cell_test[1];  id++) {
	    npoiss = rg.PoissonAhrens(meandens); 
	    cell += (double)npoiss;
	  }
	  cell /= (double)cell_test[1];
	  gcell1(i,j,k) = cell;
	}
      }
    }

    for (int i=0; i<12 ;i++) cout << gbase(0,i,123) << " "; cout << endl;
    for (int i=0; i<12 ;i++) cout << gcell0(0,i,123) << " "; cout << endl;
    for (int i=0; i<12 ;i++) cout << gcell1(0,i,123) << " "; cout << endl;
    
    bool doPixCorr = true;
    double photoZerrorMpc = 0;	
    double maxk_in_calc = 1000; // Set maximum radial k in ps calc
    bool doUnDamp = false;	// undamp Fourier components
    int nbin=175;			// Number of k bins in power spectrum
    double kmin = 0.; 
    double ratio_AngDiam = 1.;
    double ps_ref, ps_8Mpc, ps_16Mpc; // sum of Fourier coefficients (for checking)

    double kmax = PI/cell_ref;
    HProf histogram_ref(kmin, kmax, nbin);
    PowerSpec powerSpectrum_ref(gbase,cell_ref,ratio_AngDiam);
    cout << "compute PS "<< endl;
    ps_ref = powerSpectrum_ref.AccumulatePowerSpectra(histogram_ref, 
								 doPixCorr, maxk_in_calc, photoZerrorMpc, doUnDamp);
    kmax = PI/cell_test[0];
    HProf histogram_8Mpc(kmin, kmax, nbin);
    PowerSpec powerSpectrum_8Mpc(gcell0,cell_test[0],ratio_AngDiam);
    cout << "compute PS "<< endl;
    ps_16Mpc = powerSpectrum_8Mpc.AccumulatePowerSpectra(histogram_8Mpc, 
								 doPixCorr, maxk_in_calc, photoZerrorMpc, doUnDamp);
    kmax = PI/cell_test[1];
    HProf histogram_16Mpc(kmin, kmax, nbin);
    PowerSpec powerSpectrum_16Mpc(gcell1,cell_test[1],ratio_AngDiam);
    cout << "compute PS "<< endl;
    ps_16Mpc = powerSpectrum_8Mpc.AccumulatePowerSpectra(histogram_16Mpc, 
								 doPixCorr, maxk_in_calc, photoZerrorMpc, doUnDamp);


    string outfile = "/sps/lsst/data/rcecile/Simu/ref.txt";
    r_4	volcat = gbase.SizeX()*gbase.SizeY()*gbase.SizeZ()*pow(cell_ref,3); 
    cout << "write PS ref"<< endl;
    powerSpectrum_ref.Write1PS(outfile,histogram_ref,volcat);

    outfile = "/sps/lsst/data/rcecile/Simu/test8Mpc.txt";
    volcat = gcell0.SizeX()*gcell0.SizeY()*gcell0.SizeZ()*pow(cell_test[0],3); 
    powerSpectrum_8Mpc.Write1PS(outfile,histogram_8Mpc,volcat);

    outfile = "/sps/lsst/data/rcecile/Simu/test16Mpc.txt";
    volcat = gcell1.SizeX()*gcell1.SizeY()*gcell1.SizeZ()*pow(cell_test[1],3); 
    powerSpectrum_16Mpc.Write1PS(outfile,histogram_16Mpc,volcat);

    cout << "au boulot "<< endl;

  }  // End of try bloc 
  
  catch (PThrowable & exc) {
    // catching SOPHYA exceptions
    cerr << " simu_grid_syst.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  
    // catching standard C++ exceptions
    cerr << " simu_grid_syst.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  
    // catching other exceptions
    cerr << " simu_grid_syst.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of grid_data.cc program  Rc= " << rc << endl;
  return rc;	
}
