/**
 * @file  rdlss.cc
 * @brief Read over-density grid and simulate galaxy catalog
 *
 * @todo implement LF galaxy simulation properly within Mass2Gal and this program
 *       And use CumulDistM, DrawM classes instead of GalFlxTypDist
 * @todo Use LFParameters class to read in LF parameters
 * @todo why are GalFlxTypDist and extinct arguments to CreateNzHisto?
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu1
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
#include "resusage.h"

// DirectSim
#include "mydefrg.h"
#include "mass2gal.h"
#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "constcosmo.h"

void usage(void);
void usage(void) {

  cout << endl<<" Usage: rdlss [...options...]                  "<<endl<<endl;

  cout << "  Simulate a galaxy catalog from an input over-density     "<<endl;
  cout << "  grid. Over-density grid is generated by simdensity and   "<<endl;
  cout << "  supplied to the program with option -C. The simulation id"<<endl;
  cout << "  is supplied with option -i and is included in each       "<<endl;
  cout << "  galaxy's id in the following way: galaxy '1' will have a "<<endl;
  cout << "  full id number of [id]00,000,000,001, galaxy '2' will    "<<endl;
  cout << "  have a full id number of [id]00,000,000,002. This is so  "<<endl;
  cout << "  galaxies can be traced back to their original over-      "<<endl;
  cout << "  density simulations.                                     "<<endl;
  cout << endl;
        
  cout << "  To stop unnecessary simulation of extremely faint        "<<endl;
  cout << "  galaxies apply a cut on the minimum absolute magnitude   "<<endl;
  cout << "  with option -M. Only applies if simulating full galaxy   "<<endl;
  cout << "  catalog.                                                 "<<endl;
  cout << endl;
        
  cout << "  To simulate catalog within an observation cone, use      "<<endl;
  cout << "  option -a to set the maximum radius in radians a galaxy  "<<endl;
  cout << "  can be from the simulation center.                       "<<endl;
  cout << endl;
        
  cout << "  To simulate a catalog where the z dimension is made to be"<<endl;
  cout << "  the 'radial' direction use option -z.                    "<<endl;
  cout << endl;
        
  cout << "  To not simulate galaxy clustering use option -s to supply"<<endl;
  cout << "  the required galaxy density (number of galaxies per grid "<<endl;
  cout << "  cell).                                                   "<<endl;
  cout << endl;

  cout << "  Option -R ensure that the galaxies are randomly          "<<endl;
  cout << "  positioned within the cells of the over-density grid     "<<endl;
  cout << endl;
        
  // debug
        
  cout << "  The catalog is output to the file specified with option  "<<endl;
  cout << "  -O. There are 4 possible types of outputs:               "<<endl;
  cout << "   type 0: Output full galaxy catalog, either with galaxy  "<<endl; 
  cout << "           clustering or not [DEFAULT]                     "<<endl;
  cout << "   type 1: Output simple catalog: ra,dec,z only [OPTION -S]"<<endl;
  cout << "   type 2: Output catalog of just (true) redshifts         "<<endl;
  cout << "           [OPTION -Z]                                     "<<endl;
  cout << "   type 3: Output a histogram of redshifts [OPTION -H]     "<<endl;
  cout << endl;
        
        
  cout << "    EXAMPLE 1: Simulate a galaxy catalog from an over-     "<<endl;
  cout << "    density distributionin file odens.fits.  We want the   "<<endl;
  cout << "    survey to have a circular sky area with radius = pi/4  "<<endl;
  cout << "    radians, and galaxies with randomized positions in the "<<endl;
  cout << "    pixels. So not to produce too many galaxies so the     "<<endl;
  cout << "    absolute magnitude cut is set. This throws out any     "<<endl;
  cout << "    galaxies with an unobservable absolute magnitude given "<<endl;
  cout << "    their redshifts. This is simulation number 1 so the id "<<endl;
  cout << "    should be 1. The output catalog will be saved to       "<<endl;
  cout << "    out.fits:                                              "<<endl;
  cout << endl;
        
  cout << "    $ rdlss -C odens.fits -O out.fits -a 0.7854 -i 1 -R -M "<<endl;
  cout << endl;
        
  cout << "    EXAMPLE 2: Similar to example 1, but now write only    "<<endl;
  cout << "    the true redshifts to the galaxy catalog and DON'T     "<<endl;
  cout << "    set the absolute magnitude cut:                        "<<endl;
  cout <<endl;
        
  cout << "    $ rdlss -C odens.fits -O ztrue.fits -i 1 -R -Z         "<<endl;
  cout << endl;

  cout << " -O : out_file: filename that the output is written to      "<<endl;
  cout << " -C : over_dens_file: FITS file containing over-density grid"<<endl;
  cout << " -S : Only output a simple catalog of ra,dec,z              "<<endl;
  cout << " -Z : Only output a catalog containing the TRUE redshifts   "<<endl;
  cout << "      (and galaxy IDs)                                      "<<endl;
  cout << " -H : Only produce a histogram of the TRUE catalog n(z)     "<<endl;
  cout << " -i : id: Simulation identifier, integer >=1 [default=1]    "<<endl;
  cout << " -M : Don't add galaxies with M>AbsMagCut(z) [default=YES]  "<<endl; 
  cout << " -G : golden sample cut (m_i<25.3) [default=NO]           "<<endl; 
  cout << " -g : photoZ probability (mre severe than Golden sample)    "<<endl;
  cout << " -R : Randomise galaxy positions within cell [default=NO]   "<<endl;
  cout << " -a : SkyArea: Radius of survey sky area in radians (assumed"<<endl;
  cout << "       circular) [default=2PI]                              "<<endl;
  cout << " -s : ngal: No galaxy clustering in simulation, instead     "<<endl;
  cout << "            constant galaxy density equal to ngal/pixel     "<<endl;
  cout << " -z : z dimension of the over_dens_file is 'radial' [default=NO] "<<endl;
  cout << " -N : debug_out: Output mass cube, ngals cube to files      "<<endl;
  cout <<"       [debug_out]_ngals.fits and [debug_out]_mass.fits      "<<endl;
  cout << endl;
}

int main(int narg, char* arg[]) {

  cout << " ==== rdlss.cc program , simulating galaxy catalog from"<<endl;
  cout << " input over-density grid  ==== " << endl;
        
  // Make sure SOPHYA modules are initialized 
  SophyaInit();  
  FitsIOServerInit();
  InitTim();
  cout<<endl<<endl;


  string infile;            // name of over-density grid to read in
  string outfile;           // name of file to output to
  string debug_out;         // write output to here if debugging 
  double maxRadius = 6.3;   // max radius in radians gal can be from sim center
  bool doSkyAreaCut=false;  // throw out gals with theta (co-tangent wrt z axis)>maxRadius
  bool doVeryFaintCut=false;// cut on min absolute magnitude(z) 
  bool doGoldenCut=false;   // Golden sample cut
  bool doProbaGoldenCut=false;// PhotoZ probability cut
  int out_type = 0;         // output type (0=normal,1=simple,2=z only,3=n(z) histo) 
  //    bool doDebug=false;       // write out ngals, mass cubes
  bool doDebug=false;       // write out ngals, mass cubes
  bool doRandPos=false;     // randomize galaxy positions in cells
  //bool HistoOnly=false; // only histogram of z
  //bool TrueZOnly=false; // only z
  //bool SimpleSim=false; // just ra,dec,z
  bool noClustering=false;  // no clustering
  bool isZRadial=false;     // z dimension IS 'radial' direction
  int idsim = 1;
  float ngal_per_cell = 1;  // density of galaxies when no clustering
  float conv;
  
  string ProbaGoldenCutFileName="";
  
  //--- decoding command line arguments 
  char c;
  while((c = getopt(narg,arg,"hSZHMzRGC:i:O:a:s:N:g:")) != -1) {
    switch (c) {
    case 'C' :
      cout << "CC " << optarg << endl;
      infile = optarg;
      break;
    case 'i' :
      sscanf(optarg,"%d",&idsim);
      break;
    case 'O' :
      outfile   = optarg;
      break;
    case 'G' :
      doVeryFaintCut = true;
      doGoldenCut = true;
      break;
    case 'g' :
      cout << "g case " << optarg << endl;
      ProbaGoldenCutFileName   = optarg;
      doVeryFaintCut = true;
      doProbaGoldenCut = true;
      break;
    case 'S' :
      out_type = 1;
      break;
    case 'Z' :
      out_type = 2;
      break;
    case 'H' :
      out_type = 3;
      break;
    case 'M' :
      doVeryFaintCut = true;
      break;
    case 'a' :
      sscanf(optarg,"%lf",&maxRadius);
      doSkyAreaCut = true;
      break;
    case 'R' :
      doRandPos=true;
      break;            
    case 's' :
      sscanf(optarg,"%f",&ngal_per_cell);
      noClustering=true;
      break;
    case 'z' :
      isZRadial=true;
      break;
    case 'N' :
      doDebug = true;
      debug_out = optarg;
      break;
    case 'h' :
    default :
      usage(); return -1;
    }
  }

                
  cout << "     Over-density grid will be read from "<< infile <<endl;

  cout << "     Output file will be "<< outfile <<":"<<endl;
  switch (out_type) {
  case 0:
    cout << "     Creating full galaxy catalog "<<endl;
    break;
  case 1:
    cout << "     Creating simple, positions only galaxy catalog "<<endl;
    break;
  case 2:
    cout << "     Creating galaxy catalog of redshifts only"<<endl;
    break;
  case 3:
    cout << "     Creating galaxy catalog redshift histogram only"<<endl;
    break;
  default:
    throw ParmError("ERROR! out_type not understood");
    break;
  }
  
  cout << "     Simulation ID number is "<< idsim <<"00,000,000,000"<<endl; 
  if (doVeryFaintCut && out_type==0)
    cout << "     Removing galaxies with very faint absolute magnitudes"<<endl;
  if (doGoldenCut && out_type==0)
    cout << "     Removing galaxies with very faint absolute magnitudes (mi>25.3, golden cut) "<<endl;
  if (doProbaGoldenCut && out_type==0)
    cout << "    Removing galaxies using PhotoZ probabilities " << endl;
  if (isZRadial && out_type==0)
    cout << "     z-dimension is of the input cube is 'radial dimension'"<<endl;
  if (noClustering) {
    cout << "     Catalog will be simulated with a constant number"<<endl;
    cout << " density of "<< ngal_per_cell <<endl;
  }
  if (doSkyAreaCut)
    cout << "     Keeping galaxies within theta<"<< maxRadius <<endl; 
  if (doRandPos)
    cout << "     Randomising galaxy positions within pixel"<<endl;
  //-- end command line arguments
        
  int rc = 1;  
  try {  // exception handling try bloc at top level
  
    
    // Track cpu resource usage
    ResourceUsage res;
    res.Update();
    cout << " Memory size (KB):" << res.getMemorySize() << endl;
    cout << " Max memory size (KB):" << res.getMaxMemorySize() << endl;
    cout << " Maximum allowed data segment size (KB):"<< res.getMaxDataSize() <<endl;
    cout << " Resource usage info : \n" << res << endl;
        
  
    // Read over-density grid
    cout << "     Reading input file= " << infile << endl;  
    FitsInOutFile fin(infile,FitsInOutFile::Fits_RO);
    TArray<r_8> drho;
    fin >> drho;
    cout << drho.Info();
    cout << "     Print original drho array size: "<< drho.SizeX() <<"x";
    cout << drho.SizeY() <<"x"<< drho.SizeZ() <<endl;
    cout << endl;
   
    string NXtrue_s, NXwritten_s;
    int NXtrue, NXwritten;
    NXtrue_s = fin.KeyValue("NZ");
    NXwritten_s = fin.KeyValue("NAXIS1");
    NXtrue = atoi(NXtrue_s.c_str());
    NXwritten = atoi(NXwritten_s.c_str());
    int xplanes = NXwritten - NXtrue; // SimLSS simulates cube with 1 or 2 extra planes: one dimension is too long
    // xplanes should be the difference between NX and drho.SizeX()

   
    // Initialize cosmological parameters 
    cout << "     Initialise cosmology: (same as SimLSS)"<<endl;

    //////////  modif Adeline : read cosmo in Fits_RO header
    string H0_s, OmegaM_s, OmegaL_s, OmegaB_s, OmegaR_s, wDE_s, wDA_s, Sigma8_s, Ns_s;
    double h, OmegaM, OmegaL, OmegaB, OmegaR, wDE, wDA, Sigma8, n_s;
    H0_s = fin.KeyValue("H0");
    OmegaM_s = fin.KeyValue("OMEGAM0");
    OmegaL_s = fin.KeyValue("OMEGADE0");
    OmegaB_s = fin.KeyValue("OMEGAB0");
    OmegaR_s = fin.KeyValue("OMEGAR0");
    wDE_s = fin.KeyValue("DE_W0");
    wDA_s = fin.KeyValue("DE_WA");
    Sigma8_s = fin.KeyValue("SIGMA8");
    Ns_s = fin.KeyValue("N_S");

    h = atof(H0_s.c_str()) / 100.;
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
    // Cecile modif
    if (wDE != -1 or wDA !=0)  
      su.SetDarkEnergy(su.OmegaLambda(),wDE,wDA);
        
    cout << " and w0=" << su.wDE() << ", wA=" << su.waDE() << ", sigma8=" << su.Sigma8() << endl;
    cout << "Spectral index=" << su.Ns() << endl;
    cout << endl;
        
        
    // Initialize Mass2Gal and remove N extra planes 
    cout << "     Initialise Mass2Gal: remove planes" << endl;
    RandomGenerator rg;
    Mass2Gal m2g(drho, su, rg, xplanes, isZRadial);
    double mean, sig;
    TArray<r_8> mass;
    m2g.MassArray(mass);
    MeanSigma(mass, mean, sig);
    cout << endl<<"     RAW DENS CUBE STATS: Mean=" << mean << " Sigma=";
    cout << sig <<endl<<endl;
    res.Update();
    cout << " Memory size (KB):" << res.getMemorySize() << endl;                                          
    cout << " Resource usage info : \n" << res << endl;
        
        
    // Read in cube and pixel properties from fits header
    cout << "     Read in cube properties from fits header" << endl;
    m2g.ReadHeader(fin);
    // xplanes should be the difference between NX and drho.SizeX()
    int NZ=m2g.ReturnNZ();
    int diff = drho.SizeX()-NZ;
    if( xplanes!=abs(diff) ) {
      cout << " rdlss/Error: drho.SizeX()="<<drho.SizeX()<<" - NZ="<<NZ<<" != xplanes="<<xplanes<<endl;
      cout << " ... removed wrong number of planes from SimLSS cube" << endl;
      throw ParmError("ERROR: removed wrong number of planes from SimLSS cube");
    }
    m2g.SetRandPos(doRandPos);


    // Deal with negative mass cells AND drho/rho -> rho/rho^bar 
    cout << "     Clean Negative Mass Cells" << endl;
    sa_size_t nbadc = m2g.CleanNegativeMassCells(); // adds 1 to drho and sets anything <0=0
    m2g.MassArray(mass);
    cout << "     NBadCells=" << nbadc 
         << " BadPercentage=" << (double)nbadc*100./(double)mass.Size() << "%"<<endl;
    MeanSigma(mass, mean, sig);
    cout << endl<<"     CLEANED DENS CUBE STATS: Mean=" << mean << " Sigma=";
    cout << sig << endl;
    // double check there are no bad cells 
    sa_size_t nbadc2 = m2g.CheckNegativeMassCells();
    cout <<"     double check there are no bad cells ...."<< nbadc2 <<endl;
    cout <<"     check minimum value in mass array is 0"<<endl;
    double min, max;
    mass.MinMax(min,max);
    cout <<"     min of mass = "<< min <<", max of mass = "<< max <<endl;
    cout << endl<<endl;
    res.Update();
    cout << " Memory size (KB):" << res.getMemorySize() << endl;
    cout << " Resource usage info : \n" << res << endl;
        
        
    cout << "     Convert rho/rho^bar To Mean NGal"<<endl;  

    // Read in LFs: replace the below with LFParameters class
    cout << "     Set up Schechter functions for each galaxy type"<<endl;
    cout << " ... GOODS B band: Early types, Late types, Starbursts"<<endl;
    cout << " ... see Table 3 in Dahlen et al 2005"<<endl;

    //  --- Modified by Reza  , see below 
    string LFplace;
    char * plf=getenv("SIMBAOLF");
    if (plf==NULL) {
      cout <<"ERROR LF LOCATION ENVIRONMENT VARIABLE NOT DEFINED"<<endl;
      return 1;
    }
    else {
      LFplace=plf;
      cout <<"    Location of LF file is "<< LFplace <<endl;
    }
    string LFfile = LFplace + "GOODS_B_LF.txt";// add an option for this

         
    ifstream ifs;
    ifs.open(LFfile.c_str(), ifstream::in);
    if (ifs.fail())
      cout <<"  ERROR: failed to find luminosity function file"<<endl;
    //ifs.open(LFfile);
    TArray<r_4> LFTable;
    sa_size_t nr, nc;
    LFTable.ReadASCII(ifs,nr,nc);
    cout << " rdlss: from LFTable, nr="<<nr<<" nc="<<nc<<endl;
    //  cout << LFTable ;
        
    int MstarCol=2, AlphaCol=3, PhiStarCol=4;
    // ALL GALAXIES
    /*
      double MstarAz1=LFTable(MstarCol,13),alpAz1=LFTable(AlphaCol,13),
      phistarAz1=LFTable(PhiStarCol,13)*1e-4;
      double MstarAz2=LFTable(MstarCol,14),alpAz2=LFTable(AlphaCol,14),
      phistarAz2=LFTable(PhiStarCol,14)*1e-4;
      double MstarAz3=LFTable(MstarCol,15),alpAz3=LFTable(AlphaCol,15),
      phistarAz3=LFTable(PhiStarCol,15)*1e-4;
        
      // EARLY TYPES
      double MstarEz1=LFTable(MstarCol,1),alpEz1=LFTable(AlphaCol,1),
      phistarEz1=LFTable(PhiStarCol,1)*1e-4;
      double MstarEz2=LFTable(MstarCol,6),alpEz2=LFTable(AlphaCol,6),
      phistarEz2=LFTable(PhiStarCol,6)*1e-4;
      double MstarEz3=LFTable(MstarCol,10),alpEz3=LFTable(AlphaCol,10),
      phistarEz3=LFTable(PhiStarCol,10)*1e-4;
        
      // LATE TYPES
      double MstarLz1=LFTable(MstarCol,3),alpLz1=LFTable(AlphaCol,3),
      phistarLz1=LFTable(PhiStarCol,3)*1e-4;
      double MstarLz2=LFTable(MstarCol,7),alpLz2=LFTable(AlphaCol,7),
      phistarLz2=LFTable(PhiStarCol,7)*1e-4;
      double MstarLz3=LFTable(MstarCol,11),alpLz3=LFTable(AlphaCol,11),
      phistarLz3=LFTable(PhiStarCol,11)*1e-4;
        
      // STARBURST TYPES
      double MstarSz1=LFTable(MstarCol,4),alpSz1=LFTable(AlphaCol,4),
      phistarSz1=LFTable(PhiStarCol,4)*1e-4;
      double MstarSz2=LFTable(MstarCol,8),alpSz2=LFTable(AlphaCol,8),
      phistarSz2=LFTable(PhiStarCol,8)*1e-4;
      double MstarSz3=LFTable(MstarCol,12),alpSz3=LFTable(AlphaCol,12),
      phistarSz3=LFTable(PhiStarCol,12)*1e-4;
        
      /*
      //---- reading and initializing the LF (Schechter) functions 
      string LFfile = "GOODS_B_LF.txt";// add an option for this
      cout << " rdlss: reading LF params from file:"<<LFfile<<endl;
      LFParameters lfpars(LFfile, 1);
      // type=0, All LF; type=1, Early LF; type=2, Late LF; type=3, SB LF
      double MstarAz3, alpAz3, phistarAz3;
      lfpars.ReturnParsBini(MstarAz3, alpAz3,phistarAz3,2,0);
      double MstarEz3, alpEz3, phistarEz3;
      lfpars.ReturnParsBini(MstarEz3, alpEz3,phistarEz3,2,1);
      double MstarLz3, alpLz3, phistarLz3;
      lfpars.ReturnParsBini(MstarLz3, alpLz3,phistarLz3,2,2);
      double MstarSz3, alpSz3, phistarSz3;
      lfpars.ReturnParsBini(MstarSz3, alpSz3,phistarSz3,2,3);
      //----  fin modif Reza 
      */
    ////////////////////// modif Adeline (inversion between colomn and row) //////////////
    int z1row=0, z2row=1, z3row=2;
    // ALL GALAXIES
    double MstarAz1=LFTable(MstarCol,z1row),alpAz1=LFTable(AlphaCol,z1row),
      phistarAz1=LFTable(PhiStarCol,z1row)/**1e-4*/;
    double MstarAz2=LFTable(MstarCol,z2row),alpAz2=LFTable(AlphaCol,z2row),
      phistarAz2=LFTable(PhiStarCol,z2row)/**1e-4*/;
    double MstarAz3=LFTable(MstarCol,z3row),alpAz3=LFTable(AlphaCol,z3row),
      phistarAz3=LFTable(PhiStarCol,z3row)/**1e-4*/;
        
    // EARLY TYPES
    double MstarEz1=LFTable(MstarCol+3,z1row),alpEz1=LFTable(AlphaCol+3,z1row),
      phistarEz1=LFTable(PhiStarCol+3,z1row)/**1e-4*/;
    double MstarEz2=LFTable(MstarCol+3,z2row),alpEz2=LFTable(AlphaCol+3,z2row),
      phistarEz2=LFTable(PhiStarCol+3,z2row)/**1e-4*/;
    double MstarEz3=LFTable(MstarCol+3,z3row),alpEz3=LFTable(AlphaCol+3,z3row),
      phistarEz3=LFTable(PhiStarCol+3,z3row)/**1e-4*/;
        
    // LATE TYPES
    double MstarLz1=LFTable(MstarCol+6,z1row),alpLz1=LFTable(AlphaCol+6,z1row),
      phistarLz1=LFTable(PhiStarCol+6,z1row)/**1e-4*/;
    double MstarLz2=LFTable(MstarCol+6,z2row),alpLz2=LFTable(AlphaCol+6,z2row),
      phistarLz2=LFTable(PhiStarCol+6,z2row)/**1e-4*/;
    double MstarLz3=LFTable(MstarCol+6,z3row),alpLz3=LFTable(AlphaCol+6,z3row),
      phistarLz3=LFTable(PhiStarCol+6,z3row)/**1e-4*/;
        
    // STARBURST TYPES
    double MstarSz1=LFTable(MstarCol+9,z1row),alpSz1=LFTable(AlphaCol+9,z1row),
      phistarSz1=LFTable(PhiStarCol+9,z1row)/**1e-4*/;
    double MstarSz2=LFTable(MstarCol+9,z2row),alpSz2=LFTable(AlphaCol+9,z2row),
      phistarSz2=LFTable(PhiStarCol+9,z2row)/**1e-4*/;
    double MstarSz3=LFTable(MstarCol+9,z3row),alpSz3=LFTable(AlphaCol+9,z3row),
      phistarSz3=LFTable(PhiStarCol+9,z3row)/**1e-4*/;
        
    /////////////////////////////////////////////////////////////////////////////////////


      string MstarUnits="M-5log10h70";
      string phistarUnits="(Mpc/h70)^-3";
        
      cout << " ... Schechter parameters"<<endl; 
      cout << "     z range     Mstar     alpha     phistar     spec type"<<endl;
      cout << "    0.75<z<1.0  "<< MstarAz3 <<"     "<< alpAz3 <<"       ";
      cout << phistarAz3 <<"        All"<<endl;
      cout << "                "<< MstarEz3 <<"     "<< alpEz3 <<"       ";
      cout << phistarEz3 <<"         Early"<<endl;
      cout << "                "<< MstarLz3 <<"     "<< alpLz3 <<"       ";
      cout << phistarLz3 <<"          Late"<<endl;
      cout << "                "<< MstarSz3 <<"     "<< alpSz3 <<"        ";
      cout << phistarSz3 <<"        Starburst"<<endl<<endl;

      // Find conversion from mass density to galaxy density
      cout<<"     Mass to Galaxy number conversion"<<endl;
      Schechter schAz3(phistarAz3,MstarAz3,alpAz3);
      //schechter functions for each type
      Schechter schEz3(phistarEz3,MstarEz3,alpEz3);
      Schechter schLz3(phistarLz3,MstarLz3,alpLz3);
      Schechter schSz3(phistarSz3,MstarSz3,alpSz3);
        
      double schmin=-24, schmax=-13;// units of "M-5log10h70"
      int schnpt=10000;
      cout << "    ... integrating from Mbright="<< schmin <<" "<< MstarUnits;
      cout << " to Mfaint="<< schmax <<" "<< MstarUnits <<", with step=";
      cout << schnpt <<endl;
      schAz3.SetInteg(schmin,schmax,schnpt);
      double nz3=schAz3.Integrate();
      cout <<"    ... number density of galaxies: "<< nz3 <<" Mpc^-3"<<endl;
      double pixVol = m2g.ReturnPixVol();
      cout <<"    pixel volume="<< pixVol <<" Mpc^3"<<endl;
      conv = pixVol*nz3;
      cout <<"    actual gals per pixel vol="<< pixVol*nz3 << endl;
      cout <<"    gals per pixel vol="<< conv << endl;
        
        
      // turn off clustering if noClustering set
      if(noClustering)
        conv = ngal_per_cell;
        
        
      // Convert mass in each cell to a (mean) number of galaxies
      m2g.ConvertToMeanNGal(conv, noClustering); // just multiplies mass_ by conv
      res.Update();
      cout << " Memory size (KB):" << res.getMemorySize() << endl;
      cout << " Resource usage info : \n" << res << endl;
        
        
      // Write out ngals and mass cubes for debugging
      if (doDebug) {
                
        cout <<"    **** Writing out cleaned mass distribution AND ngal distribution ****"<<endl;
        TArray<int_4> ngals;
        m2g.NGalArray(ngals);
        TArray<r_8>   mass2;
        m2g.MassArray(mass2);
        cout <<"    check minimum value in mass array is 0"<<endl;
        mass2.MinMax(min,max);
        cout <<"    min of mass = "<<min<<", max of mass = "<<max<<endl;
        string outngal = debug_out +"_ngals.fits";
        FitsInOutFile fos(outngal, FitsInOutFile ::Fits_Create);
        fos << ngals;
        cout <<"    Written ngals array"<<endl;
                
        // have to do below stuff or there is a problem with 
        // writing FT'd array to a FITS file
        // now there shouldn't be a problem because of 
        // doing .PackElements() in constructor
        // TVector<r_8> gridv=m2g.ReturnGridSpec();
        // TArray<r_8> massn;
        // int ndim=3;
        // sa_size_t mydim[ndim];
        // mydim[0]=gridv(0); mydim[1]=gridv(1); mydim[2]=gridv(2);
        // massn.SetSize(ndim, mydim);
        // cout <<"    Size of mass array = "<< mass2.SizeZ() << " " << mass2.SizeY()<< " " << mass2.SizeX()<<endl;
        // for(sa_size_t iz=0; iz<mass2.SizeZ(); iz++) 
        //      for(sa_size_t iy=0; iy<mass2.SizeY(); iy++) 
        //              for(sa_size_t ix=0; ix<mass2.SizeX(); ix++) 
        //                      massn(ix, iy, iz) = mass2(ix, iy, iz);
        // cout <<" here2"<<endl;

        // lignes ci-dessus decommentees
        string outmass = debug_out +"_mass.fits";
        FitsInOutFile fos1(outmass, FitsInOutFile ::Fits_Create);
        fos1 << mass2;
        cout <<"    Written mass array"<<endl;
        cout << endl;
      }
                
      //------------------------------------------------------------------------//
      // -- update to new classes: CumulDistM and DrawM
    
      // Create galaxy distributions
      cout <<"     Set up Mb-Type 2D distribution"<<endl;
      // set distribution parameters
      int magbin=10000;
      int PrtLevel = 2;
      string IntLFUnits="(Mpc/h70)^-3";
      cout <<"    Renormalise type-specific LFs ..."<<endl;
      //get prob distribution
      // - set schechter distributions
      int type;
      type=1;
      SchechterDist schDistEz3(schAz3,schEz3,schLz3,schSz3,type);
      type=2;
      SchechterDist schDistLz3(schAz3,schEz3,schLz3,schSz3,type);
      type=3;
      SchechterDist schDistSz3(schAz3,schEz3,schLz3,schSz3,type);
        
        
      cout <<"     Find fraction of each type ..."<<endl;
      // integrate work out type fractions
      double nfE3=schDistEz3.Integrate(schmin,schmax,schnpt);
      double nfL3=schDistLz3.Integrate(schmin,schmax,schnpt);
      double nfS3=schDistSz3.Integrate(schmin,schmax,schnpt);
      double totalnr=nfE3+nfL3+nfS3;
      double fEz3=nfE3/totalnr;
      double fLz3=nfL3/totalnr;
      double fSz3=nfS3/totalnr;
      cout <<"     Type fractions from renormalised LFs are: fE="<< fEz3;
      cout <<", fL="<< fLz3 <<", fS="<< fSz3 <<endl;
      GalFlxTypDist gfdz3(rg, PrtLevel); 
      gfdz3.AddGalType(schDistEz3, schmin, schmax, fEz3, magbin, schnpt); 
      gfdz3.AddGalType(schDistLz3, schmin, schmax, fLz3, magbin, schnpt); 
      gfdz3.AddGalType(schDistSz3, schmin, schmax, fSz3, magbin, schnpt);
      //------------------------------------------------------------------------//
        
      // Apply magnitude cut if required
      if (doProbaGoldenCut)
	m2g.MaxAbsMag(ProbaGoldenCutFileName);
      else if(doGoldenCut)
	m2g.MaxAbsMag(doGoldenCut);
      else if (doVeryFaintCut)
        m2g.MaxAbsMag();
      
      // Simulate galaxies
      cout <<"     Simulate galaxies"<<endl;
      bool extinct=false;
        
      switch (out_type) {
      case 0:
        m2g.CreateGalCatalog(idsim, outfile, gfdz3, extinct, doVeryFaintCut, maxRadius, doGoldenCut, doProbaGoldenCut);
        cout << "Back to main program"<< endl;
        break;
      case 1:
        m2g.CreateSimpleCatalog(idsim, outfile, maxRadius);
        break;
      case 2:
        m2g.CreateTrueZFile(idsim, outfile, maxRadius);
        break;
      case 3:
        m2g.CreateNzHisto(outfile, gfdz3, extinct, maxRadius);
        break;
      default:
        throw ParmError("ERROR! out_type not understood");
        break;
      }
                
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " rdlss.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " rdlss.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " rdlss.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of rdlss.cc program  Rc= " << rc << endl;
  return rc;    
}
