/**
 * @file  addGausszerr.cc
 * @brief Add Gaussian photo-z errors to distances in a catalog
 *
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2009
 * @date 2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>


// SOPHYA
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
#include "ctimer.h"


// DirectSim
#include "geneutils.h"
#include "cat2grid.h"


void usage(void);
void usage(void) {

  cout << endl<<" Usage: addGausszerr [...options...]"           <<endl<<endl;
  
  cout << "  Add a Gaussian redshift error to the galaxies in a catalog."<<endl;
  cout << endl;
  
  cout << "  The size of the error added to the redshift is           "<<endl;
  cout << "  sigma_z*(1+z) which is converted into an equivalent      "<<endl;
  cout << "  comoving distance if it is to be added to the comoving   "<<endl;
  cout << "  distance. The value of sigma_z is supplied to the program"<<endl;
  cout << "  using option -E. The value of 'z' is either the redshift "<<endl;
  cout << "  of the galaxy in question, or the redshift value supplied"<<endl;
  cout << "  to the program using the -Z option                       "<<endl;
  cout << endl;
  
  cout << "  If the catalog has been simulated with the radial        "<<endl;
  cout << "  dimension parallel to the z-dimension use option -r. This"<<endl;
  cout << "  option has no effect if the error is being added directly"<<endl;
  cout << "  to the redshift.                                         "<<endl;
  cout << endl;
  
  
  cout << "  EXAMPLE: Add a photometric redshift error of size      "<<endl;
  cout << "  0.03*(1+z) to the spectroscopic redshifts in column      "<<endl;
  cout << "  labelled 'zs' in a file called cat.fits and output the   "<<endl;
  cout << "  augmented data to a file called Gausscat.fits            "<<endl;
  cout << endl;
  
  cout << "  $ addGausszerr -C cat.fits -O Gausscat.fits -E 0.03 -c zs"<<endl;
  cout << endl;
 
  cout << " -C : CatName     FITS filename containing catalog         "<<endl;
  cout << " -O : OutCatName  FITS file containing output catalog with "<<endl;
  cout << "                  Gaussian z errors                        "<<endl;
  cout << " -E : PZerr       Size of photometric redshift error:      "<<endl;
  cout << "                  PZerr*(1+zs) [DEFAULT=0.03]              "<<endl;
  cout << " -p : Error : Add an error on redshift from pdf (photo-z)  "<<endl;
  cout << " -S : Error : random seed, must be set for simulations     "<<endl;
  cout << " -d : delta_z the cata in 2 files :                        "<<endl;
  cout << "      in _cata if abs(zs-zp) > delta_z * (1+zs)            "<<endl;	
  cout << " -Z : zref        Add constant redshift error with value   "<<endl; 
  cout << "                  PZerr*(1+zref) [DEFAULT=no]              "<<endl;
  cout << " -c : ZSCol       Name of column of spec-z                 "<<endl;
  cout << endl;
}


int main(int narg, char *arg[]) {

  SophyaInit();
  FitsIOServerInit();
  
  cout << " ==== setting defaults ===="<<endl;
  string InCat, OutCat, OutCat2;
  double PZerr = 0.03;
  double delta_z = 100;
  double zref = -1; // becomes a value >0 if -Z option used
  string ZCol; 
  bool Zsp = false;
  string ZSCol = "zs"; // by default SPECTRO redshift column labelled "zs" is read in
  string pdfFileName = "";	//error on redfhift from pdf file, should not have any extansion (Adeline)
  bool errorFromPDF = false;
  bool RandomSeed = false;
  
  //--- decoding command line arguments 
  cout << " ==== decoding command line arguments ===="<<endl;
  char c;
  while((c = getopt(narg,arg,"hrSzC:O:E:d:p:Z:c:")) != -1) {
    switch (c) {
    case 'C' :
      InCat = optarg;
      break;
    case 'O' :
      OutCat = optarg;
      break;
    case 'E' :
      sscanf(optarg,"%lf",&PZerr);
      break;
    case 'd' :
      sscanf(optarg,"%lf",&delta_z);
      break;
    case 'p' :
      pdfFileName = optarg; 
      errorFromPDF = true;      
      break;     
    case 'S' :
      RandomSeed = true;
      cout << "Montecarlo mode"<< endl;
      break;
    case 'Z' :
      sscanf(optarg,"%lf",&zref);
      break;
    case 'c' :
      ZSCol = optarg; // z column names to read in
      //Zsp = true;
      break;
    case 'h' :
    default :
      usage(); return -1;
    }
  }
  
  //	// get two z column names
  //	if (Zsp)	
  //		{ 
  //		string delim=",";
  //		vector<string> results;
  //		stringSplit(ZCol,delim,results);
  //		vector<string>::iterator i;
  //		i = results.begin();
  //		ZOCol=*i;
  //		i++;
  //		if (i!=results.end())
  //			ZSCol=*i;
  //		}
  
  cout << "     Printing command line arguments ...             "<<endl<<endl;
  
  cout << "     Reading in observed catalogs "<< InCat                 <<endl;
  cout << "     Reading in spec redshifts from column "<< ZSCol        <<endl;
  cout << "     Catalog with added Gaussian redshift errors will be   "<<endl;
  cout << "     written to "<< OutCat                                  <<endl;
  string suff = "_cata.fits";
  OutCat2 = OutCat.substr(0,OutCat.length()-5)+suff;
  if (delta_z < 100) {
    cout << " and galaxies with abs(zP-zs) > " << delta_z <<  "x(1+zs) gathered in " << OutCat2 << endl; 
  }
  cout << "     Adding Gaussian photometric redshift error = "<<  PZerr <<"(1+zs)"<<endl;
  cout << "     Adding error to ";
  cout << "spectroscopic redshift"<<endl;
  if(errorFromPDF)
    cout << " from photo-z, from file "<< pdfFileName << endl;
  else
    cout << endl;  
  
  if(zref>0)
    cout << "     - but only equivalent to error at z = "<< zref <<endl;
  cout <<endl;
  
  
  try {
    
    Timer tm("AddGaussZerr");
    
    RandomGenerator rg; // random number generator
    if (RandomSeed) {
      rg.AutoInit(0);
      cout << "Seed automatically generated" << endl;
    } else {
      long seed=1;
      rg.SetSeed(seed);
    }
    
    // Read input catalog
    FitsInOutFile fin(InCat,FitsInOutFile::Fits_RO);
    fin.MoveAbsToHDU(2);
    SwFitsDataTable dt(fin,512,false);
    sa_size_t ng = dt.NEntry(); // number of galaxies
    sa_size_t nc = dt.NCols();  // number of columns
    DataTableRow row = dt.EmptyRow();
    cout << "     Number of galaxies in catalog = "<< ng ;
    cout << ", number of columns in catalog = "<< nc <<endl;
    
	     		
    // Set cosmology (read this from galaxy catalog header)
    cout << "     Initialise cosmology: (read from catalog)"<<endl;
    
    //////////  modif Adeline : read cosmo in Fits_RO header
    string H0_s, OmegaM_s, OmegaL_s, OmegaB_s, OmegaR_s, wDE_s, wDA_s, Sigma8_s, Ns_s, cell_s;
    double h, OmegaM, OmegaL, OmegaB, OmegaR, wDE, wDA, Sigma8, n_s,cell;
    H0_s = fin.KeyValue("H0");
    OmegaM_s = fin.KeyValue("OMEGAM0");
    OmegaL_s = fin.KeyValue("OMEGADE0");
    OmegaB_s = fin.KeyValue("OMEGAB0");
    OmegaR_s = fin.KeyValue("OMEGAR0");
    wDE_s = fin.KeyValue("DE_W0");
    wDA_s = fin.KeyValue("DE_WA");
    Sigma8_s = fin.KeyValue("SIGMA8");
    Ns_s = fin.KeyValue("N_S");
    cell_s = fin.KeyValue("CELL");

    h = atof(H0_s.c_str()) / 100.;
    OmegaM = atof(OmegaM_s.c_str());
    OmegaL = atof(OmegaL_s.c_str());
    OmegaB = atof(OmegaB_s.c_str());
    OmegaR = atof(OmegaR_s.c_str());
    wDE = atof(wDE_s.c_str());
    wDA = atof(wDA_s.c_str());
    Sigma8 = atof(Sigma8_s.c_str());
    n_s = atof(Ns_s.c_str());
    cell = atof(cell_s.c_str());

    SimpleUniverse su(h, OmegaM, OmegaL);
    su.SetOmegaBaryon(OmegaB);
    su.SetOmegaRadiation(OmegaR);
    su.SetSigma8(Sigma8);
    su.SetSpectralIndex(n_s);
    su.SetFlatUniverse_OmegaLambda(); // Cecile modif - no be sure that it is flat by adjusting OmegaLambda
    cout << "     OmegaK="<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
    cout << ", OmegaL="<< su.OmegaLambda() <<", OmegaB="<< su.OmegaBaryon();
    cout << ", Omega_rad=" << su.OmegaRadiation() << ", Omega_cdm=" << su.OmegaCDM() <<", H0="<< su.H0() << endl;
    cout << "check flatness: OmegaTot=" << su.OmegaTotal() << endl;
    if (wDE != -1 or wDA !=0)  
      su.SetDarkEnergy(su.OmegaLambda(),wDE,wDA);
    
    cout << " and w0=" << su.wDE() << ", wA=" << su.waDE() << ", sigma8=" << su.Sigma8() << endl;
    cout << "Spectral index=" << su.Ns() << endl;
    cout << endl;
	      		
    // Get fixed error in comoving z-dimension (if using)
    double ErrFix = 0;
    
    // Create z-dcomoving look up table incase needed
    SInterp1D z2dist, dist2z;
    int_8 nz=1000;
    vector<double> zrs, codist;
    double minz=0, maxz=10;
    double dz = (maxz-minz)/(nz-1);
    for(int kk=0; kk<nz; kk++) {
      
      double zs=minz+kk*dz;
      su.SetEmissionRedShift(zs);
      double cod =su.RadialCoordinateMpc(); // radial distance 
      zrs.push_back(zs);
      codist.push_back(cod); 
    }
    z2dist.DefinePoints(zrs,codist,minz,maxz,2*nz);
    dist2z.DefinePoints(codist,zrs,codist[0],codist[codist.size()-1],2*nz);
    
    sa_size_t nobj = 20; // number of objects to print to the screen
    sa_size_t nskip = (sa_size_t)ng/nobj;
    long j=0;    
    
    // Create new catalog 
    cout << "     Creating new catalog called "<< OutCat <<endl;
    FitsInOutFile swf(OutCat, FitsInOutFile::Fits_Create);	
    SwFitsDataTable gals(swf, 2048);
    // swf2 is created but remains empty if delta_z < 100
    FitsInOutFile swf2(OutCat2, FitsInOutFile::Fits_Create);	
    SwFitsDataTable gals2(swf2, 2048);

    for (int i=0;i<nc;i++) { // add same columns input file has
      
      cout <<"  Col number = "<< i <<", Col name = "<< dt.NomIndex(i);
      cout <<", Col type = "<< dt.GetColumnType(i) <<endl;
      if (dt.GetColumnType(i)<2)
	gals.AddIntegerColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>1) && (dt.GetColumnType(i)<3) )
	gals.AddLongColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>2) && (dt.GetColumnType(i)<4) )
	gals.AddFloatColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>3) && (dt.GetColumnType(i)<5) )
	gals.AddDoubleColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>4) && (dt.GetColumnType(i)<6) )
	gals.AddComplexColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>5) && (dt.GetColumnType(i)<7) )
	gals.AddDoubleComplexColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>6) && (dt.GetColumnType(i)<8) )
	gals.AddStringColumn(dt.NomIndex(i));
    }
    if (delta_z < 100) 
      for (int i=0;i<nc;i++) { // add same columns input file has
	if (dt.GetColumnType(i)<2)
	  gals2.AddIntegerColumn(dt.NomIndex(i));
	else if ( (dt.GetColumnType(i)>1) && (dt.GetColumnType(i)<3) )
	  gals2.AddLongColumn(dt.NomIndex(i));
	else if ( (dt.GetColumnType(i)>2) && (dt.GetColumnType(i)<4) )
	  gals2.AddFloatColumn(dt.NomIndex(i));
	else if ( (dt.GetColumnType(i)>3) && (dt.GetColumnType(i)<5) )
	  gals2.AddDoubleColumn(dt.NomIndex(i));
	else if ( (dt.GetColumnType(i)>4) && (dt.GetColumnType(i)<6) )
	  gals2.AddComplexColumn(dt.NomIndex(i));
	else if ( (dt.GetColumnType(i)>5) && (dt.GetColumnType(i)<7) )
	  gals2.AddDoubleComplexColumn(dt.NomIndex(i));
	else if ( (dt.GetColumnType(i)>6) && (dt.GetColumnType(i)<8) )
	  gals2.AddStringColumn(dt.NomIndex(i));
	
      }

    
    // then add one more that will contain z's w/ Gaussian errors
    if (errorFromPDF) 
      gals.AddFloatColumn("ZP"); else gals.AddFloatColumn("ZG"); 
    if (delta_z < 100) 
          if (errorFromPDF) 
	    gals2.AddFloatColumn("ZP"); else gals2.AddFloatColumn("ZG"); 


    DataTableRow rowout = gals.EmptyRow();    
    
    
    // Find required col names (only required if adding to z-dimension)
    sa_size_t Izs = dt.IndexNom(ZSCol);
    sa_size_t Ith,Iph;
    Iph = dt.IndexNom("PHI");
    Ith = dt.IndexNom("THETA");
    
    // loop over catalog
    cout << "     Starting loop over catalog ... "<<endl;		
    
    Cat2Grid cat(dt,su,rg,dist2z,z2dist,ZSCol,ZSCol);
    GalRecord grec;    
    //compute error from photo-z (Adeline)
    string zcat_min_s,zcat_max_s;
    double zcat_min, zcat_max;
    zcat_min_s = fin.KeyValue("ZCAT_MIN");
    zcat_max_s = fin.KeyValue("ZCAT_MAX");
    zcat_min = atof(zcat_min_s.c_str());
    zcat_max = atof(zcat_max_s.c_str());
    
    if(errorFromPDF)    
      cat.SetPhotozErrRedshift(pdfFileName,RandomSeed,zcat_min,zcat_max); // compute error from photo-z (Adeline)

	  
    size_t totcellcnt = ng;
    size_t cellcnt = 0;
    ProgressBar pgb(totcellcnt, ProgBarM_Time);
    
    cout << "*** Start loop over "<< ng << " objects"<<endl;    
    clock_t t1, t2;
    double dif=0;
    for (long i=0;i<ng;i++) {        
      dt.GetRow(i,row);
      double zs = row[Izs]; // spec-z of galaxy 
      
      t1 = clock();
      if (i == 1) {
	   cout <<" time per object = " << dif << endl;
	   cout <<" estimated time " << dif * (ng) << endl;
      }      
      
      double zG, th, ph;	
      // Adding straight to spec-z
	
      if(!errorFromPDF){
	double sigz = (zref>0) ? PZerr*(1+zref) : PZerr*(1+zs); // error at z specified OR error at z of galaxy
	zG = zs + sigz*rg.Gaussian();
	if(i<10)  cout <<"zs "<< zs << " zG="<< zG <<endl;
      }
      else {
	cat.Row2Record(row,grec);	
	zG = cat.ComputeErrorFromPhotoz(grec);
	zs = grec.zo;
	if(i<10) {
	  cout <<"    galid="<<i<<": ";
	  cout <<", zs = "<< zs  ;
	  cout <<", type="<<grec.type<<", MA="<<grec.MB <<", zP="<< zG <<endl;
	}
      
      
	// don't think I need to convert the th and ph					
      }
      
      // make sure there are no unphysical redshifts
      if (zG<0)
	zG=0.01;
      
      // fill output row with same values as input row
      for (int n=0; n<row.Size();n++)
	rowout[n]=row[n];
      // add 'redshift' with error to last column. 
      rowout[row.Size()] = zG;
      
      if (abs(zs - zG) > delta_z*(1.+zs)) {
	//	cout << "CATA "<< zs << " " <<  zG<< " "  << delta_z*(1.+zs) << " " << abs(zs - zG) <<endl;
	gals2.AddRow(rowout);
      }
      else gals.AddRow(rowout);
      

      cellcnt ++;
      pgb.update(cellcnt);
      
      t2 = clock();
      dif = (t2 - t1) / (double) CLOCKS_PER_SEC;
      
    }// end loop over file
    cout << "nCell = "<< cellcnt << "  "<< row.Size() << "  "<<rowout.Size() << endl;
    
    //modified by Adeline : write cosmo parameters in file header
    swf.WriteKey("H0", su.H0()," Cosmo.Param H0");
    swf.WriteKey("OMEGAM0", su.OmegaMatter()," Cosmo.Param OmegaMatter0 ");
    swf.WriteKey("OMEGAB0", su.OmegaBaryon()," Cosmo.Param OmegaBaryon0");
    swf.WriteKey("OMEGAR0", su.OmegaRadiation()," Cosmo.Param OmegaRadiation0");
    swf.WriteKey("OMEGAT0", su.OmegaTotal()," Cosmo.Param OmegaTot0");
    swf.WriteKey("OMEGADE0", su.OmegaLambda(),"  Cosmo.Param OmegaLambda0 (dark energy density)");
    swf.WriteKey("OMEGADK", su.OmegaCurv(),"  Cosmo.Param OmegaK ");
    swf.WriteKey("DE_W0", su.wDE(), " Cosmo.Param w0 (dark energy eq.state)");
    swf.WriteKey("DE_WA",su.waDE() , " Cosmo.Param wA (dark energy eq.state)"); 
    swf.WriteKey("SIGMA8", su.Sigma8(), " Cosmo.Param sigma8_0");
    swf.WriteKey("N_S",su.Ns()," Cosmo.Param n_s (spectral index scalar fluct.)");
    // additionnal keyword
    swf.WriteKey("CELL",cell,"cell size of input grid [Mpc]");
    
    if (delta_z < 100) {
      swf2.WriteKey("H0", su.H0()," Cosmo.Param H0");
      swf2.WriteKey("OMEGAM0", su.OmegaMatter()," Cosmo.Param OmegaMatter0 ");
      swf2.WriteKey("OMEGAB0", su.OmegaBaryon()," Cosmo.Param OmegaBaryon0");
      swf2.WriteKey("OMEGAR0", su.OmegaRadiation()," Cosmo.Param OmegaRadiation0");
      swf2.WriteKey("OMEGAT0", su.OmegaTotal()," Cosmo.Param OmegaTot0");
      swf2.WriteKey("OMEGADE0", su.OmegaLambda(),"  Cosmo.Param OmegaLambda0 (dark energy density)");
      swf2.WriteKey("OMEGADK", su.OmegaCurv(),"  Cosmo.Param OmegaK ");
      swf2.WriteKey("DE_W0", su.wDE(), " Cosmo.Param w0 (dark energy eq.state)");
      swf2.WriteKey("DE_WA",su.waDE() , " Cosmo.Param wA (dark energy eq.state)"); 
      swf2.WriteKey("SIGMA8", su.Sigma8(), " Cosmo.Param sigma8_0");
      swf2.WriteKey("N_S",su.Ns()," Cosmo.Param n_s (spectral index scalar fluct.)");
      // additionnal keyword
      swf2.WriteKey("CELL",cell,"cell size of input grid [Mpc]");
    }
    
    cout << "Check cosmo parameters : " << endl;
    cout << "  OmegaK="<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
    cout << ", OmegaL="<< su.OmegaLambda() <<", OmegaB="<< su.OmegaBaryon()  ;
    cout << ", Omega_rad=" << su.OmegaRadiation() << ", H0=" << su.H0() << ", Sig8=" << su.Sigma8() <<", n_s=" << su.Ns()<<endl; 
    cout << ", Omega_curv=" << su.OmegaCurv() << ", DE_W0=" << su.wDE() << ", DE_WA=" << su.waDE() <<endl; 
    cout << endl;
    // end modifications
    
  }
  catch(PThrowable exc ) {
    cerr << "addGausszerr.cc , Catched exception: \n" << exc.what() << endl;
  }
  catch(std::exception ex) {
    cerr << "addGausszerr.cc , Catched exception ! " << (string)(ex.what()) << endl;
  }
  catch(...) {
    cerr << "addGausszerr.cc , Catched ... ! " << endl;
  }
  
  cout << "--------------- addGausszerr.cc / END --------------------- " << endl;
}
