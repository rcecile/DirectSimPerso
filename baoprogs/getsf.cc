/**
 * @file  getsf.cc
 * @brief Compute selection functions \f$ n_{pz}(z)/n_{true}(z) \f$ and
 *        \f$ n_{sz}(z)/n_{true}(z) \f$
 *
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: Aug 2010
 * @date Aug 2010
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

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

// DirectSim
#include "geneutils.h"
#include "cat2grid.h"
#include "cosmocalcs.h"
#include "constcosmo.h"

void usage(void);
void usage(void) {

	cout << endl<<" Usage: getsf [...options...]" << endl<<endl;
	
    cout << "  Compute selection functions n_{pz}(z)/n_{true}(z) and    "<<endl;
    cout << "  n_{sz}(z)/n_{true}(z)                                    "<<endl;
    cout << endl;
	
	cout << "  Read in a catalog containing all the simulated galaxies  "<<endl;
	cout << "  (or the histogram already made from it ...)  "<<endl;
	cout << "  and a catalog of the observed galaxies. The catalog that "<<endl;
	cout << "  that contains all the simulated galaxies may be stored in"<<endl;
	cout << "  more than one file. Supply all these filenames separated "<<endl;
	cout << "  by commas to the program with the -F option. Supply the  "<<endl;
    cout << "  filename of the catalog of observed galaxies with the -O "<<endl;
    cout << "  option.                                                  "<<endl;
	cout << endl;
	
	cout << "  The -z option passes the names of the columns to read the"<<endl;
	cout << "  redshifts from, separated by commas, in the following    "<<endl;
	cout << "  order: observed photometric redshifts, observed spec     "<<endl;
	cout << "  redshifts, simulated redshifts. The default names are:   "<<endl;
	cout << "  'zp', 'z' and 'z'                                        "<<endl;
	cout << endl;

    cout << "  The selection function n_{sz}/n_{true} will be written to"<<endl;
    cout << "  the file [SFTextFile]_specz_nofz.txt and the selection   "<<endl;
    cout << "  function n_{pz}/n_{true} will be written to the file     "<<endl;
    cout << "  [SFTextFile]_nofz.txt, where SFTextFile is supplied to   "<<endl;
    cout << "  the program using the -o option.                         "<<endl;
	
	cout << "  If the -d option is used ppf files containing the n(z)   "<<endl;
	cout << "  are output to a file                                     "<<endl;
	cout << endl;
	
	cout << "  EXAMPLES: " <<endl;
	
	cout << "  $ getsf -F full.fits -O obs.fits -o selectfunc -z zp,z,z "<<endl;
    cout << "  $ getsf -F full1.fits,full2.fits,full3.fits -O obs.fits  "<<endl;
    cout << "                                  -o selectfunc -z zp,z,z  "<<endl;
    cout << endl;
	
	cout << " -F : FullCat    FITS filename containing simulated catalog(s)  "<<endl;
	cout << " -H : FullHist   text file containing the histogram from simulated catalog(s) "<<endl;
	cout << " -O : ObsCat     FITS file containing observed catalog(s)       "<<endl;
	cout << " -o : sfunc      root name of selection function file           "<<endl;
	cout << " -z : zp,zs,zz   column names to read redshifts from (see above)"<<endl;
	cout << " -d : [noarg]    Save n(z) Histos to a ppf file                 "<<endl;
	//cout << " -N : nFiles: Number of FITS files containing RDLSS output [default=1]"<<endl;
	cout << endl;
}


int main(int narg, char *arg[])	{

	SophyaInit();
	FitsIOServerInit();
  
	string FullCat, FullHist, ObsCat, SFFileName;
	//int nFiles = 1;
	bool DoDebug = false;
	bool MakeFullHist = true;
	bool isZColSpecified = false;
	string ZCol;
	string ZOCol = "zp"; // read in photo-z's as OBSERVED redshifts from OBSERVED catalog
	string ZSCol = "z";  // read in SPECTRO-z from OBSERVED catalog
	string ZFCol = "z";  // read in SPECTRO-z from FULL catalog
	
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hdRF:H:O:o:z:")) != -1) { 
	    switch (c) {
		    case 'F' :
			    FullCat = optarg;
			    break;
		    case 'H' :
			    FullHist = optarg;
			    MakeFullHist = false;
			    break;
///			case 'N' :
//			    sscanf(optarg,"%d",&nFiles);
//			    break;
		    case 'O' :
			    ObsCat = optarg;
			    break;
		    case 'o' :
			    SFFileName = optarg;
			    break;
		    case 'z' :
			    ZCol = optarg;
			    isZColSpecified = true;
			    break;
 	            case 'd' :
			    DoDebug=true;
			    break;
		    case 'h' :
		        default :
			    usage(); return -1;
		        }
	        }

    // read in names of the redshift columns
	if (isZColSpecified) { 
		string delim=",";
		vector<string> results;
		stringSplit(ZCol,delim,results);
		vector<string>::iterator i;
		i = results.begin();
		ZOCol=*i; // column name of observed redshifts
		i++;
		if (i!=results.end())
			ZSCol=*i; // column name of spec redshifts
		i++;
		if (i!=results.end())
			ZFCol=*i; // column name of spec redshifts in full catalog
		}
	
	cout << "     Printing command line arguments ... "<<endl<<endl;
	cout << "     Observed catalog in file(s) "<< ObsCat <<endl;
	cout << "     OBSERVED redshifts to be read from "<< ZOCol <<" column"<<endl;
	cout << "     SPECTRO redshifts to be read from "<< ZSCol <<" column"<<endl;
	if (MakeFullHist)
	  cout << "     True redshifts in file(s) " << FullCat << endl;
	else 
	  cout << "     True redshifts histogram in file " << FullHist << endl;

	cout << "     SPECTRO redshifts to be read from "<< ZFCol <<endl;
	cout << "     Selection function will be written to "<< SFFileName <<"_nofz.txt and ";
	cout << SFFileName <<"_specz_nofz.txt"<<endl;

	if (DoDebug)
		cout << "     Saving n(z)'s to ppf file "<< SFFileName <<"_histo.ppf"<< endl;
	cout <<endl;
	
  try {
  
		// Read in redshifts from observed catalog
		cout <<"0/ Read in observed catalog "<< ObsCat <<endl;

		// Possibility to have a list of observed catalogs (Cecile)
		string delim=",";
		vector<string> fobsnames;
		stringSplit(ObsCat,delim,fobsnames);
  
		//use the first to define the dt class
		FitsInOutFile fin(fobsnames[0], FitsInOutFile::Fits_RO);
		fin.MoveAbsToHDU(2);
		SwFitsDataTable dt(fin,512,false);
	
		// Set cosmology
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
			
		
		// Calculate selection function
		RandomGenerator rg;
		string tmp="tmptmp";
		FitsInOutFile fos(tmp,FitsInOutFile::Fits_Create);
		bool RadialZ=false;
		Cat2Grid cat(dt, su, rg, fos, ZOCol, ZSCol,RadialZ, 0., true, ObsCat);

		string OutRoot = "tmp";
		if (DoDebug)
			cat.SetDebugOutroot(OutRoot);
		
 	    cout <<" DONE" << endl;
	    if (MakeFullHist)
	      cat.SaveSelecFunc(SFFileName, FullCat, ObsCat, ZFCol, ZSCol, ZOCol,MakeFullHist);
	    else
	      cat.SaveSelecFunc(SFFileName, FullHist, ObsCat, ZFCol, ZSCol, ZOCol,MakeFullHist);

	    // both SPECTRO-z and OBSERVED-z sf's are computed here
		if( remove(tmp.c_str()) != 0 )
		  cout << "Error deleting temporary file" << endl;		       
    } 
    catch(PThrowable exc ) {
        cerr << "getsf.cc , Catched exception: \n" << exc.what() << endl;
        }
    catch(std::exception ex) {
        cerr << "getsf.cc , Catched exception ! " << (string)(ex.what()) << endl;
        }
    catch(...) {
        cerr << "getsf.cc , Catched ... ! " << endl;
        }

    cout << "--------------- getsf.cc / END --------------------- " << endl;
}
