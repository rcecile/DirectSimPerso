#include "cat2grid.h"
#include "ctimer.h"
#include "progbar.h"

using namespace std;

// autre input ng tot pour avoir ng(z)
// faire grille random comme gal : tirage des "gal-bruit" dans cellule + erreur photoZ

//Improved constructor:
Cat2Grid::Cat2Grid(SwFitsDataTable& dt, SimpleUniverse& su, RandomGenerator& rg,
                   FitsInOutFile& fos, string ZOCol, string ZSCol,
                   bool RadialZ, double PZerr, bool Print, string ObsCat, double ratio_cell) 
  : dt_(dt) , su_(su) , rg_(rg) , fos_(fos) , defselfunc_(1.), selfuncp_(&defselfunc_), defbiasfunc_(0.), biasp_(&defbiasfunc_)
  , defsigrfunc_(0.), sigrp_(&defsigrfunc_), ratio_cell_(ratio_cell)
  , ZOCol_(ZOCol) , ZSCol_(ZSCol) , RadialZ_(RadialZ) , PZerr_(PZerr), ObsCat_(ObsCat)
                  // Constructor for this class reads in the phi, theta, redshift of all galaxies in the catalog
                  // from data table dt
                  // su:          class which holds the cosmological parameters and can be used to 
                  //        calculate cosmological quantities
                  // rg:          class used to generate random numbers
                  // ZOCol: col name of OBSERVED redshifts to read in, could be spec-z, gauss-z, phot-z
                  // ZSCol: col name of SPECTRO redshifts to be read in (must be spec-z)
                  // RadialZ: if true z-dimension IS radial direction
                  // PZerr: size of photometric redshift errors: sigz = PZerr(1+z); 
                  // Print: if true prints some extra info to screen whilst in constructor
{
  if (Print)
    cout <<endl<<"    In Cat2Grid constructor ..."<<endl;

  DVList dvl;
  dvl = dt.Info();
  mean_overdensity_ = dvl["MeanOverDensity"];
  // Cecile, to deal with sim grid with a different resolution
  mean_overdensity_ /= ratio_cell;
  cout << "    Mean over-density of fudged SimLSS grid = "<< mean_overdensity_ <<endl;
        
  DoDebug_ = false; // flag to do debug
  debugoutroot_ = "tmp"; // if want to debug, files are written to filenames with this root
  sfcompute_ = false; // flag becomes true after selection function is set
  doBiasCorr_ = false; // flag becomes true after selection function is set
  PZDerr_=0;// set photoz error in Mpc equal to zero to start with
  AddGaussErr_     = false; // flag to add Gaussian photo-z err
  AddGaussErrReds_ = false; // flag to add Gaussian photo-z err on redshift
  AddPhotoErrReds_ = false; // flag to add photo-z photo-z err on redshift

  // Set counting galaxies to ZERO
  wnrand_=0; // weighted number of "galaxies" in random grid

// num "gals" INSIDE WEIGHTED grid
  // Data table column names
  Ic1_ = dt_.IndexNom("PHI");
  Ic2_ = dt_.IndexNom("THETA");
  Izs_ = dt_.IndexNom(ZSCol_); // the SPEC-Z column
        
  if (Print) {
    cout << "    Reading SPECTRO redshifts from column name "<< ZSCol_ <<endl;
    cout << "    Reading OBSERVED redshifts from column name "<< ZOCol_ <<endl;
  }
  Iz_ = dt_.IndexNom(ZOCol_);
  Iid_= dt_.IndexNom("GALID");
        
  Itype_= dt_.IndexNom("TYPE"); 
  Imag_ = dt_.IndexNom("MB"); 

  ngall_= 0;
  string delim=",";
  vector<string> fobsnames;
  stringSplit(ObsCat_,delim,fobsnames);
  int Nobsfiles = fobsnames.size();
  for (int ic=0; ic<Nobsfiles; ic++) {
    string ObsCatFile = fobsnames[ic];
    FitsInOutFile fin(ObsCatFile,FitsInOutFile::Fits_RO);
    fin.MoveAbsToHDU(2);
    SwFitsDataTable dtobs(fin,512,false);
    ngall_ +=  dtobs.NEntry();
    if (Print)    cout <<" ---------- number of galaxies in catalog "<< fobsnames[ic] << " = " << dtobs.NEntry() <<endl;
  }

  if (Print) {
    cout <<endl <<"    TOTAL number of galaxies in catalog ="<< ngall_ <<endl;
    cout <<"    Number of columns in catalog = "<< dt_.NCols() <<endl;
    cout <<"    Angle phi = column "<< Ic1_ <<", angle theta = column "<< Ic2_ <<endl;
    cout <<"    Spec-z = column "<< Izs_ <<", z used in analysis = column "<< Iz_ <<endl;
    cout <<"    (the above will be the same column when reading in spec-z)"<<endl;
  }
                
  // Create a redshift - distance look up table to save time not using su_ class computation
  int_8 nz=1000000;
  vector<double> zrs, codist;
  double minz=0, maxz=10;
  double dz = (maxz-minz)/(nz-1);
  for(int kk=0; kk<nz; kk++) {
                
    double zs=minz+kk*dz;
    su_.SetEmissionRedShift(zs);
    double cod =su_.RadialCoordinateMpc(); // radial distance 
    zrs.push_back(zs);
    codist.push_back(cod); 
  }
  z2dist_.DefinePoints(zrs,codist,minz,maxz,2*nz);
  dist2z_.DefinePoints(codist,zrs,codist[0],codist[codist.size()-1],2*nz);

  if (Print)
    cout <<"    Exit Cat2Grid constructor ..."<<endl<<endl;
};


// make interpolation table outside Cat2Grid, just use Row2Record, Rec2EuclidCoord functions
Cat2Grid::Cat2Grid(SwFitsDataTable& dt,SimpleUniverse& su,RandomGenerator& rg,
                   SInterp1D dist2z, SInterp1D z2dist,string ZOCol,string ZSCol,bool RadialZ) 
  : dt_(dt) , su_(su) , rg_(rg) , defselfunc_(1.), selfuncp_(&defselfunc_) , fos_(fosdefault_) ,
    ZOCol_(ZOCol) , ZSCol_(ZSCol) , dist2z_(dist2z) , z2dist_(z2dist) , RadialZ_(RadialZ)
{

  bool Print=false;
  if (Print)
    cout <<endl<<"    In Cat2Grid constructor ..."<<endl;
                
  PZerr_=0;
        
  DoDebug_ = false; // flag to do debug
  debugoutroot_ = "tmp"; // if want to debug, files are written to filenames with this root
  sfcompute_ = false; // flag becomes true after selection function is set
  PZDerr_=0;// set photoz error in Mpc equal to zero to start with
  AddGaussErr_     = false; // flag to add Gaussian photo-z err
  AddGaussErrReds_ = false; // flag to add Gaussian photo-z err on redshift
        
  // Data table column names
  Ic1_ = dt_.IndexNom("PHI");
  Ic2_ = dt_.IndexNom("THETA");
  
  Izs_ = dt_.IndexNom(ZSCol_); // the SPEC-Z column
  Iz_ = dt_.IndexNom(ZOCol_);
  Itype_= dt_.IndexNom("TYPE"); 
  Imag_ = dt_.IndexNom("MB"); 

  if (Print)
    cout << "    Reading SPECTRO redshifts from column name "<< ZSCol_ <<endl;
  if (Print)
    cout << "    Reading OBSERVED redshifts from column name "<< ZOCol_ <<endl;

        
  ngall_=dt_.NEntry();
  if (Print) {
    cout <<"    TOTAL number of galaxies in catalog ="<< ngall_ <<endl;
    cout <<"    Number of columns in catalog = "<< dt_.NCols() <<endl;
    cout <<"    Angle phi = column "<< Ic1_ <<", angle theta = column "<< Ic2_ <<endl;
    cout <<"    Spec-z = column "<< Izs_ <<", z used in analysis = column "<< Iz_ <<endl;
    cout <<"    (the above will be the same column when reading in spec-z)"<<endl;
  }
                
  if (Print)
    cout <<"    Exit Cat2Grid constructor ..."<<endl<<endl;
};


// copy constructor
Cat2Grid::Cat2Grid(Cat2Grid const& a)
  : dt_(a.dt_) , su_(a.su_) , rg_(a.rg_), fos_(a.fos_) ,
    defselfunc_(a.defselfunc_.sfv_), selfuncp_(a.selfuncp_)
{
  cout <<"    Cat2Grid COPY constructor"<<endl;
  Set(a);
};


// Copy Cat2Grid 
Cat2Grid& Cat2Grid::Set(const Cat2Grid& a)
{
  //Initialised in constructor
                                        
  PZerr_=a.PZerr_;
  debugoutroot_=a.debugoutroot_;        
  PZDerr_=a.PZDerr_;    
  AddGaussErr_    =a.AddGaussErr_;
  AddGaussErrReds_=a.AddGaussErrReds_;
  doBiasCorr_=a.doBiasCorr_;
  selfuncp_=a.selfuncp_;
  ZSCol_=a.ZSCol_;
  ZOCol_=a.ZOCol_;

  zsmin_=a.zsmin_;
  zsmax_=a.zsmax_;      
  xmin_=a.xmin_;
  xmax_=a.xmax_;
  ymin_=a.ymin_;
  ymax_=a.ymax_;
  zmin_=a.zmin_;
  zmax_=a.zmax_;
                 
  ngo_=a.ngo_;            
  ngall_=a.ngall_;            
  ng_=a.ng_;            
  ngout_=a.ngout_;            
  ngw_=a.ngw_;            
  wnrand_=a.wnrand_;            

  //Initialised in ComputeSFfromLF
  /*phistar_=a.phistar_;
    Mstar_=a.Mstar_;
    alpha_=a.alpha_;
    mlim_=a.mlim_;
    Mc_=a.Mc_;*/
                        
  //Initialised in ConstructGrid
  Vol_=a.Vol_;  
  Nx_=a.Nx_;
  Ny_=a.Ny_;
  Nz_=a.Nz_;
  volgrid_=a.volgrid_;
  cellsize_=a.cellsize_;
        
  Npix_=a.Npix_;
  if (a.wrgals_.Size()>0)
    {
      wrgals_=a.wrgals_;
    }
  zc_ = a.zc_;
  alph_=a.alph_;
         
};


void Cat2Grid::SetGrid(int_8 Nx, int_8 Ny, int_8 Nz, double R, double zref)
// Computes range of grid to lay over galaxy simulation
// from SimLSS cube parameters
{
  cout <<endl<<"    Cat2Grid::SetGrid()"<<endl;
        
  cout << "    Using SimLSS grid to define grid"<<endl;
  cout << "    SimLSS grid: "<< Nx <<","<< Ny <<","<< Nz <<" pixels"<<endl;
  cout << "    Centered at redshift z = "<< zref <<endl;
  cout << "    Resolution R = "<< R <<endl;
        
  cellsize_ = R;
  zref_ = zref;
        
  su_.SetEmissionRedShift(zref_);
  DCref_ = su_.RadialCoordinateMpc();
  double idzd = ceil((double)Nz/2); // index of center pixel in z direction
  double idxd = ceil((double)Nx/2);
  double idyd = ceil((double)Ny/2);     
  idx_ = idxd, idy_ = idyd,idz_ = idzd;
  cout <<"    Indices of center pixels: Ix="<< idx_ <<", Iy="<< idy_ <<", Iz="<< idz_ <<endl;
        
  Xmin_ = -(idx_-1)*cellsize_-cellsize_/2;
  Xmax_ = (idx_-1)*cellsize_+cellsize_/2;
  Ymin_ = -(idy_-1)*cellsize_-cellsize_/2;
  Ymax_ = (idy_-1)*cellsize_+cellsize_/2;
  Zmin_ = DCref_-(idz_-1)*cellsize_-cellsize_/2;
  Zmax_ = DCref_+(idz_-1)*cellsize_+cellsize_/2;
        
  // cell coords were given by:
  //X = (i-(idmidx_-1))*Dx_; }
  //(j-(idmidy_-1))*Dy_; }
  //(k-(idmidz_-1))*Dz_+DCref_; }
        
  cout <<"    Initial survey boundaries: "<< Xmin_ <<"<X<"<< Xmax_;
  cout <<", "<< Ymin_ <<"<Y<"<< Ymax_ <<", "<< Zmin_ <<"<Z<"<< Zmax_ <<endl;

  //    // round down to nearest integer SHOULD REMOVE THIS
  //    Xmin_=floor(Xmin_);
  //    Xmax_=floor(Xmax_);
  //    Ymin_=floor(Ymin_);
  //    Ymax_=floor(Ymax_);
  //    Zmin_=floor(Zmin_);
  //    Zmax_=floor(Zmax_);
  //    cout <<"    Rounded grid boundaries: "<<Xmin_<<"<X<"<<Xmax_<<", "<<Ymin_<<"<Y<"<<Ymax_<<", "<<Zmin_<<"<Z<"<<Zmax_<<endl;
        
  // calculate EXACT number of cells that reach from one side of sim to the other (including remainder)
        
  double Lx, Ly, Lz; 
  Lx = (Xmax_-Xmin_)/cellsize_;
  Ly = (Ymax_-Ymin_)/cellsize_;
  Lz = (Zmax_-Zmin_)/cellsize_;
  cout <<"    Number of cells which cover survey in each direction"<<endl;
  cout <<"    Nx="<< Lx <<", Ny="<< Ly <<", Nz="<< Lz <<endl;
  // round up to get INTEGER number of cells that cover each dimension
  // integer number cells needed
  Nx_ = (sa_size_t)ceil(Lx); Ny_ = (sa_size_t)ceil(Ly); Nz_ = (sa_size_t)ceil(Lz); 
  cout <<"    Integer number of cells: Nx="<<Nx_<<", Ny="<<Ny_<<", Nz="<<Nz_<<endl;
        
  if (Nx_!=Nx) {
    cout << "TEMP FUDGE, setting Nx to input value"<<endl;
    Nx_=Nx;
  }
  if (Ny_!=Ny) {
    cout << "TEMP FUDGE, setting Ny to input value"<<endl;
    Ny_=Ny;
  }
  if (Nz_!=Nz) {
    cout << "TEMP FUDGE, setting Nz to input value"<<endl;
    Nz_=Nz;
  }
       
  su_.SetEmissionRedShift(zref_);
  double cod = su_.RadialCoordinateMpc(); // radial distance
  cout <<"    Redshift of central pixel = "<< zref_;
  cout <<", comoving distance @zref = "<< cod <<endl;
        
  // INITIALISE ARRAYS TO GRID SIZE
  int ndim=3;
  sa_size_t mydim[ndim];
  mydim[0]=Nx_; mydim[1]=Ny_; mydim[2]=Nz_;
  if (sfcompute_) {
    wrgals_.SetSize(ndim, mydim);       // weighted galaxy density field: random catalog
  }
                
  cout <<"    EXIT Cat2Grid::SetGrid()"<<endl<<endl<<endl;
};


void Cat2Grid::SaveSelecFunc(string SFTextFile, string FullCat, string ObsCat, string ZFCol, string  ZSCol,  string ZOCol, bool MakeFullHist)
// Create Histo of observed redshifts/true redshifts and save to
// a text file called [SFTextFile]_nofz.txt
// Reads in a Fits file containing ALL the true redshifts in the
// simulation, loops over these redshifts and adds them to Histo
// object
// Loops over observed catalog and adds observed redshifts to
// Histo object
{
  cout << "    Cat2Grid::SaveSelecFunc()"<<endl;
  cout << "OBSCAT : " << ObsCat_ << endl;
  Timer tm("SaveSelecFunc");
  
  // Get full catalog file name(s)
  string delim=",";
  vector<string> fnames;
  stringSplit(FullCat,delim,fnames);
  int Nfiles = fnames.size();
  
  /*// Set name of catalog containing the "full" set of z's
    string FullCatFile;
    if(Nfiles>1)
    cout <<"    Reading in "<< Nfiles <<" full catalog files"<<endl; 
    else
    FullCatFile = FullCat;*/
  
  // Initialise random generator to specified seed
  long seed=1; 
  rg_.SetSeed(seed);
  
  // Set up Histogram bins
  double minz = 0, maxz=10;
  int nbin=5000;
  Histo nzFC(minz, maxz, nbin); // full catalog Histo
  Histo nzOC(minz, maxz, nbin); // obs catalog Histo
  // (will be same as above if spec-z used in analysis)
  Histo nzOCsz(minz, maxz, nbin); // obs catalog Histo (using spec-z) 
  
  double minzt, maxzt;

  
  if (MakeFullHist) {
    // Get full catalog file name(s)
    string delim=",";
    vector<string> fnames;
    stringSplit(FullCat,delim,fnames);
    int Nfiles = fnames.size();
    
    /*// Set name of catalog containing the "full" set of z's
      string FullCatFile;
      if(Nfiles>1)
      cout <<"    Reading in "<< Nfiles <<" full catalog files"<<endl; 
      else
      FullCatFile = FullCat;*/
    // Read in Fits file containing column with all true redshifts
    // FULL CATALOG               
    
    for (int ic=0; ic<Nfiles; ic++) {
      
      string FullCatFile = fnames[ic];
      cout <<"    Read in full catalog redshifts from "<< FullCatFile;
      cout <<" from column labelled "<< ZFCol <<endl;
      FitsInOutFile fin(FullCatFile,FitsInOutFile::Fits_RO);
      fin.MoveAbsToHDU(2);
      SwFitsDataTable dt(fin,512,false);
      cout <<endl;
      
      sa_size_t ncat = dt.NEntry();
      sa_size_t Izf = dt.IndexNom(ZFCol);
      dt.GetMinMax(Izf,minzt,maxzt); 
      DataTableRow rowin = dt.EmptyRow();
      cout <<"    Number of galaxies in this FULL catalog = "<< ncat <<endl;
      cout <<"    Min z of this FULL catalog = "<<minzt;
      cout <<", max z of this FULL catalog = "<<maxzt<<endl;
      
      cout <<"    Add to Histogram ... "<<endl;
      for(sa_size_t i=0; i<ncat; i++) { 
	dt.GetRow(i, rowin);
	double zs=rowin[Izf];
	nzFC.Add(zs);
      }
    }
    sa_size_t ngFC = nzFC.NEntries();
  }
  else {
    string sffile = FullCat + "_ZONLY.txt";
    ifstream inp;
    inp.open(sffile.c_str(), ifstream::in);
    inp.close();
    if(inp.fail()) { 
      // sffile does NOT exist
      string emsg = "ERROR! Selection function in file " + sffile;
      emsg += " does not exist";
      throw ParmError(emsg);
    }
    else {
      // sffile DOES exist
      cout <<" Histogram from all galaxies will be read from file " << sffile.c_str() <<endl;
      SInterp1D histo_full;
      vector<double> xsv, ysv;// the two columns to be read in
      
      size_t cnt=0; // count number of lines
      double cola, colb;
      inp.open(sffile.c_str(), ifstream::in);
      while(!inp.eof()) {  
        inp.clear();
        inp >> cola >> colb; // read each line into cola and colb
        cout <<" histo ZONLY = "<< cola <<", "<< colb << endl;
        
	//       if ( (!inp.good()) || inp.eof() ) {   
        if ( inp.eof() ) {   
	  break; // read until end of file
	}
	nzFC.Add(cola,colb); 
	if (cnt == 0) 
	  if (cola != minz + (maxz-minz)/nbin/2.) {
	    // sffile does NOT exist
	    string emsg = "ERROR! Histogram in file " + sffile;
	    emsg += " does not have the right minimum redshift";
	    throw ParmError(emsg);	
	  }
	
        cnt++;
      }
      inp.close();
      
      if (cola != maxz - (maxz-minz)/nbin/2.) {
	// sffile does NOT exist
	string emsg = "ERROR! Histogram in file " + sffile;
	emsg += " does not have the right maximum redshift";
	throw ParmError(emsg);	
      }
    }
  }
    
  // OBS CATALOG - loop over and add to Histo
  cout <<"    Histogram up OBSERVED catalog"<<endl;
 
  vector<string> fobsnames;
  stringSplit(ObsCat_,delim,fobsnames);
  int Nobsfiles = fobsnames.size();
  double x=1e8,y=1e8,z=-1e8,redshift=-10;
  for (int ic=0; ic<Nobsfiles; ic++) {
    
    string ObsCatFile = fobsnames[ic];
    cout <<"    Read in obs catalog redshifts from "<< ObsCatFile;
    cout <<" from column labelled "<< ZSCol <<" and " << ZOCol <<endl;
    FitsInOutFile fin(ObsCatFile,FitsInOutFile::Fits_RO);
    fin.MoveAbsToHDU(2);
    SwFitsDataTable dtobs(fin,512,false);
    cout <<endl;
    
    sa_size_t ncat = dtobs.NEntry();
    sa_size_t Izs = dtobs.IndexNom(ZSCol);
    sa_size_t Izp = dtobs.IndexNom(ZOCol);
    dtobs.GetMinMax(Izs,minzt,maxzt); 
    DataTableRow rowino = dtobs.EmptyRow();
    GalRecord grec;
    cout <<"    Number of galaxies in this OBS catalog = "<< ncat <<endl;
    cout <<"    Min z of this OBS catalog = "<<minzt;
    cout <<", max z of this OBS catalog = "<<maxzt<<endl;
    
    cout <<"    Add to Histogram ... "<<endl;
    for(sa_size_t ig=0; ig<ncat; ig++) { 
      dtobs.GetRow(ig, rowino);
      Row2Record(rowino,grec);
      if (!Filter(grec)) continue; 
      // not significantly slower than reading the redshifts and allow more possibilities
      if (RadialZ_ == true) 
	Rec2ShellCoord(grec,x,y,z,redshift);
      else 
	Rec2EuclidCoord(grec,x,y,z,redshift);
      
      nzOC.Add(redshift); // histo up OBSERVED -z
      nzOCsz.Add(grec.zs);// histogram up spec-z no matter what
      
      // print out first 10 gals
      if(ig<10) {
        cout <<"    galid="<< ig <<": ";
        grec.Print();
        if (AddGaussErr_)
          cout <<", zpG="<< redshift <<endl;
        else
          cout <<endl;
      }
    }
  }
 
    
  if (DoDebug_) {
      
    string outfile3 = SFTextFile+"_histo.ppf";
    // Let's keep the three histograms 
    POutPersist poh(outfile3);
    poh << PPFNameTag("nzFC") << nzFC;
    poh << PPFNameTag("nzOC") << nzOC;
    poh << PPFNameTag("nzOCsz") << nzOCsz;
    cout << "    Cat2Grid::SaveSelecFunc/Info nzFC,nzOC,nzOCsz";
    cout << " histos saved to file " << outfile3 << endl;
  }
    
    
  //READ HISTOGRAM VALUES INTO A FILE
  string outfile = SFTextFile+"_nofz.txt";
  string outfile2 = SFTextFile+"_specz_nofz.txt";
  string outfile0;
  if (MakeFullHist)
    outfile0 = SFTextFile+"_ZONLY.txt";
    
  cout <<"    Write n^o(z)/n^t(z) histograms to text files "<<endl;
  ifstream inp;
  ofstream outp,outp2,outp0;
    
  inp.open(outfile.c_str(), ifstream::in);
  inp.close();
  if(inp.fail()) {
      
    inp.clear(ios::failbit);
    cout << "    Writing to file ..." << outfile.c_str() << " and ";
    cout << outfile2 <<endl;
  if (MakeFullHist)
      cout << "  ...  and " << outfile0.c_str() <<endl;

    outp.open(outfile.c_str(), ofstream::out);
    outp2.open(outfile2.c_str(), ofstream::out);
    if (MakeFullHist)
      outp0.open(outfile0.c_str(), ofstream::out);

    for(int_4 i=0;i<nbin;i++) {
      r_8 bc=nzFC.BinCenter(i);
      r_8 bc2=nzFC.BinCenter(i);
      r_8 nzt=nzFC.operator()(i);
      r_8 nzo=nzOC.operator()(i);
      r_8 nzosz=nzOCsz.operator()(i);
        
      if((bc-bc2)>0.000002)
        cout <<"DIFFERENCE BETWEEN BIN CENTERS = "<< bc-bc2 <<endl;
        
      r_8 sf1=nzo/nzt; // phot-z (if ztype is phot)
      r_8 sf2=nzosz/nzt;// spec-z
        
      if (nzo<1&&nzt<1)// ie if they are both 0
        sf1=1;
      if (nzo<1&&nzt>0)// if obs is 0
        sf1=0.000001; //vsmall not zero, otherwise weight=1/SF=INF
      if (nzo>0&&nzt<1)// if true is 0 (can get this with photo-z scattering out of sim volume)
        sf1=50; //vlarge, i.e. gal downweighted to zero weight
        
      if (nzosz<1&&nzt<1)// ie if they are both 0
        sf2=1;
      if (nzosz<1&&nzt>0)// if obs is 0
        sf2=0.000001; //vsmall not zero, otherwise weight=1/SF=INF
      if (nzosz>0&&nzt<1)// if true is 0 (should not be possible if using spec-z as ztype)
        sf2=50; //vlarge, i.e. gal downweighted to zero weight
        
      //                        // if bin is outside full catalog's redshift range
      //                        if (bc<minzt||bc>maxzt)
      //                                {
      //                                sf1=1000000;//vlarge, anthing here is downweighted to zero weight
      //                                sf2=1000000;
      //                                }
        
      outp << bc <<"      "<< sf1 <<endl;
      outp2 << bc <<"      "<< sf2 <<endl;
      cout << " HISTO " << bc <<"  " << nzo <<endl;
      if (MakeFullHist)
	outp0 << bc <<"      "<< nzt <<endl;
      // cout << bc <<"      "<< sf1 <<endl;
    }
    outp.close();
    outp2.close();
    if (MakeFullHist)
      outp0.close();
  }
  else
    cout << "Error...file """ << outfile.c_str() << """ exists" << endl;
    
    
  tm.Split();
  cout <<"    Elapsed time "<< tm.TotalElapsedTime() <<endl;
  cout <<"    EXIT Cat2Grid::SaveSelecFunc()"<<endl<<endl;
};


void Cat2Grid::GalGrid(double SkyArea)
// Lays grid over data and counts the number of galaxies in each cell
// Cartesian z-direction will be radial direction from observer to center 
// of survey area
{
  cout <<endl<<"    Cat2Grid::GalGrid()"<<endl; 
  Timer tm("GalGrid");
  // Initialise random generator so specified seed (fixed or randomly generated)
        
  double minz = 0, maxz=5;
  int nbin=100;
  ngalz_.ReSize(minz, maxz, nbin); // z distribution of obs catalog   
  ngalzw_.ReSize(minz, maxz, nbin); // z distribution of weighted obs catalog 

  if (ErrRandomSeed_) {
    rg_.AutoInit(0);
    cout << "Seed automatically generated" << endl;
  } else {
    long seed=1;
    rg_.SetSeed(seed);
  }
        
  // Set sky area
  SkyArea_ = SkyArea;
        
  // Check grid has been set
        
  Npix_ = Nx_*Ny_*Nz_;
  volgrid_ = Npix_*pow(cellsize_,3);
  cout <<"    Volume of grid="<< volgrid_ <<endl;
                
  // NUMBER OF GALAXIES IN EACH CELL

  cout <<"    Cells in each direction: Nx="<< Nx_ <<", Ny="<< Ny_ <<", Nz=" << Nz_ <<endl;
  cout <<"    Grid boundaries: "<< Xmin_ <<"<X<"<< Xmax_ <<", " << Ymin_ <<"<Y<"<< Ymax_ <<", "<< Zmin_ <<"<Z<"<< Zmax_ <<endl;
  cout <<"    Grid cell size = "<< cellsize_ <<endl;
  cout <<"    Start loop over galaxies..."<<endl;
  cout <<"    applying filters ..."<<endl;
  cout <<"    Reading in OBSERVED redshifts from column "<< ZOCol_ <<endl;
  if (AddGaussErr_)
    cout <<"    Adding Gaussian photo-z error of size "<< PZDerr_ <<" Mpc to these"<<endl;
  //for photo-z error (Adeline)
  if (AddPhotoErrReds_){
    cout <<"    Adding photo-z error  "  <<endl;
    if(gRandom) 
      delete gRandom;
    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);
  }
        
  //Loop over ObsCat list of files - Cecile & Adeline modification
        
  string delim=",";
  vector<string> fobsnames;
  stringSplit(ObsCat_,delim,fobsnames);
  int Nobsfiles = fobsnames.size();

  // Take each galaxy and find which cell it lies in
  GalRecord grec;
  double x=1e8,y=1e8,z=-1e8,redshift=-10;
  double invselfunc=1.;
  double bias=0.;
  ngo_ = 0;
  ng_ = new sa_size_t[vgrids_.size()];
  ngout_ = new sa_size_t[vgrids_.size()];
  ngw_ = new r_8[vgrids_.size()];
  meanz_ = new r_8[vgrids_.size()];
  for (int igd=0; igd < vgrids_.size(); igd ++) {
    ng_[igd] = 0;
    ngout_[igd] = 0;
    ngw_[igd] = 0.;
    meanz_[igd] = 0.;
 }
  for (int ic=0; ic<Nobsfiles; ic++) {
    string ObsCatFile = fobsnames[ic];
    cout << "Read file #" <<  ic <<": "<< ObsCatFile << endl;
    FitsInOutFile fin(ObsCatFile,FitsInOutFile::Fits_RO);
    fin.MoveAbsToHDU(2);
    SwFitsDataTable dtobs(fin,512,false);
 
    sa_size_t n_gal = dtobs.NEntry();
    ProgressBar pgb(n_gal, ProgBarM_Time);  // ProgBarM_None, ProgBarM_Percent, ProgBarM_Time

    DataTableRow rowin = dtobs.EmptyRow();
    for(long ig=0; ig<n_gal; ig++) {
      // get row values from data table
      dtobs.GetRow(ig, rowin); 
            
      // add row values to GalRecord class
      Row2Record(rowin,grec);
      // add possible filtering on e.g. magnitude here
      // no filtering is implemented yet
      if (!Filter(grec)) continue; 
            
      // count galaxy as observed 
      if(!AddPhotoErrReds_)
        ngo_++; 
      else if (redshift != 10)
        ngo_++; 

      // convert galaxy position into Euclid coord - it is here that z error really matters
      if (RadialZ_ == true) 
	Rec2ShellCoord(grec,x,y,z,redshift);
      else 
	Rec2EuclidCoord(grec,x,y,z,redshift);
            
      // return selection function value at gal redshift
      invselfunc = 1./(*selfuncp_)(redshift);
            
      // basic check that selection function value isn't crazy
      if (invselfunc<0)
        cout <<"phi<0, invselfunc ="<< invselfunc <<", z="<< redshift <<endl;
            
      // print out first 10 observed gals, just as a basic check
      if(ig<10) {
        cout <<"    galid="<< ig <<": ";
        grec.Print();
        if (AddGaussErr_)
          cout <<", zpG="<< redshift <<", 1/selfunc="<< invselfunc <<endl;
        else if(AddPhotoErrReds_)
          cout <<", type="<<grec.type<<", MA="<<grec.MB <<", redshift="<< redshift <<endl;
        else
          cout <<endl;
      }
            
      //fill z distribution histo
      ngalz_.Add(redshift);
      ngalzw_.Add(redshift,invselfunc);

      // add galaxy to correct grid pixel 
      int Ngrid_per_gal = 0;
      for (size_t igd=0; igd<vgrids_.size(); igd++) {
	if (vgrids_[igd].AddGal(x,y,z,invselfunc)) {
	  cout << " gal OK " <<endl;
	  ngw_[igd] += invselfunc; 
	  ng_[igd] ++;
	  Ngrid_per_gal ++;
	  meanz_[igd] += redshift;
	} else ngout_[igd] ++;

      }
      pgb.update(ig); 
      if (Ngrid_per_gal > 1) {
	cout <<"    galid="<< ig <<": ";
        grec.Print();
	cout << "falls in " << Ngrid_per_gal << "grids."<< endl;
	throw ParmError("ERROR! this galaxy falls in more than 1 grid");
      }
    }


    cout <<"    Number of galaxies actually observed = "<< ngo_;
    cout <<", TOTAL number of galaxies in catalog = "<< n_gal <<endl;
    for (size_t igd=0; igd<vgrids_.size(); igd++) {
      cout <<  " GRID n." << igd << endl;
      cout <<"    Mean redshift of the galaxies in the grid = " << meanz_[igd] / ng_[igd] << endl;
      cout <<"    Number of galaxies inside = "<< ng_[igd] <<endl;
      cout <<"    Number of 'galaxies' inside WEIGHTED = "<< ngw_[igd] <<endl;
      cout <<"    Number of galaxies outside  = "<< ngout_[igd] <<endl;
    }
      
    cout <<"============== Analysis of catalog #" << ic << " done  ================================="<< endl << endl;
  }

  if (ngo_>ngall_)
    throw ParmError("ERROR! galaxy number discrepancy: ngo_!=ngall_");

  // write the arrays here!
  WriteGalArrays();

  tm.Split();
  cout <<"    Elapsed time "<< tm.TotalElapsedTime() <<endl;
  cout <<"    EXIT Cat2Grid::GalGrid()"<<endl<<endl; 
  
}

void Cat2Grid::computeMeanNg_z(double NormNgalMean)
{
  vector<double> zbin(ngalz_.NBins());
  vector<double> mean_ng(ngalz_.NBins());
  for (int i=0; i <(ngalz_.NBins()); i++) {
    double z=ngalz_.BinCenter(i);
    zbin[i] = z;
    double n=ngalz_(i);
    su_.SetEmissionRedShift(z);
    double d = su_.RadialCoordinateMpc();

    // shell thickness
    su_.SetEmissionRedShift(z-ngalz_.BinWidth()/2.);
    double d1 = su_.RadialCoordinateMpc();
    su_.SetEmissionRedShift(z+ngalz_.BinWidth()/2.);
    double d2 = su_.RadialCoordinateMpc();
    // shell surface & volume
    double surf = PI * d*d * 2. * (1. - cos(SkyArea_));
    double dvol = surf * (d2 - d1);

    mean_ng[i] = n / dvol * cellsize_ * cellsize_ * cellsize_ * NormNgalMean;
  }
  ngalz_cell_.DefinePoints(zbin, mean_ng);

}

void Cat2Grid::Row2Record(DataTableRow& rowin, GalRecord& rec)
{
  
  rec.xcoo = -9999 ; rec.ycoo= -9999;
  rec.alpha = rowin[Ic1_];   rec.delta=rowin[Ic2_];
  
  rec.zs = rowin[Izs_];  
  rec.zo = rowin[Iz_];
  rec.type = rowin[Itype_];
  rec.MB = rowin[Imag_];
  //    cout << "Row2Record "<< rec.zs << " " <<  rec.zo << endl;
  // rec.zs = SPECTRO Z
  // rec.zo = OBSERVED Z (could be spec-z, gauss-z, photo-z)
  
  // Fill in the magnitudes if available 
};


bool Cat2Grid::Filter(GalRecord& rec)
// Cut on magnitudes if applicable
{
  return true;
};

void Cat2Grid::readRootFile(string pdfFileName){
  /*
  ** Read root file and save TF1 function in an array
  */
  cout << "read  root file : "<< pdfFileName << " from the redshift slice "<< nz_0_ << endl;
  PDFrootFile_ =  new TFile(pdfFileName.c_str(), "open");
  //ZBIN    photozErrorArray = new TF1***[nz];
  int nz_reduced = 5;
  photozErrorArray = new TF1***[nz_reduced];
  string funcName;
  //ZBIN    for (int iz=0; iz<nz; iz++){
  for (int iz=0; iz<nz_reduced; iz++){
    photozErrorArray[iz] = new TF1**[ntype];

    for (int it=0; it<ntype; it++){
      photozErrorArray[iz][it] = new TF1*[nmag];
      for (int imag=0; imag<nmag; imag++){
	photozErrorArray[iz][it][imag] = NULL;
        funcName = Form( "fpdf%d_%d_%d", iz+nz_0_, it, imag );
        photozErrorArray[iz][it][imag] = (TF1*) PDFrootFile_ -> Get(funcName.c_str());
      }
    }
  }
}

double Cat2Grid::ComputeErrorFromPhotoz(GalRecord& rec){
  int iz, itype, imag, iz_reduced;
  double zs, type, mag, prob;
  bool plot_all_pb=false;
  zs = rec.zo;
  type = rec.type;
  mag = rec.MB;
  //cout << " zs "<< zs << "  type "<< type << " MB " << mag <<endl;
                
  //ZBIN  iz = int( (zs - z_min)/binz + 0.5);

  iz = int( (zs - z_min )/binz + 0.5);
  iz_reduced = int( (zs - (nz_0_*binz + z_min) )/binz + 0.5);
  itype = int(type-1);
  imag = int( (mag-mag_min)/binMag +0.5);
   
  if(iz>nz || itype>ntype || imag>nmag)         return 50.;
  
  if(photozErrorArray[iz_reduced][itype][imag] == NULL) return 60.;
  
  if (qualcut){
    double prob = QualCutProb[iz][itype][imag];
    
    //cout << iz << " " << itype << " " << imag << " " << prob << endl;
    
    if (prob==0)
      return 20.; 
    if (prob < 1){
      double rand_prob = rand.Uniform();
      if (rand_prob>prob)
	return 10.;
    }
  }
  return zs + photozErrorArray[iz_reduced][itype][imag] -> GetRandom();
}

void Cat2Grid::Rec2EuclidCoord(GalRecord& rec, double& x, double& y, double& z, double& redshift)
// From GalRecord object which contains:
// double alpha,delta,glong,glat,xcoo,ycoo,zcoo;
// double zs,zo;
// int type;
// double u,g,r,i,z,y;
// Calculates: x,y,z comoving coordinates, and returns which redshift is being analysed
// Choice of 3 redshifts: photometric, spectroscopic and spectroscopic + Gaussian error
{

  redshift = rec.zo;// zo is redshift to be used in analysis (could be spec-z,phot-z,gauss-z)
  
  // error in redshift will be translated in change in x,y,z. So to do first.
  if (AddGaussErrReds_) {
    redshift += PZerr_ * (1.+redshift) * rg_.Gaussian();
  }

  if(AddPhotoErrReds_){
    redshift = ComputeErrorFromPhotoz(rec);
  }

  if (doBiasCorr_) 
    redshift += (*biasp_)(redshift) * (1.+redshift);


  double dc = z2dist_(redshift);
  double zref;

  double ph=rec.alpha;
  double th=rec.delta;
  if(std::isnan(th)) {
    cout <<"    Cat2Grid::Rec2EuclidCoord/PB theta=nan -> set theta=0"<<endl;
    th=0;
  }
  
  x=dc*cos(ph)*sin(th);
  y=dc*sin(ph)*sin(th);
  z=dc*cos(th);  
  //Photo-z error (Adeline)
  //redshift = 10 means galaxies don't pass the quality cut 
};


void Cat2Grid::Rec2ShellCoord(GalRecord& rec, double& x, double& y, double& z, double& redshift)
// From GalRecord object which contains:
// double alpha,delta,glong,glat,xcoo,ycoo,zcoo;
// double zs,zo;
// int type;
// double u,g,r,i,z,y;
// Calculates: x,y,z comoving coordinates, and returns which redshift is being analysed
// Choice of 3 redshifts: photometric, spectroscopic and spectroscopic + Gaussian error
{

  // x,y must computed with the spectroZ : not affected by the redshift error (Cecile)
  redshift = rec.zs;
  double dc = z2dist_(redshift);
  double ph=rec.alpha;
  double th=rec.delta;
  //Commenter en attendant Cecile
  //  if(isnan(th)) {
  //    cout <<"    Cat2Grid::Rec2ShellCoord/PB theta=nan -> set theta=0"<<endl;
  //    th=0;
  //  }
  double r = dc*tan(th);
  r = th * dc;
  x = r*cos(ph);
  y = r*sin(ph);
  
  // zo is redshift to be used in analysis (could be spec-z,phot-z,gauss-z). It can be spectro or already with error or with error computed bellow
  redshift = rec.zo;
  z = z2dist_(redshift);

  // error in redshift won't change angles. So to do after x,y computation.
  if (AddGaussErrReds_) {
    redshift += PZerr_ * (1.+redshift) * rg_.Gaussian();
  }

  if(AddPhotoErrReds_){
    redshift = ComputeErrorFromPhotoz(rec);
  }

  if (doBiasCorr_) 
    redshift += (*biasp_)(redshift) * (1.+redshift);
};


void Cat2Grid::GetCellCoord(sa_size_t i, sa_size_t j, sa_size_t k, double& x, double& y, double& z)
// Given grid pixel index (i,j,k) return comoving cartesian coordinate of pixel center (x,y,z)
{

  x=GetCellX(i);
  y=GetCellY(j);
  z=GetCellZ(k);
  return;
};


void Cat2Grid::RandomGrid(double NormNgalMean, bool SaveArr, bool seed, bool SigRandom)
// Compute random weighted grid with same selection function as data
{
  ErrRandomSeed_  = seed;
  cout <<endl<<"    Cat2Grid::RandomGrid()"<<endl;
  Timer tm("RandomGrid");
  
  if (!sfcompute_)
    cout <<"    Selection function should be CONSTANT with z, check this ..."<<endl;
  else
    cout <<"    Check selection function .... "<<endl;
  cout <<"    Only filling pixels with theta <= "<< SkyArea_ <<endl;
  
  
  if (ErrRandomSeed_) {
    rg_.AutoInit(0);
    cout << "Seed automatically generated for random grid" << endl;
  } else {
    long seed=1;
    rg_.SetSeed(seed);
  }          

  
  double NgalWObsTot = 0.;
  for(size_t i=0; i<vgrids_.size(); i++) NgalWObsTot += ngw_[i];
  double NgalWObsPerCell = NormNgalMean * NgalWObsTot /  (double)Npix_ / (double)vgrids_.size();
  cout <<"    Random catalog mean density = "<< NgalWObsPerCell <<endl;    

  computeMeanNg_z(NormNgalMean);

  for (int i=0;i<ngalz_.NBins(); i++) 
    if (ngalz_(i) > 0) { 
      cout << "For ngalz_ bin "<< i << " (z = " << ngalz_.BinCenter(i) << ")  n gal=" << ngalz_(i) << ", n gal weighted="<< ngalzw_(i) ;
      cout << "  and normalized ng per cell=" << ngalz_cell_(ngalz_.BinCenter(i)) << endl;
    }
 

  for(size_t iv=0; iv<vgrids_.size(); iv++) {
    if (wrgals_.NbDimensions()<2) {
      int ndim=3;
      sa_size_t mydim[ndim];
      mydim[0]=Nx_; mydim[1]=Ny_; mydim[2]=Nz_;
      wrgals_.SetSize(ndim, mydim);       // weighted galaxy density field: random catalog
    }
  
    int Nmargin = 0;
    TArray<r_8> superwrgals;
    if (SigRandom) {
      int ndim=3;
      sa_size_t superdim[ndim];
      double margin = 400.; // 400 = 8 x 50 Mpc with 50 Mpc roughly 1 sigma for err = 0.03(1+z) so 4 sigma on each side
      Nmargin = (int)(margin / cellsize_ + 0.5) *2;
      cout << Nmargin << endl;
      superdim[0]=Nx_ + Nmargin; superdim[1]=Ny_ + Nmargin; superdim[2]=Nz_ + Nmargin;
      superwrgals.SetSize(ndim, superdim);       // weighted galaxy density field: random catalog
      cout << "Random grid will be extracted from a super-grid of size "<< superdim[0] <<"x"<< superdim[1] <<"x"<< superdim[2] <<" pixels" <<endl;
    } 
    
    if (SaveArr) {
      int ndim=3;
      sa_size_t mydim[ndim];
      mydim[0]=Nx_; mydim[1]=Ny_; mydim[2]=Nz_;
      zc_.SetSize(ndim, mydim); // redshifts of pixel centers
    }
    
    cout <<"    Compute weights and random catalog ..."<<endl;
    ProgressBar pgb((Nx_ + Nmargin)*(Ny_ + Nmargin), ProgBarM_Time);  // ProgBarM_None, ProgBarM_Percent, ProgBarM_Time
    size_t ccnt=0;
    int print_ing=0;
    int offsetX_rdm = ceil((Nx_ + Nmargin)/2. -1./2.);
    int offsetY_rdm = ceil((Ny_ + Nmargin)/2. -1./2.);
    int offsetZ_rdm = ceil((Nz_ + Nmargin)/2. -1./2.);

    int Xmin_rdm =  -ceil((Nx_ + Nmargin)/2)*cellsize_;	   
    int Ymin_rdm =  -ceil((Ny_ + Nmargin)/2)*cellsize_;	   
    int Zmin_rdm =  DCref_ - ceil((Nz_ + Nmargin)/2)*cellsize_;
 
    for(int i=0; i<Nx_ + Nmargin; i++) {
      for(int j=0; j<Ny_ + Nmargin; j++) {
	for(int k=0; k<Nz_ + Nmargin; k++)  {
	  
	  double z_rdm,rx,ry,rz;
	  double rphi,rtet,dc,zoverr;
	  double xc,yc,zc;
	  double selfunc, redshift_rdm, dcell_rdm,xc_rdm,yc_rdm,zc_rdm;
	  sa_size_t ix, jy, kz;

	  // Average number of gals expected in cell
	  xc = Xmin_rdm + i*cellsize_ + cellsize_/2 ;
	  yc = Ymin_rdm + j*cellsize_ + cellsize_/2;
	  zc = Zmin_rdm + k*cellsize_ + cellsize_/2;
	  double dcell =  ((RadialZ_) ? zc : sqrt(xc*xc+yc*yc+zc*zc));
	  double redshift = dist2z_(dcell);
	  uint_8 npoiss = rg_.PoissonAhrens(ngalz_cell_(redshift)); // Poisson fluctuate

	  if (SaveArr && iv == vgrids_.size()-1) 
	    zc_(i,j,k) = redshift;

	  if(SigRandom) { // case error on the redshift
	    for(int ing=0; ing<npoiss; ing++) { // from gal 1 to gal npoiss-1 in this cell...
	      if (print_ing <10)  cout << "D avant " << xc << " " << yc << " " << zc << endl;
	      // no random inside the cell: would add an additionnal error
	      // theta, phi using "true" position
	      double phi=atan2(yc,xc);
	      double zoverr=zc/dcell;
	      double theta=acos(zoverr);
	      // add an error on the redshift
	      z_rdm = redshift + (1.+redshift) * rg_.Gaussian() * ((*sigrp_)(redshift));
	      // compute new cell using theta, phi, z_rdm
	      dc = z2dist_(z_rdm);
	      rx=dc*cos(phi)*sin(theta);
	      ry=dc*sin(phi)*sin(theta);
	      rz=dc*cos(theta);  	      
	      ix=(sa_size_t)floor( rx/cellsize_)     +offsetX_rdm;
	      jy=(sa_size_t)floor( ry/cellsize_)     +offsetY_rdm;
	      kz=(sa_size_t)floor((rz-DCref_ -1)/cellsize_  )  +offsetZ_rdm;
	      if (print_ing <10)  {
		cout << "RANDOM SIG "<< redshift << "   Gaussian "<< ((*sigrp_)(redshift)) << " sig = " << ((*sigrp_)(redshift)) << " new z " << z_rdm ;
		cout << " old cell "<< i << " "<< j << " "<<  k << " et new cell "<<  ix << " "<<  jy << " "<<  kz << "(old d "<< dcell << ", new d " << dc << ")"<< endl;
	      }
	      if (print_ing <10) cout << cellsize_ << " "<<  DCref_ << " "<< offsetX_rdm<< " "<< offsetY_rdm<< " "<<  offsetZ_rdm << endl;
	      if ((ix>=0) && (ix<Nx_ + Nmargin) && (jy>=0) && (jy<Ny_ + Nmargin) && (kz>=0) && (kz<Nz_ + Nmargin))   { 
		selfunc = (*selfuncp_)(z_rdm);
		superwrgals(ix,jy,kz) += 1. / selfunc ;
	      } 
	      print_ing ++;

	    }
	    
	  } else { // case no error on the redshift
	    selfunc = (*selfuncp_)(redshift);
	    
	    wrgals_(i,j,k)= (double)npoiss / selfunc; 
	    
	  }
	  
	}   // end of loop over k (z-direction)
	ccnt++;  pgb.update(ccnt); 
      } // end of loop over j (y-direction)
    } // end of loop over  (x-direction)
    
    if(SigRandom) {
      // we extract the central cube of size NxNyNz from the super-cube with margin
      for(int i=0; i<Nx_; i++) {
	for(int j=0; j<Ny_ ; j++) {
	  for(int k=0; k<Nz_; k++)  {
	    wrgals_(i,j,k)= superwrgals(i + Nmargin/2,j + Nmargin/2,k + Nmargin/2);
	  }
	}
      }
    }
    wnrand_=wrgals_.Sum();
    cout <<"    Number of galaxies in RANDOM weighted grid = "<< wnrand_;
    cout <<", should be approximately input density*npixels = "<< NgalWObsPerCell*Npix_ <<endl;
    cout <<"    Calculated weighted random grid density = "<< wnrand_/Npix_;
    cout <<", input density = "<< NgalWObsPerCell <<" (should be approximately the same)"<<endl;
    VarianceRandomGrid();
    cout <<"    Mean of weighted random grid should be roughly "<< NgalWObsPerCell << endl;
    double minv,maxv;
    
    // write the arrays here!
    fos_ << wrgals_;
    
    if (SaveArr && iv == vgrids_.size()-1) {
      fos_ << zc_;
      cout << "WRITE Z array" << endl;
    }
  }
  // write histo of number of galaxies (weighted or not)
  fos_ << ngalz_;
  fos_ << ngalzw_;
  
  tm.Split();
  cout <<"    Elapsed time "<< tm.TotalElapsedTime()<<endl;
  cout <<"    EXIT Cat2Grid::RandomGrid()"<<endl<<endl; 
  
};

void Cat2Grid::SetGaussErrRedshift(double Err, double zref, bool seed) {
  AddGaussErr_    = true;
  AddGaussErrReds_= true;
  ErrRandomSeed_  = seed;
  PZerr_ = Err; 
  PZDerr_ = ZErr2CoDistErr(su_,Err,zref);
  cout <<"    Gaussian errors WILL be added to redshift "<< PZerr_   ; 
  if (ErrRandomSeed_) cout << " in Montecarlo mode" << endl<<endl; 
  else cout << endl<<endl; 
}

        
void Cat2Grid::SetPhotozErrRedshift(string pdfFileName, bool seed, double zcat_min, double zcat_max){
  InitParam(); 
  ifstream qualcutfile(Form("%s.txt",pdfFileName.c_str()));
  string header;
  qualcut = 1;
  
  for (int ihead = 0; ihead<3; ihead++){
    qualcutfile >> header;
    
    if (qualcutfile.eof() || header.find("#")==std::string::npos)
      cout << "parameters not defined in ... " << pdfFileName << ".txt" << endl;
    
    if (header.find("#z")!=std::string::npos){
      qualcutfile >> z_min >> z_max >> binz;
      //cout << "z found : " << z_min << " " << z_max << " " << binz << endl;
    }
    if (header.find("#type")!=std::string::npos){
      qualcutfile >> type_min >> type_max >> binType;
    }
    if (header.find("#mag")!=std::string::npos){
      qualcutfile >> mag_min >> mag_max >> binMag;
    }
  }
  
  nz = int( (z_max-z_min)/binz )+1;
  ntype = int( (type_max-type_min)/binType)+1;
  nmag = int( (mag_max-mag_min)/binMag)+1;
  
  
  cout << endl << "######################################" << endl;
  cout << "#Z   : " << nz << " " << z_min << " " << z_max << " " << binz << endl;
  cout << "#T   : " << ntype << " " << type_min << " " << type_max << " " << binType << endl;
  cout << "#MAG : " << nmag << " " << mag_min << " " << mag_max << " " << binMag << endl;
  cout << "######################################" << endl << endl;
  
  double allprob = 1, p;
  QualCutProb = new double**[nz];
  int iz=0, it=0, imag=0;
  for (iz=0; iz<nz; iz++){
    QualCutProb[iz] = new double*[ntype];
    for (it=0; it<ntype; it++){
      QualCutProb[iz][it] = new double[nmag];
      for (imag=0; imag<nmag; imag++){
	QualCutProb[iz][it][imag] = 0;
	//cout << "Set " << iz << " " << it << " " << imag << endl;
      }
    }
  }
  
  cout << "reading qualcutfile" << endl;
  
  while(1){
    qualcutfile >> iz;
    if (qualcutfile.eof())
      break;
    qualcutfile >> it >> imag >> p;
    //cout << "Read " << iz << " " << it << " " << imag << endl;
    allprob*=p;
    QualCutProb[iz][it][imag] = p;
  }
  if (TMath::Abs(allprob-1)<1e-2)
    qualcut = 0;
  qualcutfile.close();
  
  cout << allprob << " -> qualcut treatment = " << qualcut << endl;
  
  AddPhotoErrReds_ = true;
  ErrRandomSeed_ = seed;        
  pdfFileName_ = Form("%s.root", pdfFileName.c_str());
  //PDFrootFile_ =  new TFile(Form("%s.root", pdfFileName.c_str()), "open");
  int i_zcat_min = int( (zcat_min - z_min)/binz + 0.5);
  int i_zcat_max = int( (zcat_max - z_min)/binz + 0.5);
  cout << "**********************************************************************************************************************"<< endl;
  cout << "REDSHIFT range : from " << zcat_min << " to "<< zcat_max << " so from bin "<< i_zcat_min << " to "<< i_zcat_max << ". " << endl;
  cout << "*********************************************************************************************************************"<< endl;
  if (i_zcat_max-i_zcat_min > 4 ) {
    cout << "catalogue redshift is too large to support the photoZ method, please cut it."<< endl;
    exit(-1);
  }
  nz_0_ = i_zcat_min;
  readRootFile(pdfFileName_);
  rand.SetSeed(0);
  cout <<"    photo-z errors WILL be added to redshift, from root file "<<  pdfFileName_; 
  if (ErrRandomSeed_) cout << " in Montecarlo mode" << endl<<endl; 
  else cout << endl<<endl; 
  
}

/** Parameters over which pdf are computed (Adeline) */
void  Cat2Grid::InitParam(){                
  /*
    z_min = 0.1;
    z_max = 3.5;
    binz = 0.1;
    
    type_min = 1;
    type_max = 4;
    binType = 1;
    
    mag_min = -24;
    mag_max = -12;
    binMag = 0.2;
    
    nz = int( (z_max-z_min)/binz )+1;
    ntype = int( (type_max-type_min)/binType )+1;
    nmag = int( (mag_max-mag_min)/binMag )+1;
  */
    
  if(gRandom) 
    delete gRandom;
  gRandom = new TRandom3(0);
  gRandom->SetSeed(0);
}

void Cat2Grid::VarianceRandomGrid() {
  double mean,sigma;
  MeanSigma(wrgals_,mean,sigma);
  cout <<"    Weighted random array has mean = "<< mean;
  cout <<", variance = "<< sigma*sigma <<endl;
}


void  Cat2Grid::WriteHeader(string IncatName) {

  fos_.WriteKey("NX",Nx_," axe transverse 1");
  fos_.WriteKey("NY",Ny_," axe transverse 2");
  fos_.WriteKey("NZ",Nz_," axe longitudinal (redshift)");
  fos_.WriteKey("DX",cellsize_," Mpc");
  fos_.WriteKey("DY",cellsize_," Mpc");
  fos_.WriteKey("DZ",cellsize_," Mpc");
  fos_.WriteKey("ZREF",zref_," reference redshift");  
  fos_.WriteKey("ThRAD",SkyArea_," radius of circular sky area");
  fos_.WriteKey("MeanOverDensity",mean_overdensity_," mean dens SimLSS delta-fudge");
  fos_.WriteKey("InCat",IncatName," original catalog");
  fos_.WriteKey("NGALCAT",ngall_," N gals in InCat");
  fos_.WriteKey("NumberOfGrids",int(vgrids_.size())," number of grids");

  for(size_t i=0; i<vgrids_.size(); i++) {
    fos_.WriteKey("NGRID",ng_[i]," N gals in grid");
    fos_.WriteKey("NOGRID",ngout_[i]," N gals outside grid");
    fos_.WriteKey("NWGRID",ngw_[i]," N weighted gals in grid");   
  } 
  fos_.WriteKey("NRGRID",wnrand_," N (weighted if appl.) gals in random grid");
          
  //modified by Adeline : write cosmo parameters in file header
  fos_.WriteKey("H0", su_.H0()," Cosmo.Param H0");
  fos_.WriteKey("OMEGAM0", su_.OmegaMatter()," Cosmo.Param OmegaMatter0 ");
  fos_.WriteKey("OMEGAB0", su_.OmegaBaryon()," Cosmo.Param OmegaBaryon0");
  fos_.WriteKey("OMEGAR0", su_.OmegaRadiation()," Cosmo.Param OmegaRadiation0");
  fos_.WriteKey("OMEGAT0", su_.OmegaTotal()," Cosmo.Param OmegaTot0");
  fos_.WriteKey("OMEGADE0", su_.OmegaLambda(),"  Cosmo.Param OmegaLambda0 (dark energy density)");
  fos_.WriteKey("OMEGADK", su_.OmegaCurv(),"  Cosmo.Param OmegaK ");
  fos_.WriteKey("DE_W0", su_.wDE(), " Cosmo.Param w0 (dark energy eq.state)");
  fos_.WriteKey("DE_WA",su_.waDE() , " Cosmo.Param wA (dark energy eq.state)"); 
  fos_.WriteKey("SIGMA8", su_.Sigma8(), " Cosmo.Param sigma8_0");
  fos_.WriteKey("N_S",su_.Ns()," Cosmo.Param n_s (spectral index scalar fluct.)");
          
  cout << "Check cosmo parameters : " << endl;
  cout << "  OmegaK="<< su_.OmegaCurv() <<", OmegaM="<< su_.OmegaMatter();
  cout << ", OmegaL="<< su_.OmegaLambda() <<", OmegaB="<< su_.OmegaBaryon()  ;
  cout << ", Omega_rad=" << su_.OmegaRadiation() << ", H0=" << su_.H0() << ", Sig8=" << su_.Sigma8() <<", n_s=" << su_.Ns() <<endl; 
  cout << ", Omega_curv=" << su_.OmegaCurv() << ", DE_W0=" << su_.wDE() << ", DE_WA=" << su_.waDE() <<endl; 
  cout << endl;
  // end modifications
};


void Cat2Grid::WriteGalArrays() {
  if (vgrids_.size()>0) {
    cout << "Cat2Grid::WriteGalArrays() - writing vector of grids with NGrids="<<vgrids_.size()<<endl;
    for(size_t i=0; i<vgrids_.size(); i++) fos_ << vgrids_[i].getGrid();
  }
}

