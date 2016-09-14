#include "cat2grid.h"
#include "ctimer.h"
#include "progbar.h"




//Improved constructor:
Cat2Grid::Cat2Grid(SwFitsDataTable& dt, SimpleUniverse& su, RandomGenerator& rg,
                   FitsInOutFile& fos, string ZOCol, string ZSCol,
		   double PZerr, bool Print, string ObsCat) 
  : dt_(dt) , su_(su) , rg_(rg) , fos_(fos) , defselfunc_(1.), selfuncp_(&defselfunc_)
  , ZOCol_(ZOCol) , ZSCol_(ZSCol) , PZerr_(PZerr), ObsCat_(ObsCat)
                  // Constructor for this class reads in the phi, theta, redshift of all galaxies in the catalog
                  // from data table dt
                  // su:          class which holds the cosmological parameters and can be used to 
                  //        calculate cosmological quantities
                  // rg:          class used to generate random numbers
                  // ZOCol: col name of OBSERVED redshifts to read in, could be spec-z, gauss-z, phot-z
                  // ZSCol: col name of SPECTRO redshifts to be read in (must be spec-z)
                  // PZerr: size of photometric redshift errors: sigz = PZerr(1+z); 
                  // Print: if true prints some extra info to screen whilst in constructor
{
  if (Print)
    cout <<endl<<"    In Cat2Grid constructor ..."<<endl;

  DVList dvl;
  dvl = dt.Info();
  mean_overdensity_ = dvl["MeanOverDensity"];
  cout << "    Mean over-density of fudged SimLSS grid = "<< mean_overdensity_ <<endl;
        
  DoDebug_ = false; // flag to do debug
  debugoutroot_ = "tmp"; // if want to debug, files are written to filenames with this root
  sfcompute_ = false; // flag becomes true after selection function is set
  PZDerr_=0;// set photoz error in Mpc equal to zero to start with
  AddGaussErr_     = false; // flag to add Gaussian photo-z err
  AddGaussErrReds_ = false; // flag to add Gaussian photo-z err on redshift
  AddPhotoErrReds_ = false; // flag to add photo-z photo-z err on redshift

  // Set counting galaxies to ZERO
  ng_=0,ngo_=0; // num gals INSIDE grid, num gals OBSERVED
  ngw_=0;// num "gals" INSIDE WEIGHTED grid
  ngout_=0; // number of galaxies OUTSIDE grid
  ngzok_=0; // number of galaxies OUTSIDE grid
  wnrand_=0,nrand_=0; // weighted/unweighted number of "galaxies" in random grid

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


// make interpolation table outside Cat2Grid, just use Row2Record, Rec2ShellCoord functions
Cat2Grid::Cat2Grid(SwFitsDataTable& dt,SimpleUniverse& su,RandomGenerator& rg,
                   SInterp1D dist2z, SInterp1D z2dist,string ZOCol,string ZSCol)
		   : dt_(dt) , su_(su) , rg_(rg) , defselfunc_(1.), selfuncp_(&defselfunc_) , fos_(fosdefault_) ,
		   ZOCol_(ZOCol) , ZSCol_(ZSCol) , dist2z_(dist2z) , z2dist_(z2dist)
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
  // All in capital letters to avoid any pb with IDL
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
  selfuncp_=a.selfuncp_;
  ZSCol_=a.ZSCol_;
  ZOCol_=a.ZOCol_;

  //Initialised in FindMinMaxCoords 
  zsmin_=a.zsmin_;
  zsmax_=a.zsmax_;      
  xmin_=a.xmin_;
  xmax_=a.xmax_;
  ymin_=a.ymin_;
  ymax_=a.ymax_;
  zmin_=a.zmin_;
  zmax_=a.zmax_;
                 
  ng_=a.ng_;            
        
  //Initialised in ComputeSFfromLF
  /*phistar_=a.phistar_;
    Mstar_=a.Mstar_;
    alpha_=a.alpha_;
    mlim_=a.mlim_;
    Mc_=a.Mc_;*/
        
  //Initialised in ComputeSFfromLF or ComputeSF
  if (a.dL_.size()>0)
    {
      dL_=a.dL_;
      phi_=a.phi_;
      dc_=a.dc_;
    }
        
  //Initialised in ComputeSF
  if (a.zs_.size()>0)
    zs_=a.zs_;
        
  //Initialised in ConstructGrid
  Vol_=a.Vol_;  
  Nx_=a.Nx_;
  Ny_=a.Ny_;
  Nz_=a.Nz_;
  volgrid_=a.volgrid_;
  cellsizeXY_=a.cellsizeXY_;
  cellsizeZ_=a.cellsizeZ_;
        
  Npix_=a.Npix_;
  if (a.ngals_.Size()>0)
    ngals_=a.ngals_;
  if (a.randomcat_.Size()>0)
    {
      randomcat_=a.randomcat_;
      wngals_=a.wngals_;
      weights_=a.weights_;
      wrgals_=a.wrgals_;
    }
  alph_=a.alph_;
         
};


double Cat2Grid::FindMinMaxCoords()
// Finds min and max of cartesian x,y,z coords
// of all galaxies in simulation
// Find min/max redshift
{
  cout <<endl<<"    Cat2Grid::FindMinMaxCoords()"<<endl;
        
  if (AddGaussErr_)  {
    cout <<"    Remember: adding Gaussian errors, equivalent distance error = ";
    cout << PZDerr_ <<" Mpc"<<endl;
  }
                
  long seed=1; 
  rg_.SetSeed(seed);
        
  sa_size_t ng=0; // reset to zero to count up all OBSERVED galaxies
        
  // artificially set to very large min, very small max
  double minx=1e10,miny=1e10,minz=1e10;
  double maxx=0,maxy=0,maxz=0;
  double minzs=20,maxzs=0;
        
  DataTableRow rowin = dt_.EmptyRow();
  GalRecord grec;
  double x=1e8,y=1e8,z=-1e8,redshift=-10;
        
  ProgressBar pgb(ngall_, ProgBarM_Percent);
  for(long ig=0; ig<ngall_; ig++) {
                
    dt_.GetRow(ig, rowin);
    Row2Record(rowin,grec);
    if (!Filter(grec)) continue;
    Rec2ShellCoord(grec,x,y,z,redshift);
    // the redshift returned here is the redshift to be used
    // in the analysis, could be spec-z, phot-z, gauss-z
        
    // print out first 10 gals
    if(ig<10)   
      {
        cout <<"    galid test="<<ig<<": ";
        grec.Print();
        cout <<", x="<<x<<", y="<<y<<", z="<<z;
        if (AddGaussErr_)
          cout <<", zpG="<<redshift<<endl;
        else
          cout <<endl;
      } 
    ng++; 
                
    // perform check to find minimums/maximums  
    if(x<minx)
      minx = x;
    if(y<miny)
      miny = y;
    if(z<minz)
      minz = z;
                        
    if(x>maxx)
      maxx = x;
    if(y>maxy)
      maxy = y;
    if(z>maxz)
      maxz = z;
                        
    if (redshift<minzs)
      minzs = redshift;
    if (redshift>maxzs)
      maxzs = redshift;

    pgb.update(ig);
  }
                
  cout <<"    Number of observed galaxies = "<< ng;
  cout <<", Total number of galaxies = "<< ngall_ <<endl;
        
  if (ng>ngall_)
    throw ParmError("ERROR! galaxy number discrepancy");
        
  // set min and max to class variables
  xmin_ = minx;
  xmax_ = maxx;
  ymin_ = miny;
  ymax_ = maxy;
  zmin_ = minz;
  zmax_ = maxz;
  zsmin_ = minzs;
  zsmax_= maxzs;
                
  // find maximum (obs) luminosity distance in survey
  su_.SetEmissionRedShift(zsmax_);
  double maxdL=su_.LuminosityDistanceMpc();
                
  cout <<"    compute approx SURVEY volume .... "<<endl;
  cout <<"    xmin="<< xmin_ <<", xmax="<< xmax_ <<endl;
  cout <<"    ymin="<< ymin_ <<", ymax="<< ymax_ <<endl;
  cout <<"    zmin="<< zmin_ <<", zmax="<< zmax_ <<endl;
  Vol_ = (xmax_-xmin_)*(ymax_-ymin_)*(zmax_-zmin_);
  cout <<"    -> Volume = "<< Vol_ <<" Mpc"<<endl;
  cout <<"    OBSERVED redshift range (from column "<< ZOCol_ <<"): ";
  cout << zsmin_ <<"<z<"<< zsmax_ <<endl;
  cout <<"    Maximum dL="<< maxdL <<endl;
        
  cout <<"    EXIT Cat2Grid::FindMinMaxCoords()"<<endl<<endl;
        
  return maxdL;
        
};


void Cat2Grid::SetGrid(double Theta, double zl, double zh,double R_XY, double R_Z)
// Computes range of grid to lay over galaxy simulation
// from specified survey volume
{
  cout <<endl<<"    Cat2Grid::SetGrid()"<<endl;
  cellsizeXY_=R_XY;
  cellsizeZ_=R_Z;

  cout << "    Using survey geometry to define grid"<<endl;
  cout << "    Grid to cover: "<< zl <<"<z<"<< zh <<", Theta = "<< Theta <<" radians"<<endl;
  cout << "    Pixel size XY = "<< cellsizeXY_ <<endl;
  cout << "    Pixel size Z = "<< cellsizeZ_ <<endl;
  // To find minx,maxx,miny,maxy need to find size Theta extends
  // at the maximum redshift (zh).  Since x=y=0 at center of survey:
  su_.SetEmissionRedShift(zh);
  Xmin_ = -su_.RadialCoordinateMpc()*tan(Theta);
  Ymin_ = Xmin_;
  Xmax_ = -Xmin_;
  Ymax_ = -Ymin_;
  // Min z value is just distance to min redshift (zl) at CENTER of survey
  su_.SetEmissionRedShift(zl);
  Zmin_=su_.RadialCoordinateMpc();
  // Max z value is distance to max redshift at EDGE of survey
  su_.SetEmissionRedShift(zh);
  Zmax_=su_.RadialCoordinateMpc()/cos(Theta);
  cout <<"    CHECK Cartesian survey boundaries: "<< Xmin_ <<"<X<"<< Xmax_;
  cout <<", "<< Ymin_ <<"<Y<"<< Ymax_ <<", "<< Zmin_ <<"<Z<"<< Zmax_ <<endl;
        
  // calculate EXACT number of cells that reach from one side of sim to the other (including remainder)
  double Lx, Ly, Lz; 
  Lx = (Xmax_-Xmin_)/cellsizeXY_;
  Ly = (Ymax_-Ymin_)/cellsizeXY_;
  Lz = (Zmax_-Zmin_)/cellsizeZ_;
  cout <<"    Number of cells which cover survey in each direction"<<endl;
  cout <<"    Nx="<< Lx <<", Ny="<< Ly <<", Nz="<< Lz <<endl;
  // round DOWN to get INTEGER number of cells that cover each dimension
  // round down because don't want a pixel going outside of survey area
  // integer number cells needed
  Nx_ = (sa_size_t)floor(Lx); Ny_ = (sa_size_t)floor(Ly); Nz_ = (sa_size_t)floor(Lz); 
  cout <<"    Integer number of cells: Nx="<< Nx_ <<", Ny="<< Ny_ <<", Nz="<< Nz_ <<endl;
        
  double idzd = ceil((double)Nz_/2); // index of center pixel in z direction
  double idxd = ceil((double)Nx_/2);
  double idyd = ceil((double)Ny_/2);    
  idx_ = idxd, idy_ = idyd,idz_ = idzd;
  cout <<"    Indices of center pixels: Ix="<< idx_ <<", Iy="<< idy_ <<", Iz="<< idz_ <<endl;
  DCref_ = Zmin_+cellsizeZ_/2+(idz_-1)*cellsizeZ_;
  cout <<"CHECK: comoving distance = "<< DCref_ <<endl;
        
  // Need to recompute Xmin_ etc to EXACTLY match grid edges
  Xmin_ = -(idx_-1)*cellsizeXY_-cellsizeXY_/2;
  Xmax_ = (idx_-1)*cellsizeXY_+cellsizeXY_/2;
  Ymin_ = -(idy_-1)*cellsizeXY_-cellsizeXY_/2;
  Ymax_ = (idy_-1)*cellsizeXY_+cellsizeXY_/2;
  Zmin_ = DCref_-(idz_-1)*cellsizeZ_-cellsizeZ_/2;
  Zmax_ = DCref_+(idz_-1)*cellsizeZ_+cellsizeZ_/2;
        
  cout <<"    Initial survey boundaries: "<< Xmin_ <<"<X<"<< Xmax_;
  cout <<", "<< Ymin_ <<"<Y<"<< Ymax_ <<", "<< Zmin_ <<"<Z<"<< Zmax_ <<endl;

  //    // round down to nearest integer SHOULD REMOVE THIS
  //    Xmin_=floor(Xmin_);
  //    Xmax_=floor(Xmax_);
  //    Ymin_=floor(Ymin_);
  //    Ymax_=floor(Ymax_);
  //    Zmin_=floor(Zmin_);
  //    Zmax_=floor(Zmax_);
  //    cout <<"    Rounded survey boundaries: "<<Xmin_<<"<X<"<<Xmax_<<", "<<Ymin_<<"<Y<"<<Ymax_<<", "<<Zmin_<<"<Z<"<<Zmax_<<endl;

  DCref_ = Zmin_+R_Z/2+(idz_-1)*R_Z;
        
  // convert DCref_ to a redshift
  int_8 nz=1000;
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
  double mind = codist[0];
  double maxd = codist[codist.size()-1];
  SInterp1D dist2z(codist,zrs,mind,maxd,2*nz);
  zref_ = dist2z(DCref_);
  cout <<"    Redshift of central pixel = "<< zref_;
  cout <<", found from comoving distance = "<< DCref_ <<endl;
        
  // INITIALISE ARRAYS TO GRID SIZE
  int ndim=3;
  sa_size_t mydim[ndim];
  mydim[0]=Nx_; mydim[1]=Ny_; mydim[2]=Nz_;
  ngals_.SetSize(ndim, mydim);  // number of galaxies in each cell: data
  wngals_.SetSize(ndim, mydim); // weighted galaxy density field: data
  if (sfcompute_) {
    randomcat_.SetSize(ndim, mydim);// number of galaxies in each cell: random catalog
    weights_.SetSize(ndim, mydim);      // weights of galaxy density field (per cell)
    wrgals_.SetSize(ndim, mydim);       // weighted galaxy density field: random catalog
  }
        
  cout <<"    EXIT Cat2Grid::SetGrid()"<<endl<<endl<<endl;
};


void Cat2Grid::SetGrid(int_8 Nx, int_8 Ny, int_8 Nz, double R_XY, double R_Z, double zref)
// Computes range of grid to lay over galaxy simulation
// from SimLSS cube parameters
{
  cout <<endl<<"    Cat2Grid::SetGrid()"<<endl;
        
  cout << "    Using SimLSS grid to define grid"<<endl;
  cout << "    SimLSS grid: "<< Nx <<","<< Ny <<","<< Nz <<" pixels"<<endl;
  cout << "    Centered at redshift z = "<< zref <<endl;
  cout << "    Resolution R_XY = "<< R_XY <<endl;
  cout << "    Resolution R_Z = "<< R_Z <<endl;

  cellsizeXY_ = R_XY;
  cellsizeZ_ = R_Z;
  zref_ = zref;
        
  su_.SetEmissionRedShift(zref_);
  DCref_ = su_.RadialCoordinateMpc();
  double idzd = ceil((double)Nz/2); // index of center pixel in z direction
  double idxd = ceil((double)Nx/2);
  double idyd = ceil((double)Ny/2);     
  idx_ = idxd, idy_ = idyd,idz_ = idzd;
  cout <<"    Indices of center pixels: Ix="<< idx_ <<", Iy="<< idy_ <<", Iz="<< idz_ <<endl;
        
  Xmin_ = -(idx_-1)*cellsizeXY_-cellsizeXY_/2;
  Xmax_ = (idx_-1)*cellsizeXY_+cellsizeXY_/2;
  Ymin_ = -(idy_-1)*cellsizeXY_-cellsizeXY_/2;
  Ymax_ = (idy_-1)*cellsizeXY_+cellsizeXY_/2;
  Zmin_ = DCref_-(idz_-1)*cellsizeZ_-cellsizeZ_/2;
  Zmax_ = DCref_+(idz_-1)*cellsizeZ_+cellsizeZ_/2;
        
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
  Lx = (Xmax_-Xmin_)/cellsizeXY_;
  Ly = (Ymax_-Ymin_)/cellsizeXY_;
  Lz = (Zmax_-Zmin_)/cellsizeZ_;
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
  ngals_.SetSize(ndim, mydim);  // number of galaxies in each cell: data
  wngals_.SetSize(ndim, mydim); // weighted galaxy density field: data (same as ngal if no SF)
  if (sfcompute_) {
    randomcat_.SetSize(ndim, mydim);// number of galaxies in each cell: random catalog
    weights_.SetSize(ndim, mydim);      // weights of galaxy density field (per cell)
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
  cout <<"    Cat2Grid::SaveSelecFunc()"<<endl;
  cout << "OBSCAT : " << ObsCat_ << endl;
  Timer tm("SaveSelecFunc");
  
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
  if(! inp.fail()) {  
    cout << "Error...file """ << outfile.c_str() << """ exists" << endl;
    exit(-1);
  }
    
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
    double minzt, maxzt;
    
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
 
  string delim=",";
  vector<string> fobsnames;
  double minzt, maxzt;
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
      Rec2ShellCoord(grec,x,y,z,redshift);
      
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
      
  inp.clear(ios::failbit);
  cout << "    Writing to file ..." << outfile.c_str() << " and " << outfile2.c_str() <<endl;
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
    
    outp << bc <<"      "<< sf1 <<endl;
    outp2 << bc <<"      "<< sf2 <<endl;
    if (MakeFullHist)
      outp0 << bc <<"      "<< nzt <<endl;
    // cout << bc <<"      "<< sf1 <<endl;
  }
  outp.close();
  outp2.close();
  if (MakeFullHist)
      outp0.close();

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
  if (!ngals_.IsAllocated())
    throw ParmError("ERROR! Grid has not been set");
        
  Npix_ = Nx_*Ny_*Nz_;
  volgrid_ = Npix_*cellsizeXY_*cellsizeXY_*cellsizeZ_;
  cout <<"    Volume of grid="<< volgrid_ <<endl;

  // Find total number of observed pixels
  cout <<"    Finding number of observed pixels ... "<<endl;
  n_obs_pixels_ = ObsPixels();
  cout <<"    ... "<< n_obs_pixels_ <<" are observed out of "<< Npix_ <<" total pixels"<<endl;
                
  // NUMBER OF GALAXIES IN EACH CELL
  // first set array to zero
  ngals_=0.;

  cout <<"    Cells in each direction: Nx="<< Nx_ <<", Ny="<< Ny_ <<", Nz=";
  cout << Nz_ <<endl;
  cout <<"    Grid boundaries: "<< Xmin_ <<"<X<"<< Xmax_ <<", ";
  cout << Ymin_ <<"<Y<"<< Ymax_ <<", "<< Zmin_ <<"<Z<"<< Zmax_ <<endl;
  cout <<"    Grid cell XY size = "<< cellsizeXY_ <<endl;
  cout <<"    Grid cell Z size  = "<< cellsizeZ_ <<endl;
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
  double selphi=1.;

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

      // convert galaxy position into Shell coord - it is here that z error really matters
      Rec2ShellCoord(grec,x,y,z,redshift);
            
      // return selection function value at gal redshift
      selphi = (*selfuncp_)(redshift);
            
      // basic check that selection function value isn't crazy
      if (selphi<0)
        cout <<"phi<0, phi="<< selphi <<", z="<< redshift <<endl;
            
      // print out first 10 observed gals, just as a basic check
      if(ig<10) {
        cout <<"    galid="<< ig <<": ";
        grec.Print();
        if (AddGaussErr_)
          cout <<", zpG="<< redshift <<", phi="<< selphi <<endl;
        else if(AddPhotoErrReds_)
          cout <<", type="<<grec.type<<", MA="<<grec.MB <<", redshift="<< redshift <<endl;
        else
          cout <<endl;
      }
            
      // add galaxy to correct grid pixel 
      // ng_ and ngout_ are summed up in here:
      AddToCell(x,y,z,selphi); 
      pgb.update(ig); 
    }

    cout <<"    Number of galaxies actually observed = "<< ngo_;
    cout <<", TOTAL number of galaxies in catalog = "<< n_gal <<endl;
    cout <<"    Number of galaxies inside grid = "<< ng_ <<endl;
          
    cout <<"    Number of 'galaxies' inside WEIGHTED grid = "<< ngw_ <<endl;
    cout <<"    .... (should be about equal to the number of galaxies in FULL";
    cout <<" catalog - not printed here)"<<endl;
          
    cout <<"    Number of galaxies outside grid = "<< ngout_ <<endl;
    cout <<"============== Analysis of catalog #" << ic << " done  ================================="<< endl << endl;
  }

  if (ngo_>ngall_)
    throw ParmError("ERROR! galaxy number discrepancy: ngo_!=ngall_");
  if (ng_!=ngals_.Sum())
    throw ParmError("ERROR! galaxy number discrepancy: ng_!=ngals_.Sum()");
  if (ngo_!=ngals_.Sum())
    cout <<"    Number of galaxies outside grid = "<< ngout_ <<endl;

  // normalise the arrays
  NormNArrays();

  // write the arrays here!
  WriteGalArrays();

  tm.Split();
  cout <<"    Elapsed time "<< tm.TotalElapsedTime() <<endl;
  cout <<"    EXIT Cat2Grid::GalGrid()"<<endl<<endl; 

};


sa_size_t Cat2Grid::ObsPixels()
{

  sa_size_t nobspix = 0;
  double thetac;
  for(int i=0; i<Nx_;i++)
    for(int j=0; j<Ny_;j++)
      for(int k=0; k<Nz_;k++) {
	
        double xc,yc,zc;
        GetCellCoord(i,j,k,xc,yc,zc);
        double dcell  = zc;
        double thetac = atan(sqrt(xc*xc + yc*yc)/dcell);
	if (thetac<=SkyArea_)
	  nobspix++;
      }
  return nobspix;
};


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
  //  cout << " zs "<< zs << "  type "<< type << " MB " << mag <<endl;
                
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
  ngzok_++;
  return zs + photozErrorArray[iz_reduced][itype][imag] -> GetRandom();
}

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
  if(isnan(th)) {
    cout <<"    Cat2Grid::Rec2ShellCoord/PB theta=nan -> set theta=0"<<endl;
    th=0;
  }
  double r = dc*tan(th);
  x = r*cos(ph);
  y = r*sin(ph);
  
  // zo is redshift to be used in analysis (could be spec-z,phot-z,gauss-z). It can be spectro or already with error or with error computed bellow
  redshift = rec.zo;
  z = z2dist_(redshift);

  if (AddGaussErrReds_) {
    redshift += PZerr_ * (1.+redshift) * rg_.Gaussian();
  }

  if(AddPhotoErrReds_){
    redshift = ComputeErrorFromPhotoz(rec);
  }
};


void Cat2Grid::AddToCell(double x, double y, double z,double phi)
// Find which index galaxy with (x,y,z) coordinates lies in
// Fill ngals and wngals arrays accordingly
{

  sa_size_t indexx=(sa_size_t)floor((x-Xmin_)/cellsizeXY_);
  sa_size_t indexy=(sa_size_t)floor((y-Ymin_)/cellsizeXY_);
  sa_size_t indexz=(sa_size_t)floor((z-Zmin_)/cellsizeZ_);
 
  if (indexx<Nx_&&indexy<Ny_&&indexz<Nz_&&indexx>=0&&indexy>=0&&indexz>=0) {
        
    ngals_(indexx,indexy,indexz)++;
    if (sfcompute_) { // NOTE weight does not have to be 1/SF
      wngals_(indexx,indexy,indexz)+=1/phi;
      ngw_+=1/phi; }
    else                // AND if cat has no SF, wngals_=ngals_
      wngals_(indexx,indexy,indexz)++;
    ng_++;
  }
  else
    ngout_++;
};


void Cat2Grid::GetCellCoord(sa_size_t i, sa_size_t j, sa_size_t k, double& x, double& y, double& z)
// Given grid pixel index (i,j,k) return comoving cartesian coordinate of pixel center (x,y,z)
{

  x=GetCellX(i);
  y=GetCellY(j);
  z=GetCellZ(k);
  return;
};


void Cat2Grid::RandomGrid(double nc, bool SaveArr)
// Compute random weighted grid with same selection function as data
{
  cout <<endl<<"    Cat2Grid::RandomGrid()"<<endl;
  Timer tm("RandomGrid");
        
  cout <<"    Random catalog mean density = "<< nc <<endl;
  if (!sfcompute_)
    cout <<"    Selection function should be CONSTANT with z, check this ..."<<endl;
  else
    cout <<"    Check selection function .... "<<endl;
  cout <<"    Only filling pixels with theta <= "<< SkyArea_ <<endl;
        
                
  //    throw ParmError("ERROR! selection function has NOT been computed");
                
  if (wrgals_.NbDimensions()<2) {
    int ndim=3;
    sa_size_t mydim[ndim];
    mydim[0]=Nx_; mydim[1]=Ny_; mydim[2]=Nz_;
    randomcat_.SetSize(ndim, mydim);// number of galaxies in each cell: random catalog
    weights_.SetSize(ndim, mydim);      // weights of galaxy density field (per cell)
    wrgals_.SetSize(ndim, mydim);       // weighted galaxy density field: random catalog
  }
  if (SaveArr) {
    int ndim=3;
    sa_size_t mydim[ndim];
    mydim[0]=Nx_; mydim[1]=Ny_; mydim[2]=Nz_;
    zc_.SetSize(ndim, mydim); // redshifts of pixel centers
  }
  cout <<"    Compute weights and random catalog ..."<<endl;
        
  ProgressBar pgb(Nx_*Ny_, ProgBarM_Time);  // ProgBarM_None, ProgBarM_Percent, ProgBarM_Time
  size_t ccnt=0;
  for(int i=0; i<Nx_; i++) {
          
    //tm.Split();
    //cout<<"outer loop = "<<i<<", time elapsed = "<<tm.TotalElapsedTime()<<endl;
    for(int j=0; j<Ny_; j++) {
      for(int k=0; k<Nz_; k++)  {
              
        double xc,yc,zc;
        GetCellCoord(i,j,k,xc,yc,zc);
        double dcell = zc;
        double thetac = atan(sqrt(xc*xc + yc*yc)/dcell);
        double redshift = dist2z_(dcell);
        double selphi = (*selfuncp_)(redshift);
              
        if (thetac<=SkyArea_)   {
                
          weights_(i,j,k) = 1/selphi; // NOTE weight does not have to be 1/SF
          // Average number of gals expected in cell
          double mu = selphi*nc; // just Poisson selphi
          uint_8 npoiss = rg_.PoissonAhrens(mu); // Poisson fluctuate
          randomcat_(i,j,k) = (double)npoiss;
                
          wrgals_(i,j,k)=weights_(i,j,k)*randomcat_(i,j,k); // deleted alpha
        }
        else {
          randomcat_(i,j,k) = 0;
          weights_(i,j,k) = 0;
          wrgals_(i,j,k) = 0;
        }
              
        if (SaveArr)
          zc_(i,j,k) = redshift;
              
      }   // end of loop over k (z-direction)
      ccnt++;  pgb.update(ccnt); 
    } // end of loop over j (y-direction)
  } // end of loop over  (x-direction)
  nrand_=randomcat_.Sum(); 
  wnrand_=wrgals_.Sum();
  cout <<"    number of gals in grid = "<< ng_;
  cout <<", number of gals in random grid = "<< nrand_ <<endl;
  cout <<"    Number of galaxies in RANDOM weighted grid = "<< wnrand_;
  cout <<", should be approximately input density*npixels = "<< nc*n_obs_pixels_ <<endl;
  cout <<"    Calculated weighted random grid density = "<< wnrand_/n_obs_pixels_;
  cout <<", input density = "<< nc <<" (should be approximately the same)"<<endl;
        
  // Ratio of number of observed WEIGHTED gals to WEIGHTED gals in random catalog
  if (ngw_>0)
    alph_=(double)(ngw_/wnrand_); 
  else
    alph_=(double)(ng_/wnrand_); 
  cout <<"    alpha=ngw/nsw="<< alph_ <<endl; // alpha should be ntrue/nsyn
  VarianceRandomGrid();
  cout <<"    Mean of weighted random grid should be roughly "<< nc;
  cout <<" if all pixels are observed, "<< nc*n_obs_pixels_/Npix_ <<" if not"<<endl;
  double minv,maxv;
  randomcat_.MinMax(minv,maxv);
  cout <<"    Minimum value in random cat grid = "<< minv <<endl;
  cout <<"    Maximum value in random cat grid = "<< maxv <<endl;
  weights_.MinMax(minv,maxv);
  cout <<"    Minimum value in weights grid = "<< minv <<endl;
  cout <<"    Maximum value in weights grid = "<< maxv <<endl;

  // normalise the array
  NormRArray();

  // write the arrays here!
  fos_ << wrgals_;
  if (SaveArr)
    fos_ << zc_;

  tm.Split();
  cout <<"    Elapsed time "<< tm.TotalElapsedTime()<<endl;
  cout <<"    EXIT Cat2Grid::RandomGrid()"<<endl<<endl; 

};

void Cat2Grid::NormNArrays(){           
  double normm=(double)n_obs_pixels_/ng_;
  ngals_ *= normm;
  double normw = (double)n_obs_pixels_/ngw_;
  wngals_*= normw;
  cout <<"     Normalising n-gals array by "<< normm <<endl;
  cout <<"     Normalising weighted n-gals array by "<< normw <<endl;
}

void Cat2Grid::NormRArray(){            
  double normr = (double)(n_obs_pixels_/wnrand_)*( 1/sqrt(alph_));
  cout <<"     Normalising weighted random array by "<< normr <<endl;
  wrgals_*= normr; 
  VarianceRandomGrid();
}

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
  
  nz = int( (z_max-z_min)/binz +0.5)+1;
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
  if (i_zcat_max-i_zcat_min > 12 ) {
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



vector<double> Cat2Grid::XtractSubArray(TArray<r_8>& nsub, sa_size_t dNp, 
                                        sa_size_t Np, double theta, int arrayflag)
// Np    : the starting pixel number of sub-array
// dNp   : width of sub-array in z-direction
// theta : MAX angle grid should cover, if want whole grid just make sure theta
//                 is any value greater than angular extent of whole grid
// arrayflag : specifies which array: ngals,wngals,wrgals to return in nsub
{
  cout <<"    Cat2Grid::XtractSubArray() Z-axis split version"<<endl;
        
  cout <<"    Sub-array beginning at array pixel index "<< Np <<endl;
  sa_size_t Nx = GetNTransPix(Np,theta);
        
  // min and max pixel indices
  sa_size_t minxi = idx_-Nx;
  sa_size_t minyi = idy_-Nx;
  sa_size_t maxxi = idx_+Nx;
  sa_size_t maxyi = idy_+Nx;
  sa_size_t minzi = Np;
  sa_size_t maxzi = Np+dNp-1;
        
  // check min,max pixel indices
  if (maxxi>ngals_.SizeX()-1) {
    cout <<"    maxx="<< maxxi <<", so setting maxx="<< ngals_.SizeX()-1 <<endl;
    maxxi = ngals_.SizeX()-1;
  }
  if (maxyi>ngals_.SizeY()-1) {
    cout <<"    maxy="<< maxyi <<", so setting maxy="<< ngals_.SizeY()-1 <<endl;
    maxyi = ngals_.SizeY()-1;
  }
  if (maxzi>ngals_.SizeZ()-1) {
    cout <<"    maxz="<< maxzi <<", so setting maxz="<< ngals_.SizeZ()-1 <<endl;
    maxzi = ngals_.SizeZ()-1;
  }
  if (minxi<0) {
    cout <<"    minx="<< minxi <<", so setting minx=0"<<endl;
    minxi = 0;
  }
  if (minyi<0) {
    cout <<"    miny="<< minyi <<", so setting miny=0"<<endl;
    minyi=0;
  }
  if (minzi<0) {
    cout <<"    minz="<< minzi <<", so setting minz=0"<<endl;
    minzi=0;
  }
                
  if (minxi<0||minyi<0||minzi<0)
    throw ParmError("ERROR! either minx<0||miny<0||minz<0");
  if (maxxi>ngals_.SizeX()-1||maxyi>ngals_.SizeY()-1||maxzi>ngals_.SizeZ()-1)
    throw ParmError("ERROR! either maxx>ngals_.SizeX()-1||maxy>ngals_.SizeY()-1||maxz>ngals_.SizeZ()-1");

  cout <<"    Ranges of sub array: (xmin,xmax) = ("<< minxi <<","<< maxxi <<"),";
  cout <<" (ymin,ymax) = ("<< minyi <<","<< maxyi <<"), (zmin,zmax) = (";
  cout << minzi <<","<< maxzi <<")"<<endl;
        
  if(arrayflag<0)
    nsub = ngals_(Range(minxi,maxxi), Range(minyi,maxyi), Range(minzi,maxzi)).PackElements();
  else if (arrayflag==0)
    nsub = wngals_(Range(minxi,maxxi), Range(minyi,maxyi), Range(minzi,maxzi)).PackElements();
  else if (arrayflag>0)
    nsub = wrgals_(Range(minxi,maxxi), Range(minyi,maxyi), Range(minzi,maxzi)).PackElements();
                
  // need to return min and max redshifts of this sub grid
  // min redshift is redshift corresponding to: Zmin_+Np*cellsizeZ_;
  // max redshift is redshift corresponding to:  = Nx*cellsizeZ_/sin(theta)
  double mind = Zmin_+Np*cellsizeZ_;
  double maxd = mind + (maxzi-minzi)*cellsizeZ_;//(Nx*cellsize_)/sin(theta);
  double centerd = (mind + maxd)/2;
        
  // We have to initialize the distance to redshift conversion interpolator
  int_8 nz=1000;
  vector<double> zrs, codist;
  double minz=0, maxz=10;
  double dz = (maxz-minz)/(nz-1);
  for(int kk=0; kk<nz; kk++) 
    {
      double zs=minz+kk*dz;
      su_.SetEmissionRedShift(zs);
      double cod =su_.RadialCoordinateMpc(); // radial distance 
      zrs.push_back(zs);
      codist.push_back(cod); 
    }
  SInterp1D dist2z(codist,zrs,codist[0],codist[codist.size()-1],2*nz);
        
  double minza = dist2z(mind);
  double maxza = dist2z(maxd);
  double cenza = dist2z(centerd);
  vector<double> zedge;
  zedge.push_back(minza);
  zedge.push_back(cenza);
  zedge.push_back(maxza);
  cout << "    min,center,max redshifts of array = "<< minza <<","<< cenza;
  cout <<","<< maxza <<endl;
        
  cout <<"    EXIT Cat2Grid::XtractSubArray()"<<endl<<endl;

  return zedge;
};


vector<double> Cat2Grid::XtractSubArray(TArray<r_8>& nsub, long x1, long x2,
                                        long y1, long y2, long z1, long z2, int arrayflag)
// Extract sub array uses pixel ranges in each dimension
{

  cout <<"    Cat2Grid::XtractSubArray() Sub-array range version"<<endl;
        
  if(arrayflag<0)
    nsub = ngals_(Range(x1,x2), Range(y1,y2), Range(z1,z2)).PackElements();
  else if (arrayflag==0)
    nsub = wngals_(Range(x1,x2), Range(y1,y2), Range(z1,z2)).PackElements();
  else if (arrayflag>0)
    nsub = wrgals_(Range(x1,x2), Range(y1,y2), Range(z1,z2)).PackElements();
        
  // need to return min and max redshifts of this sub grid
  // min redshift is redshift corresponds to distance to first pixel in z-dim
  // max redshift is redshift corresponds to distance to last pixel in z-dim
  double mind = Zmin_+z1*cellsizeZ_;
  double maxd = mind + (z2-z1)*cellsizeZ_;
  double centerd = (mind + maxd)/2;
        
  // We have to initialize the distance to redshift conversion interpolator
  int_8 nz=1000;
  vector<double> zrs, codist;
  double minz=0, maxz=10;
  double dz = (maxz-minz)/(nz-1);
  for(int kk=0; kk<nz; kk++) 
    {
      double zs=minz+kk*dz;
      su_.SetEmissionRedShift(zs);
      double cod =su_.RadialCoordinateMpc(); // radial distance 
      zrs.push_back(zs);
      codist.push_back(cod); 
    }
  SInterp1D dist2z(codist,zrs,codist[0],codist[codist.size()-1],2*nz);
        
  double minza = dist2z(mind);
  double maxza = dist2z(maxd);
  double cenza = dist2z(centerd);
  vector<double> zedge;
  zedge.push_back(minza);
  zedge.push_back(cenza);
  zedge.push_back(maxza);
  cout << "    min,center,max redshifts of array = "<< minza <<","<< cenza;
  cout <<","<< maxza <<endl;
        
  cout <<"    EXIT Cat2Grid::XtractSubArray()"<<endl<<endl;

  return zedge;
};


sa_size_t Cat2Grid::GetNTransPix(sa_size_t Np, double theta)
// return number of pixels in x or y dimension that cover
// an angle of size theta at distance Dp
{

  double Dp = Zmin_ + Np*cellsizeXY_; // min distance to center pixel on x-y plane
  double xw = Dp*tan(theta)-cellsizeXY_/2;
  sa_size_t Nx=(sa_size_t)floor(xw/cellsizeXY_);
        
  cout <<"    Distance to sub array = "<< Dp <<" Mpc"<<endl;
  cout <<"    Number of pixels which cover "<< theta <<" radians";
  cout <<" at distance "<< Dp <<" Mpc, = "<< Nx <<endl;
        
  return Nx;
};



void  Cat2Grid::WriteHeader(string IncatName) {
  // Commented by Cecile: we can write the header in any case
  // if ( ngo_<1||nrand_<1 )
  //          throw ParmError("ERROR! Not ready to write Fits Header yet!"); 

  fos_.WriteKey("NX",Nx_," axe transverse 1");
  fos_.WriteKey("NY",Ny_," axe transverse 2");
  fos_.WriteKey("NZ",Nz_," axe longitudinal (redshift)");
  fos_.WriteKey("DX",cellsizeXY_," Mpc");
  fos_.WriteKey("DY",cellsizeXY_," Mpc");
  fos_.WriteKey("DZ",cellsizeZ_," Mpc");
  fos_.WriteKey("ZREF",zref_," reference redshift");  
  fos_.WriteKey("ThRAD",SkyArea_," radius of circular sky area");
  fos_.WriteKey("MeanOverDensity",mean_overdensity_," mean dens SimLSS delta-fudge");
  fos_.WriteKey("InCat",IncatName," original catalog");
  fos_.WriteKey("NGALCAT",ngall_," N gals in InCat");
  fos_.WriteKey("NOBSPIX",n_obs_pixels_," original catalog");
  fos_.WriteKey("NGRID",ng_," N gals in grid");
  fos_.WriteKey("NOGRID",ngout_," N gals outside grid");
  fos_.WriteKey("NWGRID",ngw_," N weighted gals in grid");    
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


void Cat2Grid::ZeroGalArrays(){
  wngals_.ZeroSize();
  ngals_.ZeroSize();
}

void Cat2Grid::WriteGalArrays() {
  ngals_.PackElements();
  wngals_.PackElements();
  fos_ << ngals_; fos_ << wngals_;
}

void Cat2Grid::SaveNgArray(string outfileroot, string IncatName) {
  
  string outfile;
  outfile = "tmp_" + outfileroot;//outfileroot +"_ngals.fits";
  cout <<"    Writing NGALS array to "<< outfile <<endl;
  FitsInOutFile fos(outfile, FitsInOutFile::Fits_Create);
  
  fos << ngals_;
  //WriteHeader(fos,IncatName);
  fos.WriteKey("NGRID",ng_," N gals in grid");
  fos.WriteKey("NOGRID",ngout_," N gals outside grid");
}

        
