/**
 * @file  cat2grid.h
 * @brief grid a galaxy catalog for power spectrum analysis
 *
 * @todo move some methods out to a selection function class eg SaveSelecFunc
 *
 * @author Alex Abate, Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2009
 * @date 2009
 *
 */
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <typeinfo>

// sophya
#include "machdefs.h"
#include "sopnamsp.h"
#include "array.h"
#include "fabtwriter.h"
#include "hisprof.h"
#include "histerr.h"
#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "stsrand.h"
#include "mydefrg.h"
#include "swfitsdtable.h"
#include "fitsmanager.h"
#include "resusage.h"
#include "timing.h"

#include "cosmocalcs.h"
#include "constcosmo.h"
#include "slininterp.h"
#include "schechter.h"
#include "selectfunc.h"
#include "mass2gal.h"
// add new class to loop over grids with rotation
#include "projgrid.h"

//Root
#include "TFile.h"
#include "TF1.h"
#include "TRandom3.h"
using namespace std;

/** @class GalRecord

GalRecord class holds galaxy variables

*/
class GalRecord {
 public:

  /** Constructor */
  GalRecord() { alpha=delta=glong=glat=xcoo=ycoo=zcoo=0.; zs=zo=0.; type=1; u=g=r=i=z=y=24.; MB=-22; }

  /** Print galaxy properties */
  void Print()
  { cout <<"Printing GalRecord: ";
    cout <<"alpha="<< alpha <<", delta="<< delta <<", x="<< xcoo <<", y="<< ycoo <<", zs="<< zs <<", zo="<< zo;}

  double alpha;    /**< right ascension coodinate                           */
  double delta;    /**< declination coodinate                               */
  double glong;    /**< galactic longitude coordinate                       */
  double glat;     /**< galactic latitude coordinate                        */
  double xcoo;     /**< Euclidean x coordinate                              */
  double ycoo;     /**< Euclidean y coordinate                              */
  double zcoo;     /**< Euclidean z coordinate                              */
  double zs;       /**< spectroscopic redshift                              */
  double zo;       /**< observed redshift                                   */
  int type;        /**< galaxy type                                         */
  double u;        /**< u magnitude                                         */
  double g;        /**< g magnitude                                         */
  double r;        /**< r magnitude                                         */
  double i;        /**< i magnitude                                         */
  double z;        /**< z magnitude                                         */
  double y;        /**< y magnitude                                         */
  double MB;     /**< absolute magnitude                                  */
};


/** @class Cat2Grid
 *
 * Reads in a galaxy catalog (ra,dec,z) from a SwFitsDataTable which must have column names:
 * "phi" and "theta", but the z column name can be set via an argument
 *
 * Lays a grid across the galaxy catalog with grid cell size (default 8 Mpc) and 
 * grid position specified in SetGrid function 
 * 
 * OPTIONALLY: adds statistical photo-z error sig = PZerr*(1+z) to z-dimension (if PZerr>0)
 * OPTIONALLY: corrects for selection effects (if SetSelectionFunction is called, which sets sfcompute_=true)
 */
class Cat2Grid
{
 public:

  /** Constructor - Normal usage
      @param dt       data table containing galaxy catalog
      @param su       cosmology calculator
      @param rg       random number generator
      @param fos      FITS file containing gridded galaxy data
      @param ZOCol    name of column in galaxy catalog containing observed z
      @param ZSCol    name of column in galaxy catalog containing spec z
      @param RadialZ  if true sets the z-dimension to be radial direction 
      @param PZerrAxis    size of Gaussian photo-z error on Z axis
      @param PzerrReds    size of Gaussian photo-z error on redshift
      @param Print    if true prints extra to the screen                    */
  Cat2Grid(SwFitsDataTable& dt, SimpleUniverse& su, RandomGenerator& rg,
           FitsInOutFile& fos, string ZOCol="zp", string ZSCol="z",
           bool RadialZ=false, double PZerr=0,bool Print=true, string ObsCat="", double ratio_cell=1.);
                 
  /** Constructor - Takes the interp z->d functions as arguments rather than
      calculating them every time so you can use the Row2Record, 
      Rec2EuclidCoord functions
      @param dt       data table containing galaxy catalog
      @param su       cosmology calculator
      @param rg       random number generator
      @param dist2z   distance to redshift conversion table
      @param z2dist   redshift to distance conversion table
      @param ZOCol    name of column in data table containing observed z
      @param ZSCol    name of column in data table containgin spec z
      @param RadialZ  if true sets the z-dimension to be radial direction   */
  Cat2Grid(SwFitsDataTable& dt, SimpleUniverse& su, RandomGenerator& rg,
           SInterp1D dist2z, SInterp1D z2dist, 
           string ZOCol="zp", string ZSCol="zs", bool RadialZ=false);

  /** Copy constructor */
  Cat2Grid(Cat2Grid const& a);
        
  /** Destructor */
  virtual ~Cat2Grid() {};
        
  /** Copy method */
  virtual Cat2Grid& Set(const Cat2Grid& a);
        
        
  /* MAIN FUNCTIONS */
        
  // computes equivalent distance error to the photoz error (PZerr) given to the constructor
  //double CompDistErr(); 
        
  /** Find the minimum and maximum x,y,z,zs (zs=phot or spec) of galaxies in 
      catalog                                                               */
  double FindMinMaxCoords(); 
                
  /** Set grid over the galaxy data by specifying redshift of central pixel, 
      size of the pixels, and (approx) number of pixels along each dimension
      @param Nx    (approx) number of pixels along x-dimension
      @param Ny    (approx) number of pixels along y-dimension
      @param Nz    (approx) number of pixels along z-dimension
      @param R     pixel size in Mpc
      @param zref  redshift of central pixel                                */
  void SetGrid(int_8 Nx, int_8 Ny, int_8 Nz, double R, double zref); 

  /** Define a set of grids to which galaxies would be projected */
  inline void SetGrids(vector<ProjGrid>& vgrids) { 
    vgrids_=vgrids; 
  }   

  /** Create a histogram of observed redshifts divided by true redshifts 
      and save it to a text file 
      @param SFTextFile    selection function is saved to file called [SFTextFile]_nofz.txt 
      @param FullCat       name of FITS file containing z of all simulated gals
      if >1 file to read files are named FullCat_#ofn.fits
      @param ZCol          name of column in FullCat containing true z
      @param Nfiles        number of files n to read in  */
  void SaveSelecFunc(string SFTextFile, string FullCat, string ObsCat, string ZFCol="zs", string  ZSCol="zs",  string ZOCol="zp", bool MakeFullHist=true);
        
  /** Set selection function                                                */
  inline void SetSelectionFunction(SelectionFunctionInterface& sf) { 
    selfuncp_=&sf; sfcompute_=true; 
    cout <<"    Set selection function"<<endl;};
        
  /** Set bias function                                                */
  inline void SetBiasFunction(SelectionFunctionInterface& b) { 
    biasp_=&b; 
    cout <<"    Set bias function"<<endl;};
        
  /** Set sigma function                                                */
  inline void SetSigrFunction(SelectionFunctionInterface& b) { 
    sigrp_=&b; 
    cout <<"    Set sigma for random function"<<endl;};
            
  /** Project the galaxy distribution into the grid, fill grid arrays
      @param SkyArea    angle covered by observation cone                   */
  void GalGrid(double SkyArea = 999);
          
  /** Compute <ng(z)> form the histogram of the observed number of galaxies  */
  void computeMeanNg_z(double NormNgalMean = 1.);

  /** Check if all pixels are seen. Error if not                         */
  void ObsPixels();
        
  /** Convert row from catalog data table into a GalRecord                  */
  void Row2Record(DataTableRow& rowin, GalRecord& rec);
        
  /** Return true if galaxy accepted to be used in analysis. All galaxies are
      currently accepted                                                    */
  bool Filter(GalRecord& rec); 
        
  /** Convert galaxy position in GalRecord (phi,theta,z) into Shell 
      coordinates
      @param rec        galaxy record
      @param x          Shell x coordinate
      @param y          Shell y coordinate
      @param z          Shell z coordinate
      @param redshift   redshift of galaxy                                  */
  void Rec2ShellCoord(GalRecord& rec, double& x, double& y, double& z, double& redshift);
        
  /** Convert galaxy position in GalRecord (phi,theta,z) into Euclidean 
      coordinates
      @param rec        galaxy record
      @param x          Euclidean x coordinate
      @param y          Euclidean y coordinate
      @param z          Euclidean z coordinate
      @param redshift   redshift of galaxy                                  */
  void Rec2EuclidCoord(GalRecord& rec, double& x, double& y, double& z, double& redshift);
        
  /** Add galaxy with position x,y,z and weighting 1/phi to correct pixel in
      galaxy grid
      @param x          galaxy Euclidean x coordinate
      @param y          galaxy Euclidean y coordinate
      @param z          galaxy Euclidean z coordinate
      @param phi        weight of galaxy                                    */
  void AddToCell(double x, double y, double z,double phi=1);
        
  /** Return cartesian coord (x,y,z) of pixel cell given pixel cell index (i,j,k) 
      @param i    pixel index in x-dimension
      @param j    pixel index in y-dimension
      @param k    pixel index in z-dimension
      @param x    pixel coordinate in x-dimension
      @param y    pixel coordinate in y-dimension
      @param z    pixel coordinate in z-dimension                           */

   void AddToCellMultiGrids(double x, double y, double z,double phi=1);
        
  /** Return cartesian coord (x,y,z) of pixel cell given pixel cell index (i,j,k) 
      @param i    pixel index in x-dimension
      @param j    pixel index in y-dimension
      @param k    pixel index in z-dimension
      @param x    pixel coordinate in x-dimension
      @param y    pixel coordinate in y-dimension
      @param z    pixel coordinate in z-dimension                           */
 void GetCellCoord(sa_size_t i, sa_size_t j, sa_size_t k, double& x, double& y, double& z);
        
  /** Return cartesian coord of pixel cell in the x-dimenstion given a pixel 
      cell index 
      @param i    pixel index in x-dimension                                */
  inline double GetCellX(sa_size_t i)  
    { return Xmin_ + cellsize_/2 + i*cellsize_; };
        
  /** Return cartesian coord of pixel cell in the y-dimenstion given a pixel 
      cell index 
      @param j    pixel index in y-dimension                                */
  inline double GetCellY(sa_size_t j)  
    { return Ymin_ + cellsize_/2 + j*cellsize_; };
        
  /** Return cartesian coord of pixel cell in the z-dimenstion given a pixel 
      cell index 
      @param k    pixel index in z-dimension                                */
  inline double GetCellZ(sa_size_t k)  
    { return Zmin_ + cellsize_/2 + k*cellsize_; };
        
  /** Compute random weighted grid with same selection function as data 
      @param nc         mean density of random grid
      @param SaveArr    if true fill an array of pixel center redshifts     */
  void RandomGrid(double NormNgalMean=1, bool SaveArr=true, bool seed=false, bool SigRandom=false);

  /** Add Gaussian errors to redshifts  (Cecile)
      @param Err   redshift error size (Err*(1+z))                          */
  void SetGaussErrRedshift(double Err, double zref, bool seed);
                
              
  /** Add photo-z errors to redshifts compute from pdf (Adeline)
      @param Err   pdf File Name (root file)                          */        
  void SetPhotozErrRedshift(string pdfFileName, bool seed, double zcat_min, double zcat_max);

  /** Parameters over which pdf are computed (Adeline) */
  void InitParam();
        
  /** Read root file, creat an array with TF1  (Adeline)
      @param Name        */  
  void readRootFile(string rootFileName);

  /** Compute errors from photoz (Adeline)
   */
  double ComputeErrorFromPhotoz(GalRecord& rec);        

  /** Compute variance of random grid                                       */
  void VarianceRandomGrid();

  /** Write FITS header */
  void WriteHeader(string IncatName);
        
  /** Zero the size of the galaxy grids                                     */
  void ZeroGalArrays();
        
  /** Write the galaxy grids to a FITS file                                 */
  void WriteGalArrays();

       
  /* MINOR FUNCTIONS */
        
  /** Convert comoving distance to redshift                                 */
  double ConvCoD2z(double dc) { return dist2z_(dc); };
        
  /** Return the approx volume of the survey: (maxx-minx)*(maxy-miny)*(maxz-minz)*/
  r_4 ReturnSurveyVolume(){return Vol_;};
        
  /** Return the volume of the grid                                         */
  r_4 ReturnGridVolume(){return volgrid_;};
        
  /** Return number of galaxies in galaxy number grid                       */
  sa_size_t ReturnNg(sa_size_t i){return ng_[i];}; 
        
  /** Return number of galaxies in WEIGHTED galaxy number grid              */
  r_8 ReturnWg(sa_size_t i){return ngw_[i];}; 
        
  /** Return total number of galaxies in the simulation                     */
  sa_size_t ReturnNgAll(){return ngall_;};
                
  /** Return weighted random grid (may or may not be normalised)            */
  void ReturnWrgals(TArray<r_8>& wrgals) { wrgals = wrgals_; };
                
  /** Return redshift grid            */
  void ReturnZc(TArray<r_8>& zc) { zc = zc_; };
        
  /** Return alpha, number of galaxies in weighted galaxy number grid divided
      by weighted random number grid                                        */
  double ReturnAlpha(){ return alph_;};
        
  /** Return photometric redshift error in comoving distance units          */ 
  double ReturnPZDerr(){return PZDerr_;};
        
  /** Return grid specification                                             */
  TVector<r_8> ReturnGridSpec(){
    int row=4; TVector<r_8> gridspec(row); 
    gridspec(0)=Nx_; gridspec(1)=Ny_; 
    gridspec(2)=Nz_; gridspec(3)=cellsize_;
    return gridspec; };
                
  /** Set output filename root to output debugging files to                 */
  void SetDebugOutroot(string debugoutroot)
  { debugoutroot_=debugoutroot; DoDebug_=true; };
        
                
 protected:
        
        
  SwFitsDataTable& dt_;                   /**< data table containing galaxy catalog */
  SimpleUniverse& su_;                    /**< cosmology                            */
  RandomGenerator& rg_;                   /**< random number generator              */
  SelectionFunctionInterface* selfuncp_;    /**< selection function           */
  SelectionFunctionInterface* biasp_;      /**< bias function          */
  SelectionFunctionInterface* sigrp_;      /**< sigma for random function          */
  SInterp1D z2dist_;                /**< redshift to distance look up table   */
  SInterp1D dist2z_;                      /**< distance to redshift look up table   */
        
  bool DoDebug_;                  /**< true if debugging                              */
  bool AddGaussErr_;        /**< true if adding Gaussian error  */
  bool AddGaussErrReds_;    /**< true if adding Gaussian error to redshifts (Cecile)   */
  bool AddPhotoErrReds_;          /**< true if adding photo-z error on redshift (Adeline)   */
  bool ErrRandomSeed_ ;     /**< true if the seed for the error generation is new at each execution */
  bool sfcompute_;        /**< true if selection function has been set        */
  bool doBiasCorr_;        /** add a correction on the redshift corresponding to the bias */
  bool RadialZ_;                  /**< if true z-dimension IS radial direction        */

  FitsInOutFile& fos_;  /**< FITS file containing gridded galaxy data     */
  string debugoutroot_;/**< root file name to save things to when debugging */
        
  double PZerr_;      /**< size of photo-z error along z-axis               */
  double PZDerr_;           /**< size of photo-z error in comoving distance units */
  TFile *PDFrootFile_;/**< root file to read pdf for photo-z error computation (Adeline)                */
        
  string pdfFileName_ ;/**< name of rootfile containing the photo-z errors (Adeline)    */
  double z_min;     /**< parameters over which pdf ared computed (Adeline)                      */
  double z_max; 
  double binz;
  double type_min;
  double type_max;
  double binType;
  double mag_min;
  double mag_max;
  double binMag;
  int nz;
  int ntype;
  int nmag;
  double ***QualCutProb;      /**< Probability to pass the quality cut(Adeline)   */
  bool qualcut; // if no quality cut treatment qualcut=0 //
  TF1 ****photozErrorArray;      /**< Array of TF1 function, give the photo-z error (Adeline)           */
  TRandom3 rand;

  string ZOCol_;      /**< column name containing observed redshifts        */
  string ZSCol_;            /**< column name containing spec redshifts            */
  double volgrid_;    /**< volume of grid in Mpc^3                          */
  double cellsize_;   /**< pixel size of grid in Mpc                        */
  r_4 Vol_;           /**< approx volume of galaxy catalog                  */
  double SkyArea_;    /**< sky area: not clear on what this does?           */
  double xmin_;       /**< min x-coord of galaxies                          */
  double xmax_;       /**< max x-coord of galaxies                          */
  double ymin_;       /**< min y-coord of galaxies                          */
  double ymax_;       /**< max y-coord of galaxies                          */
  double zmin_;       /**< min z-coord of galaxies                          */
  double zmax_;       /**< max z-coord of galaxies                          */
  double Xmin_;       /**< min x-coord of grid                              */
  double Xmax_;       /**< max x-coord of grid                              */
  double Ymin_;       /**< min y-coord of grid                              */
  double Ymax_;       /**< max y-coord of grid                              */
  double Zmin_;       /**< min z-coord of grid                              */
  double Zmax_;       /**< max z-coord of grid                              */
  double zsmin_;      /**< min observed redshift of gals in catalog (could be spec or photo-z) */
  double zsmax_;                /**< max observed redshift of gals in catalog (could be spec or photo-z) */
  sa_size_t Nx_;      /**< number of pixels of grid in x-dimension          */    
  sa_size_t Ny_;      /**< number of pixels of grid in y-dimension          */ 
  sa_size_t Nz_;      /**< number of pixels of grid in z-dimension          */ 
  sa_size_t Npix_;    /**< total number of pixels in grid                   */ 
  double ratio_cell_;  /** cell size grid / cell size initial simulation    */
  sa_size_t ngo_;     /**< observed number of gals in sim                   */
  sa_size_t ngall_;   /**< total number of gals in sim                      */
  sa_size_t *ng_;      /**< number of gals inside grid               */
  sa_size_t *ngout_;   /**< number of gals outside grid              */ 
  r_8 *ngw_;           /**< number of gals inside weighted grid      */
  r_8 *meanz_;         /**< mean redshift of a grid      */
  r_8 wnrand_;        /**< number of gals in weighted random grid           */

  int_4 nz_0_;        /** for PDF proba, where to start the table filling    */
  int_8 idx_;         /**< index of center pixel in x-dimension             */
  int_8 idy_;         /**< index of center pixel in y-dimension             */
  int_8 idz_;         /**< index of center pixel in z-dimension         */              
  double DCref_;      /**< comoving distance to center pixel                */
  double zref_;     /**< redshift of center pixel                         */  
  double alph_;       /**< number of gals in weighted galaxy number grid/by weighted random number grid */
  sa_size_t Ic1_;     /**< index of phi or x column in data table           */
  sa_size_t Ic2_;     /**< index of theta or y column in data table         */
  sa_size_t Izs_;     /**< index of spec-z column in data table             */
  sa_size_t Iz_;      /**< index of obs-z column in data table              */
  sa_size_t Iid_;     /**< index of gal id column in data table             */
  sa_size_t Iug_;     /**< index of u-g color column in data table          */
  sa_size_t Igr_;     /**< index of g-r color column in data table          */
  sa_size_t Iri_;     /**< index of r-i color column in data table          */
  sa_size_t Iiz_;     /**< index of i-z color column in data table          */
  sa_size_t Izy_;     /**< index of z-y color column in data tabl           */    
  sa_size_t Itype_;   /**< index of type column in data tabl              */  
  sa_size_t Imag_;    /**< index of MA column in data tabl                */  
  double mean_overdensity_; /**< mean over-density of simlss grid AFTER setting cells with <-1 to =-1 */
  string ObsCat_;     /**< optionnal, list of catalog files, useful if more than 1 catalog must be read */
	  
  vector<ProjGrid> vgrids_;   /**< Vector of grids filled with weighted number of galaxies */
  
  TArray<r_8> wrgals_;  /**< weighted random number grid                   */
  TArray<r_8> zc_;          /**< redshift value at grid pixel centers          */
  Histo ngalz_;       /**< Histogram of the number   of galaxies per redshift bin */
  Histo ngalzw_;       /**< Histogram of the number   of galaxies per redshift bin weighted by inverse of selfunc*/
  SLinInterp1D ngalz_cell_; /** mean observed nb of galaxies per cell corrected cells set to 0 by argument*/
  SelectionFunctionInterface defselfunc_;   // Default selection function 
  SelectionFunctionInterface defbiasfunc_;   // Default bias function 
  SelectionFunctionInterface defsigrfunc_;   // Default sigma function 
  FitsInOutFile fosdefault_;
};
