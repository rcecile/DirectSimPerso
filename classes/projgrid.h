/* ----
   Project   LSST/BAO/PhotoZ
   C.Renault, R.Ansari , Jan 2017
                                                     -------  */

#ifndef PROJGRID_SEEN
#define PROJGRID_SEEN

#include <iostream>
#include <fstream>
#include <string>

#include <vector>

//----- sophya includes 
#include "machdefs.h"
#include "sopnamsp.h"
#include "rotation3d.h"

/*
 #include <slininterp.h>
using namespace SOPHYA;
using namespace std;
*/

/* Quelques classes utilitaires - pour transmettre les arguments, etc ... */
class GridCenter {
public:
  GridCenter(double rc=3000., double t0=0., double p0=0.)
  { r_center=rc; theta0=t0;  phi0=p0;  }
  GridCenter(GridCenter const &a)
  { r_center=a.r_center; theta0=a.theta0;  phi0=a.phi0; }
  double r_center;
  double theta0, phi0;
};


class GridDef {
public:
  GridDef()
  { Nx=Ny=Nz=256; dx=dy=dz=6.; r_center=3000.; theta0=phi0=0.; }
  GridDef(long int nx, long int ny, long int nz, double ddx, double ddy, double ddz, double rc, double t0=0., double p0=0.)
  { Set(nx,ny,nz,ddx,ddy,ddz,rc,t0,p0); }
  void Set(long int nx, long int ny, long int nz, double ddx, double ddy, double ddz, double rc, double t0=0., double p0=0.)
  { Nx=nx; Ny=ny; Nz=nz; dx=ddx; dy=ddy; dz=ddz; r_center=rc; theta0=t0;  phi0=p0;  return; }
  GridDef(GridDef const& a)
  { Nx=a.Nx; Ny=a.Ny; Nz=a.Nz; dx=a.dx; dy=a.dy; dz=a.dz;  r_center=a.r_center; theta0=a.theta0;  phi0=a.phi0;  }
  GridDef& operator = (GridDef const& a)
  { Nx=a.Nx; Ny=a.Ny; Nz=a.Nz; dx=a.dx; dy=a.dy; dz=a.dz;  r_center=a.r_center; theta0=a.theta0;  phi0=a.phi0; return *this; }
  void SetCenter(GridCenter const &a)
  { r_center=a.r_center; theta0=a.theta0;  phi0=a.phi0;  }
  
  double GetCellVolume() const
  { return (dx*dy*dz); }
  double GetVolume() const
  { return (double)(Nx*Ny*Nz)*(dx*dy*dz); }
  
  long int Nx,Ny,Nz;
  double dx, dy, dz;
  double r_center;
  double theta0, phi0;
};

// on definit le type de contenu des grilles, r_8 ( ou r_4) 
 
#define TFG r_8
// -- ProjGrid class: projecting into a grid, defined with a rotation with respect to the catalog Oxyz reference system
class ProjGrid {
public:
// Constructor
 ProjGrid( GridDef& outg) :
  outg_(outg), out_grid(outg_.Nx, outg_.Ny, outg_.Nz),
    out_rot((fabs(outg_.theta0)>1.e-9)?outg_.phi0+Angle::PioTwoCst():0., outg_.theta0, 0.)  {
    // ---
  ocx_= (double)OutNx()*0.5;
  ocy_= (double)OutNy()*0.5;
  ocz_= (double)OutNz()*0.5;
  dx_ = (double)getOutdX();
  dy_ = (double)getOutdY();
  dz_ = (double)getOutdZ();
  Nx_ = (double)OutNx();
  Ny_ = (double)OutNy();
  Nz_ = (double)OutNz();
  }

  // Constructeur de copie
 ProjGrid(ProjGrid const& a) :
  outg_(a.outg_), out_grid(a.out_grid), out_rot(a.out_rot), 
    ocx_(a.ocx_), ocy_(a.ocy_), ocz_(a.ocz_), dx_(a.dx_), dy_(a.dy_), dz_(a.dz_) {
  } 
  // destructeur 
  ~ProjGrid() { }

  ProjGrid& operator=(ProjGrid const& a)
    {
      outg_=a.outg_; out_grid=a.out_grid; out_rot=a.out_rot; 
      ocx_=a.ocx_; ocx_=a.ocy_; ocz_=a.ocz_;
      dx_=a.dx_; dx_=a.dy_; dz_=a.dz_;
   }
    // Return the number of cells alons each direction  
  inline sa_size_t OutNx()   const { return out_grid.SizeX(); }
  inline sa_size_t OutNy()   const { return out_grid.SizeY(); }
  inline sa_size_t OutNz()   const { return out_grid.SizeZ(); }

  // return the cell size along each direction
  inline double getOutdX() const {  return outg_.dx; } 
  inline double getOutdY() const {  return outg_.dy; } 
  inline double getOutdZ() const {  return outg_.dz; } 

  inline bool AddGal(double x, double y, double z, double invselfunc) 
  {
    sa_size_t ix, jy, kz;
    // CECILE aucune gal ne passe le cut, plein de x,y, negatif : y a un probleme qui reste 
    // parfois ne compile que si make clean ...
    Vector3d vout=out_rot.Rotate( Vector3d(x,y,z) );
    ix=(sa_size_t)floor(vout.X()/dx_+ocx_);
    jy=(sa_size_t)floor(vout.Y()/dy_+ocy_);
    kz=(sa_size_t)floor((vout.Z()-outg_.r_center)/dz_+ocz_);
    cout << ix << " " << jy << " " << kz << " et "  << out_grid.SizeX() << " " << out_grid.SizeY() << " " << out_grid.SizeZ()<< endl ;
    if ((ix<0) || (ix>=out_grid.SizeX()) || (jy<0) || (jy>=out_grid.SizeY()) || (kz<0) || (kz>=out_grid.SizeZ()))  return false;
    out_grid(ix,jy,kz) += invselfunc;
    return true;
  }

  // Acces a la grille remplie 
  TArray< TFG > & getGrid() { return out_grid; }

  protected:
  GridDef& outg_;              // input grid definition 
  TArray< TFG > out_grid;       // output grid

  EulerRotation3D out_rot;
  double ocx_, ocy_, ocz_, dx_, dy_, dz_;  // normalized coordinates (cell size=1) of the central cell for the output grid 
  double Nx_, Ny_, Nz_;  // normalized coordinates (cell size=1) of the central cell for the output grid 
  //  sa_size_t ocx_, ocy_, ocz_, dx_, dy_, dz_, Nx_, Ny_, Nz_;  // normalized coordinates (cell size=1) of the central cell for the output grid 
};




#endif
