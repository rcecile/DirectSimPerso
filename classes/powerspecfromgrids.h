/**
 * @file  powerspecfromgrids.h
 * @brief Compute power spectrum from an array of over-densities
 *
 *
 * @author Alex Abate and Reza Ansari
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2008
 * @date 2008
 *
 */

#ifndef POWERSPECFROMGRIDS_H_SEEN
#define POWERSPECFROMGRIDS_H_SEEN


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
#include "fabtwriter.h"
#include "array.h"
#include "hisprof.h"
#include "histerr.h"
#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "stsrand.h"

// DirectSim
#include "sinterp.h"
#include "constcosmo.h"

/** @class PowerSpecFromGrids
  *
  * Reads in galaxy fluctuation grid and computes power spectrum 
  */
class PowerSpecFromGrids
{
public:
 	
	/** Constructor 
        @param delta_gal      galaxy fluctuation grid
        @param dx             grid pixel size in x-dimension (Mpc)
        @param ratio_AngDiam  ratio angular diameter with different cosmologies)              */
  PowerSpecFromGrids (TArray<r_8> delta_gal, double dr, double ratio_AngDiam); 


	/** Compute Fourier transform of galaxy fluctuation grid
	    @param real_grid      real value array 
	    @param fourier_grid   Fourier space array                             */
	void ComputeFourier(TArray<r_8>& real_grid, TArray< complex< r_8 > >& fourier_grid);  
	
	/** Compute Fourier transform of Fourier space grid back to real space 
	    @param fourier_grid   Fourier space array
	    @param real_grid      real value array                                */
	void ComputeFourierBack(TArray< complex<r_8> >& fourier_grid, TArray<r_8>& real_grid);
	
	/** Compute spacing of Fourier transformed grid ie \f$ dk_x =  2\pi/(N_xdx) \f$*/
	// ratio_AngDiam = AngularDiameter(w0,wa in parameters) / AngularDiameter(w0, wa in grid fits header)
	void SetDKDR(); 
	
	/** Compute power spectrum from Fourier transformed galaxy fluctuation grid.
	    and fill histogram. Power spectrum is a one-dimensional function of
	    wavevector length (universe is isotropic)
	    @param hp        histogram of power spectrum values (to be filled)  */
	    //modif Adeline : add  hnode_ and hkeepMode_ 
	double AccumulatePowerSpectra(HProf& hp); 	
	
	/** Write power spectra to a file. None of the power spectra read in as
	    arguments are normalised yet. Columns written to the file are (1) k-values, 
	    (2) estimated P(k), (3) simulation P(k)/ simulation with distortion P(k),
	    (4) shot-noise grid PS (5) error on PS
	    @param fname       file to write power spectra to
	    @param Pdata       raw estimated power spectrum
	    @param volData     volume of catalog used to estimate power spectrum
	    @param Pdata_noise shot noise	
	    @param nGalGrid    number of galaxies in weight grid	
	**/
	void WritePS(string fname, HProf& Pdata, r_4 Voldata, HProf& Pdata_noise, double h, double nGalGrid=0,int nGridData=1, bool comp_SN=false); 
	
	// Minor functions
	//void SetMaxKrad(double maxk){maxk_=maxk;};
	
	/** Return maximum wavevector in the x-dimension                          */
	double ReturnKxMax(){kxmax=four_.SizeX()*Dkx_; return kxmax;};
	
	/** Return maximum wavevector in the y-dimension                          */
	double ReturnKyMax(){kymax=four_.SizeY()*Dky_/2; return kymax;};
	
	/** Return maximum wavevector in the z-dimension                          */
	double ReturnKzMax(){kzmax=four_.SizeZ()*Dkz_/2; return kzmax;};
	
	/** Return galaxy fluctuation grid                                        */
	inline TArray<r_8>& ReturnRealDist() { return drho_; };
	
	/** Return Fourier transformed galaxy fluctuation grid                    */
	inline TArray< complex< r_8 > >& ReturnFour() { return four_; };
	
	/** Return galaxy fluctuation grid convolved with 1D Gaussian             */
	void ReturnConvDist(TArray<r_8>& drhoconv) { drhoconv = drhoconv_; };
	
	/** Return estimated power spectrum                                       */
	TVector<r_8> ReturnPk() { return Pobs_; };
	
	/** Return k values of estimated power spectrum                           */
	TVector<r_8> Returnk() { return kobs_; };
	
	/** Zero the size of the galaxy fluctuation grid and FT of galaxy 
	    fluctuation grid                                                      */
	void ZeroSizeArrays() { drho_.ZeroSize(); four_.ZeroSize(); };
	
	/** Zero the size of the FT of galaxy fluctuation grid                    */
	void ZeroFourArray() { four_.ZeroSize(); };
	
	/** Set redshift of power spectrum */
	void Setzc(double zc) { zc_=zc; };

protected:
	TArray<r_8> drho_;	            /**< galaxy fluctuation grid              */
	TArray<r_8> drhoconv_;	        /**< galaxy fluctuation grid convolved with 1D Gaussian  */
	TArray< complex< r_8 > > four_; /**< FT of galaxy fluctuation grid                       */
	TArray< complex< r_8 > > fourconv_; /**< FT of galaxy fluc. grid convolved with Gaussian */
	TVector<r_8> Pobs_;             /**< estimated power spectrum             */
	TVector<r_8> kobs_;             /**< k values of estimated power spectrum */
	r_4 Vol_;			            /**< volume of survey                     */
	sa_size_t Nx_;                  /**< number of grid pixels in x-dimension */
	sa_size_t Ny_;                  /**< number of grid pixels in y-dimension */
	sa_size_t Nz_;                  /**< number of grid pixels in z-dimension */
	int_8 NRtot_;		            /**< total number of pixels               */
	long NCx_;			            /**< Nz_/2+1                              */
	double Dx_;                     /**< size of grid pixel in x-dimension    */
	double Dy_;                     /**< size of grid pixel in y-dimension    */
	double Dz_;	                    /**< size of grid pixel in z-dimension    */
	double ratio_AngDiam_; /** ratio between angular diameters for different cosmologies**/
	double Dkx_;         /**< size of grid pixel in Fourier space x-dimension */
	double Dky_;         /**< size of grid pixel in Fourier space y-dimension */
	double Dkz_;         /**< size of grid pixel in Fourier space z-dimension */
	double kxmax;        /**< (2PI)/(Nx_*Dx_)                                 */
	double kymax;        /**< (2PI)/(Ny_*Dy_)                                 */
	double kzmax;        /**< (2PI)/(Nz_*Dz_)                                 */
	bool erronredshift_; /**< if true error is on redshit, if false error is on Z-axis                 */
	double zc_;          /**< redshift of power spectrum                      */
};

#endif
