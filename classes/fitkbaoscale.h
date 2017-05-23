/**
 * @file  fitkbaoscale.h
 * @brief Fit BAO scale by fitting decaying sinoid to the "wiggles only"
 *        power spectrum
 *
 * @todo implement more sophisticated fit method
 *
 * @author Alex Abate, Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2008
 * @date 2008
 *
 */
#ifndef FITBAO_H_SEEN
#define FITBAO_H_SEEN


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
#include "slininterp.h"

// DirectSim
#include "pkspectrum.h"
#include "cosmocalcs.h"
#include "chisqstats.h"
#include "constcosmo.h"


/** @class DecaySineFunc
  *
  * Computes Decaying sinusoid, See Blake & Glazebrook 2003 eqn 3
  *
  * The power of 1.4 in the decay term originates from the Silk damping fitting
  * formula in Eisenstein & Hu 1998.  Varying this decay length as well as the
  * amplitude will probably not have a significant effect on the fitted values of 
  * s (the BAO scale) 
  */
class DecaySineFunc
{
public:

    /** Constructor
        @param amp    amplitude of decaying sinusoid function 
        @param h      Hubble parameter in units of 100 km/s/Mpc
        @param decay  decay length                                            */
	DecaySineFunc(double h, double decay = 1.4) 
    	: h_(h) , decay_(decay) {};
    	
    /** Return decaying sinusoid function \f$ 1+Ak\exp(-(k/0.1h)^d)\sin(2\pi k/k_a) \f$
        @param k    
        @param s                                                             */
  	virtual double operator()(double k, double amp, double s) {
	  return ( 1+amp*k*exp(-  (pow(k/0.1/h_,decay_))  )*sin(s*k) );
	  // return ( 1+amp*sqrt(k)*exp(-  (pow(k/0.1,decay_))  )*sin(s*k) ); // suggestion Stephane mais pas convaincue (Cecile)
    		};  

    /** Print parameters of decaying sinusoid function                        */
   	void PrintParas()
		{ cout <<"    decay length = "<< decay_ <<", h = "<< h_ <<endl; };
	

protected:
    double h_;          	/**< Hubble parameter in units of 100 km/s/Mpc          */
    double decay_;      	/**< decay length                                       */
};


/** @class FitBAOScale
  *
  * Computes chi-square between wiggles-only power spectrum and decaying sinosoid 
  *
  */
class FitBAOScale
{
public:

	/** Constructor
	    @param pspectrum    observed power spectrum: columns(1,2,3) = (k,P(k),sigma_P(k))
	    @param su           cosmology 
	    @param ref_file     file with the ref spectrum without oscillation 
	    @param zref         redshift of power spectrum
	    @param simu_mode    if false, read shot noise and sigma */
  FitBAOScale(TArray<r_8> pspectrum, SimpleUniverse& su, double zref, string ref_file, bool simu_mode);

	//modif Adeline : psmooth is read in a file
  //	FitBAOScale(TArray<r_8> pspectrum, TArray<r_8> psmooth, SimpleUniverse& su, double zref, double sig8, double n);

	/** Initialise observed power spectrum variables and other vector variable 
	    sizes using pspectrum
	    @param pspectrum    observed power spectrum                           */
	void InitVect(TArray<r_8> pspectrum);
	
	/** Compute smooth (no wiggles) power spectrum
	    @param OmegaM    matter density
	    @param OmegaL    dark energy/cosmological constant density
	    @param OmegaB    baryon density
	    @param h         Hubble parameter in units of 100 km/s/Mpc
	    @param sig8      sigma8 (amplitude of fluctuations)
	    @param n         spectra index
	    @param R         scale of sigma8                                   
	void ComputeSmoothPS(double OmegaM, double OmegaL, double OmegaB, double h, 
	                                 double sig8=0.8, double n=1, double R = 8);
	not used anymore because sim spectrum read
	------------------------------------------
	   */
	/** Set grid of trial s values 
	    mins    start s value
	    maxs    end s value
	    ns      number of s values to try                                   */
	void Sets(double mins, double maxs, int ns, double mina, double maxa, int na)
		{ Smin_=mins; Smax_=maxs; nS_=ns;  Amin_=mina; Amax_=maxa; nA_=na; };
		
	/** Compute ratio between observed and reference power spectra            */
	void Pratio() {
	    for (int i=0; i<kobs_.Size(); i++){
		    Pratio_(i) = Pobs_(i)/Pref_(i); 
		    sig_(i)    = sig_(i) /Pref_(i); 
            }
	};
		    
	/** Compute the chi-square values as a function of s using ratio of power
	    spectra up to some maximum k
	    @param maxk    maximum k of power spectra to use in chi-square calc   */
	void ComputeChisq(double maxk);
	
	/** Return redshift of power spectrum                                     */
	double ReturnZPS() { return zref_; };
	
	/** Return observed power spectrum 
	    @param ko    k values
	    @param Po    observed P(k)
	    @param sig   errors on observed P(k)                                  */
	void ReturnPS(TVector<r_8>& ko, TVector<r_8>& Po, TVector<r_8>& sig)
		{ ko = kobs_; Po = Pobs_; sig = sig_; };
		
	/** Return ratio between observed and reference power spectra             */
	TVector<r_8> ReturnPSRatio() { return Pratio_; };
	
	/** Compute best-fit s with error with given n-sigma precision            */
	void BestfitStdDev(double& bestfit_S, double& bestfit_A, double& bestfit_Chi, double& siglow, double& sighigh, double& sigbest, int nsig = 1);
       
	/** Analytical approximation of the error on the power spectrum due to the 
	    sample variance 
	    @param VolCat    volume of galaxy catalog used to compute power spectrum */
	void CalcSigSampVar(double VolCat);
	
	/** Write results: redshift of power spectrum, best-fit s, n-sigma of 
	    error ranges, upper error range, lower error range to a file 
	    @param outfile    file to write results to                            */
	void WriteResults(string outfile);
	
	/** Write marginalized chi2 to a file 
	    @param outfile    file to write functions to   			  */
	void WriteChisq(string outfile);
	

protected:
	SimpleUniverse& su_;	 /**< cosmology                                   */
	TVector<r_8> kobs_;      /**< k values of observed power spectrum         */
	TVector<r_8> Pobs_;      /**< observed power spectrum                     */
	TVector<r_8> sig_;       /**< error on observed power spectrum            */
	TVector<r_8> Pref_;      /**< smooth reference power spectrum             */
	TVector<r_8> Pratio_;    /**< ratio between observed and reference power spectra 		  */
	TVector<r_8> Svals_list_;/**< grid of trial s values                     */
	TVector<r_8> Avals_list_;/**< grid of trial A values                      */
	TArray<r_8> Chisq_;     /**< chi-square value at each A,S value           */
	TVector<r_8> Chisq1D_;  /**< chi-square value at each S value, 1 row of A           */
	
	double Smin_;           /**< min s value in trial s grid               */   
	double Smax_;           /**< max s value in trial s grid               */ 
	int nS_;                /**< number of s values in trial s grid        */
	double Amin_;		 /**< min Amp value in trial Amp grid             */  
  	double Amax_; 		 /**< max Amp value in trial Amp grid             */  
   	int nA_;		 /**< number of Amp value in trial Amp grid       */  
   
	
	double zref_;            /**< redshift of power spectrum                  */
	double h_;               /**< Hubble parameter in units of 100 km/s/Mpc   */
	double bestfitS_;         /**< best fit s                                  */
	double bestfitA_;     /**< best fit amp (adeline)      		  */
	double bestChisq_;     /**< best fit amp (adeline)      		  */
	double errup_;           /**< upper error range                           */
	double errdown_;         /**< lower error range                           */
	int nsig_;               /**< n-sigma of error ranges                     */
	bool simu_mode_;         /**< if false, read shot noise and sigma         */
};

#endif


