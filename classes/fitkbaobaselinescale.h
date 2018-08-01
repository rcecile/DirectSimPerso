/**
 * @file  fitkbaobaselinescale.h
 * @brief Fit BAO scale by fitting decaying sinoid to the "wiggles only"
 *        power spectrum
 *
 * Created on: 2017
 * @date 2017
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
#include "generalfit.h"
#include "generaldata.h"
#include "bruit.h"

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
class DecaySineFuncbaseline
{
public:

    /** Constructor
        @param amp    amplitude of decaying sinusoid function 
        @param h      Hubble parameter in units of 100 km/s/Mpc
        @param decay  decay length                                            */
	DecaySineFuncbaseline(double h, double decay = 1.4) 
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
class FitBAObaselineScale
{
public:

	/** Constructor
	    @param pspectrum    observed power spectrum: columns(1,2,3) = (k,P(k),sigma_P(k))
	    @param simu_mode    if false, read shot noise and sigma */
  FitBAObaselineScale(TArray<r_8> pspectrum, bool simu_mode, double maxk, double mink, double h);

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
	void ComputeChisq(int degmax);
	
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

	//----------------RAJOUT PERSO------------------//
	/** Return baseline from power spectrum with inputs s and polynomial degree
	    @param SpecTheo                                                       */
	void PSpectSmoothing(double s, int deg);

	/** Write table containing the smoothed spectrum to a file 
	    @param outfile   file to write table to                           */
	void WriteTable(string outfile);

	/** Write table containing the avrgd points at each oscillation to a file 
	    @param outfile   file to write table to                           */
	void WriteTableAvrg(string outfile);

	void ReadBaseline(string ref_file);

	void FitBaseline(double s, int deg);

	double PSFunc(double const* x, double const* p);

	/** Destructor */
	virtual ~FitBAObaselineScale() {};

	//--------FIN RAJOUT PERSO--------------------//

protected:
	TVector<r_8> kobs_;      /**< k values of observed power spectrum         */
	TVector<r_8> Pobs_;      /**< observed power spectrum                     */
	TVector<r_8> sig_;       /**< error on observed power spectrum            */
	TVector<r_8> Pref_;      /**< smooth reference power spectrum             */
	TVector<r_8> Pratio_;    /**< ratio between observed and reference power spectra 		  */
	TVector<r_8> Svals_list_;/**< grid of trial s values                     */
	TVector<r_8> Avals_list_;/**< grid of trial A values                      */
	TArray<r_8> Chisq_;     /**< chi-square value at each A,S value           */
	TVector<r_8> Chisq1D_;  /**< chi-square value at each S value, 1 row of A           */
	TArray<r_8> SpecTheo_;
	TArray<r_8> pspectrum_;
	TVector<r_8> Errpoly_;
	TArray<r_8> avrg_;     
	
	double Smin_;           /**< min s value in trial s grid               */   
	double Smax_;           /**< max s value in trial s grid               */ 
	int nS_;                /**< number of s values in trial s grid        */
	double Amin_;		 /**< min Amp value in trial Amp grid             */  
  	double Amax_; 		 /**< max Amp value in trial Amp grid             */  
   	int nA_;		 /**< number of Amp value in trial Amp grid       */  
	int degree_;              /**< degree of the polynome used for fitting     */
	double mink_;
	double maxk_;
   
	Poly myfitpoly_;          /**< polynomial fit with coeff and errors        */
	Poly mybestfitpoly_;      /**< best (min chi2) polynomial fit with coeff and errors        */
	double h_;               /**< Hubble parameter in units of 100 km/s/Mpc   */
	double bestfitS_;         /**< best fit s                                  */
	double bestfitA_;     /**< best fit amp (adeline)      		  */
	double errA_;          /**< error on fit amp (adeline)      		  */
	double bestChisq_;     /**< best fit amp (adeline)      		  */
	double errup_;           /**< upper error range                           */
	double errdown_;         /**< lower error range                           */
	int nsig_;               /**< n-sigma of error ranges                     */
	bool simu_mode_;         /**< if false, read shot noise and sigma         */
};



#endif


