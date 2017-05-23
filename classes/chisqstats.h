#ifndef CHISQSTAT_H_SEEN
#define CHISQSTAT_H_SEEN

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
#include "geneutils.h"
#include "constcosmo.h"

/** @class ChisqStats
  * Simple chisq calculation
  * For a 1D probability function finds the best fit value and the error
  */
class ChisqStats 
{
public:
	
  /** Constructor  find best fit when ka AND A are variable
      @param Svals     parameter values of chi-square (must be evenly spaced)
      @param Avals     parameter values of chi-square (must be evenly spaced)
      @param chisq     chi-square value for each parameter
      @param dof       degree of freedom                                    */
 ChisqStats(TVector<r_8> svals, TVector<r_8> ampvals,  TArray<r_8> chisq, double dof=1)
   : Svals_list_(svals), Avals_list_(ampvals), Chisq_(chisq), dof_(dof) {
    dof_ = dof;
    cout << "init chisq Array" << Chisq_.SizeX() << " " << Chisq_.SizeY()<<endl;
    // for (int i=0;i<10;i++) cout << "TEST 2D "<< i << " " << Chisq_(i,10) << endl;
  };
 ChisqStats(TVector<r_8> svals, TVector<r_8> chisq, double dof=1)
   : xvals_(svals), Chisq1D_(chisq), dof_(dof) {
    dof_ = dof;
    cout << "init chisq Vector " << Chisq1D_.Size() << endl;
    // for (int i=0;i<10;i++) cout << "TEST 1D "<< i << " " << Chisq1D_(i) << endl;

  };
  
  
  /** Find the best fit parameter value ie parameter value at the minimum 
      chi-square                                                            */
  double BestFit(double& bestamp, double& bests); 
  
  /** Estimate the parameter error 
      @param siglow    lower error range
      @param sighigh   upper error range
      @param clevel    confidence level of error range
      @param npt       number of points to interpolate chi-square func with */
  void ErrSig(double& siglow, double& sighigh, double& sigbest, double clevel=0.683,int npt=100); // find error
  //int NearestIndex(double,TVector<r_8>);
  
  /** Return chi-square function                                            */
  TArray<r_8> ReturnChisq(){return Chisq_;};
  
  /** Return parameter values of chi-square                                 */
  TVector<r_8> ReturnXvals(){return xvals_;};
  
  /** Return lower error range                                              */
  double ReturnMinus(){return minus_;};
  
  /** Return upper error range                                              */
  double ReturnPlus(){return plus_;};
  
  /** Marginalize chisq over Amp and ka (Adeline) 	 		  */
  void GetMarg(double *step_, TArray<r_8> MargChisq);
  
  
 protected:
	TArray<r_8> Chisq_; /**< chi-square value for each parameter         */
	TVector<r_8> Chisq1D_; /**< chi-square value for each parameter         */
	TVector<r_8> Svals_list_;/**< grid of trial s values                     */
	TVector<r_8> Avals_list_;/**< grid of trial A values                      */
	TVector<r_8> xvals_;     /**< grid of trial values                     */

	double dof_;             /**< degree of freedom                           */
	double minus_;           /**< lower error range                           */
	double plus_;            /**< upper error range                           */

};

#endif


