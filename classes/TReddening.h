

#ifndef TREDDENING_H_
#define TREDDENING_H_

//_______________________________
//
// TReddening.h
//
//  Created on: 3 aout 2011
//      Author: alexiagorecki
//
//      Redddening extinction laws
//      Cardelli law for early and late type galaxies
//      Calzetti law for starbursts galaxies Im, SB2, SB3
//
//      J.�A. Cardelli, G.�C. Clayton, and J.�S. Mathis. The relationship between
//		infrared, optical, and ultraviolet extinction. Astrophys. J., 345:245�256, 1989.
//      D.�Calzetti et�al. The Dust Content and Opacity of Actively Star-Forming Galaxies.
//	 	Astrophys. J., 533:682�695, 2000.
//
//	ex:
//		TReddening*	reddening = new TReddening();
//		reddening -> SetRVCardelli();
//		reddening -> SetRVCalzetti();

#include <string>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <vector>

#include "TMath.h"
using namespace std;

class TReddening {
public:
	//constructor
	TReddening();
	//destructor
	~TReddening();
	// set value of Rv for Cardelli Law, default is 3.1
	void SetRVCardelli(double RV = 3.1);
	// set value of Rv for Calzetti Law, default is 4.05
	void SetRVCalzetti(double RV = 4.05);
	// set Prevot law
	double Late(double l);
	double Prevot(double l);
	double Prevot_law(double l);
	double Prevot_file(double l);
	// set Calzetti law
	double Cardelli(double l);
	// set Calzetti law
	double Calzetti(double l);
	double Calzetti_newSB1(double l);
	double Calzetti_newSB2(double l);
	double Calzetti_newSB3(double l);
	// set Calzetti law
	double Calzetti_old(double l);
	// return the value of  pow(10., -0.4 * ebv * law) where law is given
	// by Cardelli or Calzetti
	double AttenuationFactor(double ebv, double lambda, string sedext);

	double fRV_card;
	double fRV_calz;
	double tempeval;

private:
	
	double a_lowl(double *, double *);
	double b_lowl(double *, double *);
	double a_lowVL(double *, double *);
	double b_lowVL(double *, double *);
	double a_FarUV(double *, double *);
	double b_FarUV(double *, double *);
	double cardelli_highl(double *lambda, double *par);
	
	void ReadTabCalz1();
	void ReadTabCalz2();
	void ReadTabCalz3();
	void ReadTabPrev();
	

	vector<double> tab_ext_calz1;
	vector<double> tab_lambda_calz1;
	vector<double> tab_ext_calz2;
	vector<double> tab_lambda_calz2;
	vector<double> tab_ext_calz3;
	vector<double> tab_lambda_calz3;

	vector<double> tab_ext_prev;
	vector<double> tab_lambda_prev;

};

#endif /* TREDDENING_H_ */
