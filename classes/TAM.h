
#ifndef TAM_h
#define TAM_h
//___________________________
/*
 * TAM.h
 *
 *  Created on: Sep 27, 2010
 *      Author: gorecki
 *
 *      cf Madau  1995/2000
 *
 *      Modified on: Sep 19 2014 by JS Ricol
 *

 */

#include <vector>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
/*   CHECK-REZA-JS  : on mets les includes ROOT ds le .cc uniquement 
     TAM.h peut etre inclus sans conflit avec SOPHYA 
#include <string.h>
#include <TROOT.h>
#include "TRandom3.h"
#include "TDataCard.h"
#include "TApparentMagnitude.h"
#include "TKcorrection.h"
*/


// Forward declarations for ROOT classes 
class TRandom3;
class TApparentMagnitude;

using namespace std;

class TAM {

public:
  // constructor
  TAM(); 
  // destructor
  ~TAM();
  bool PassGoldenCut(double mag_abs, int type, double z);
  
 private:
  TRandom3 *rand;
  TApparentMagnitude *tam;
  double *appmag;
  double *err_appmag;
  
};
#endif /* TAM_h*/
