/*
 * TAM.C
 *
 *  Created on: Sep 27, 2010
 *      Author: gorecki
 */

#include "TAM.h"

#include <string.h>

#include "TRandom3.h"

#include "TDataCard.h"
#include "TApparentMagnitude.h"
#include "TKcorrection.h"


TAM::TAM(){
  rand = new TRandom3();
  rand->SetSeed(0);
  
  // Read DataCard
  TDataCard *datacard = new TDataCard();
  if (datacard->Read("/afs/in2p3.fr/home/j/jricol/PhotoZ/dev/data/DataCards/bao")<=0){
    cout << "problem with datacard ... exit program" << endl;
  }
  else 
    datacard->Print();
  
  TKcorrection *kcorr = new TKcorrection();
  string kfile = Form("%s/%s",datacard->tabledir.c_str(),datacard->kcorrfile.c_str());
  kcorr -> LoadTable(kfile);
  
  tam = new TApparentMagnitude();
  tam -> LoadDataCard(datacard);
  tam -> SetKcorrectionTable(kcorr);
  tam -> SetParRange(kcorr->GetParRange());
  appmag = new double[6];
  err_appmag = new double[6];
}

TAM::~TAM() {
  delete appmag;
  delete err_appmag;
}

bool TAM::PassGoldenCut(double mag_abs, int type, double z){
  double ebv = tam->GetEBV(type);
  tam->EvalApparentMagnitudes(appmag, err_appmag, mag_abs, type, z, ebv, 3); // Dahlen etc: LF in mag_B
  // tam->EvalApparentMagnitudes(appmag, err_appmag, mag_abs, type, z, ebv, -1); // Ramos: LF in mag_i
  // cout << "AM : " << appmag[3] << " (type = " << type << ", z = " << z << ")" << endl;
  bool pass = appmag[3]<25.3;
  // bool pass = appmag[3]<22.5; // TEST PROVISOIRE pour comparaison Zucca
  // if (pass) cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!! PASS = " << pass*1 << endl;
  return pass;
}

