/*
 * TBFilter.C
 *
 *  Created on: Sep 19 2014
 *      Author: JS Ricol
 */

#include "TBFilter.h"
#include "TDataCard.h"
#include "TFilter.h"

TBFilter::TBFilter(TDataCard *card)
{
  SetConversionUnitLambda(card->bfilterlambdaunit);
  SetFileType("tf1");
  string *bname = new string[1];
  bname[0] = card->Getbfiltername();

  //cout << "B Filter : " << card->Getbfiltername() << endl;

  LoadFile(1, bname);
}

TBFilter::~TBFilter() {
  ;
}
