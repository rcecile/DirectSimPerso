
#ifndef TBFilter_h
#define TBFilter_h

//_________________________
///
// TBFilter.h
//
//  Created on: Mar 19, 2014
//      Author: JS Ricol
//

#include "TFilter.h"
#include "TDataCard.h"

using namespace std;

class TBFilter : public TFilter{
 public:
  //constructor
  TBFilter(TDataCard *card);
  // destructor
  ~TBFilter();
};

#endif /*TBFilter_h*/
