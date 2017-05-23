#include "chisqstats.h"

double ChisqStats::BestFit(double& bestamp, double& bests)
{
  int Sindexmin, Aindexmin;
  double minchisq = findMinimumPosition(Chisq_, Aindexmin, Sindexmin);
  cout <<"    Value of min(chisq) = "<< minchisq <<endl;
  cout <<"    Value of s at min(chisq) = "<< Svals_list_(Sindexmin) <<endl;
  cout <<"    Value of amp at min(chisq) = "<< Avals_list_(Aindexmin) <<endl;
  bestamp = Avals_list_(Aindexmin);
  bests = Svals_list_(Sindexmin);
  return minchisq;
};

void ChisqStats::GetMarg(double *step_, TArray<r_8> MargChisq)
{ 
    cout << "  marginalize chisq" << endl;
    int nA = Chisq_.SizeX();
    int nS = Chisq_.SizeY();
    
    cout << "   "<< nA << " x "<< nS <<" val in chisq "<< endl; 
    cout << "  dS "<< step_[0] << " dA "<<step_[1] << endl; 
    
    int nVar = int(dof_);
    if(dof_>2)
    	cout << "ERROR dof > 2" <<endl;
    int dim[nVar];
    double step;
    
    for(int ivar=0; ivar<nVar; ivar ++){
      if(ivar == 0){	//marginalization over A => chi(S)
	dim[0] = Chisq_.SizeY();
	dim[1] = Chisq_.SizeX();
	step = step_[1]; //dA
      }if(ivar == 1){	//marginalization over S => chi(A)
	dim[0] = Chisq_.SizeX();
	dim[1] = Chisq_.SizeY();
	step = step_[0]; //dS
      }
      
      double content = 0;
      for(int i=0; i<dim[0]; i++){ 
	content = 0;
	for(int j=0; j<dim[1]; j++){ 
	  if(ivar == 0)
	    content += (Chisq_(j,i) * step);	 
	  if(ivar == 1){
	    content += (Chisq_(i,j) * step);
	  }
	}
	MargChisq(ivar, i) = content;
	
      }
    }
}

void ChisqStats::ErrSig(double& siglow, double& sighigh, double& sigbest, double clevel, int npt)
{
// xvals_= abscissa for logP (must be evenly spaced)
// logP  = vector of probabilities
// clevel= e.g. 0.68 (is default)
// npt is number of points to interpolate with
// using log(prob) instead of prob
// avoids more NaNs and therefore gives a more accurate answer
    int ylevelstep=100;
    TVector<r_8> logPraw = -0.5*Chisq1D_;

    // REMOVE ANY XVALUES WHERE LOGP=INF
    
    // count number of NOT Infs
    int cnt=0;
    for (int i=0; i<xvals_.Size(); i++) {
	    if(!isinf(logPraw(i)))
	        cnt++;
	    else 
	        cout <<"Inf at x="<< xvals_(i) <<endl;
	    }
    
    // make and fill new vectors which don't have any Inf values for log(P)
    TVector<r_8> logP(cnt), xv(cnt);
    int ii=0;
    for (int i=0; i<xvals_.Size(); i++) {
    
	    if(!isinf(logPraw(i))) {
	        xv(ii)=xvals_(i);
		logP(ii)=logPraw(i); 
		ii++;
            }
    }
	    
    // CHECK LOGP IS NOT FLAT OR NAN AND RENORMALISE
    double logPmax, logPmin;
    logP.MinMax(logPmin,logPmax);
    //cout <<"logPmax="<<logPmax<<"logPmin="<<logPmin<<endl;

    // check that logP is not flat
    if (logPmax==logPmin) return;
    
    // check that logP is not all NaN's
    int nan_cnt = 0;
    for (int i=0; i<logP.Size(); i++)
        if (my_isnan(logP(i)))
            nan_cnt++;
    if (nan_cnt == logP.Size())
        return;

    // Some problem occurs when logP values are vv small
    // Fix this by using normalised pdf where possible:
    logP=logP - logPmax;

    // INTERPOLATE XV ONTO FINER GRID X
    int ncoarse = xv.Size(); // old size of xvals
    int n = ncoarse*npt; // new size of xvals: renamed x
    double xstep = (xv(ncoarse-1) - xv(0))/n;
    //cout <<"xstep="<<xstep<<", n="<<n<<endl;

    SInterp1D interplogP(xv,logP,xv(0),xv(ncoarse-1),n);
    TVector<r_8> x(n), y(n); // interped xv and logP vectors
    for (int i=0; i<n; i++) {
	    x(i)=xv(0) + i*xstep;
        y(i)=exp(interplogP(x(i))); // interpolate logP=-chisq/2 then do exp
	    }

    // FIND MAX PROBABILITY AND INDEX POSITION
    double ymin, ymax;
    y.MinMax(ymin,ymax);
    int imax = findClosestElement(y, ymax);
    if(imax==(n-1))
	    cout << "    IMAX IS LAST INDEX! min y="<< ymin <<", max y="<< ymax 
	         << ", (imax="<< imax <<")"<<endl;

    // check none of y is a nan
    //double sum=0;
    /*for (int i; i<n; i++)
	{
	if ( isnan(y(i)) ) {cout <<"y has NaNs"<<endl;}
	sum+=y(i);
	}*/
    //cout <<"sum="<<sum<<", xstep"<<xstep<<", sum2="<<y.Sum()<<endl;


    // STEP THROUGH LEVELS IN Y, EVALUATING CONFIDENCE IN EACH
    double norm = y.Sum()*xstep;
    //cout << "     norm = "<<norm<<endl;
    double ylevel=ymax;
    double conf=0;
    //make a vector to integerate which is the same as y, but zero where y.lt.lim
    TVector<r_8> tempvect;
    while (conf<clevel) {
    
        ylevel = ylevel-ymax/ylevelstep;  
        tempvect = y;
        // where y<ylevel set to 0
        for (int i=0; i<n ; i++)
	        if(y(i)<ylevel) tempvect(i)=0;
 
        // integrate it, normalise, compare with confidence level reqd
        conf = (tempvect.Sum()*xstep)/norm;
        }

    // START AT PEAK AND STEP BACKWARDS ALONG X AXIS UNTIL y<ylevel
    TVector<r_8> y2peak(imax+1);
    for (int k=0; k<(imax+1); k++)
	    y2peak(k)=y(k);
	
    int i = findClosestElement(y2peak, ylevel); 
    minus_ = x(i);

    // START AT PEAK AND STEP FORWARD ALONG X AXIS UNTIL y<ylevel
    if(imax<(n-1)) {
	
	    TVector<r_8> yfrompeak(n-(imax+1));
	    for (int k=0; k<(n-(imax+1)); k++)
		    yfrompeak(k)=y(k+(imax+1));
		
	    i = findClosestElement(yfrompeak, ylevel);
	    plus_=x(i+(imax+1));
	    }
	else {
	    plus_=-9999;
		cout <<"    WARNING: probability peaks at last index"<<endl;
		}

    siglow = minus_;
    sighigh = plus_;
    sigbest = x(imax);
};


/*int ChisqStats::NearestIndex(double val,TVector<r_8> vect)
{

double diff=1e9;
int index;
for (int i=0; i<vect.Size(); i++)
	{
	if (abs(val-vect(i))<diff)
		{
		diff = abs(val-vect(i));
		index=i;
		}
	}
return index;

};*/
