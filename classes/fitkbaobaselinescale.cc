#include "fitkbaobaselinescale.h"


FitBAObaselineScale::FitBAObaselineScale(TArray<r_8> pspectrum, bool simu_mode, double maxk, double mink, double h)
  : simu_mode_(simu_mode)
{
  // fill power spectrum vectors
  InitVect(pspectrum);
  sa_size_t nc=pspectrum_.SizeX();
  sa_size_t nr=pspectrum_.SizeY();
  sa_size_t sz[2];
  sz[0]=nc;sz[1]=nr;
  SpecTheo_.SetSize(2,sz); // Define the SpecTheo object's size, protected var in .h
  maxk_=maxk; mink_=mink; 
  h_=h;

};

void FitBAObaselineScale::InitVect(TArray<r_8> pspectrum)
{
  pspectrum_=pspectrum;
  int nk = pspectrum_.SizeY(); // 2nd dim is k value direction
  //cout << "     Power spectrum has "<< nk <<" k values"<<endl;
  
  // set sizes of vectors
  kobs_.SetSize(nk); Pobs_.SetSize(nk); sig_.SetSize(nk);
  Pref_.SetSize(nk);
  Pratio_.SetSize(nk);
  
  double dPmc; //dPmc is added to the Pk value read in the table

  // col 0 = k ; col 1 = P(k), col 3 = shot noise, col 4 = err
  for (int i=0; i<nk; i++) {
    kobs_(i) = pspectrum_(0,i);

    if (simu_mode_) {
      Pobs_(i) = pspectrum_(1,i);
      sig_(i)  = 1.;
    } else {
      //Pobs_(i) = pspectrum(1,i) - pspectrum(6,i); 
      Pobs_(i) = pspectrum_(1,i)-pspectrum_(2,i); // RAJOUT PERSO
      //sig_(i)  = pspectrum(7,i);
      sig_(i) = pspectrum_(3,i); // RAJOUT PERSO
    }
    //cout <<	kobs_(i) << "    "<<Pobs_(i) << "    "<< sig_(i) <<endl;
  }
  
};

void FitBAObaselineScale::ComputeChisq( int degmax)
{
  int compt1=0,compt2=0;
  // set default s range
  Smin_=120.0; Smax_=180.; nS_=200; // no need to very fine grid: interpolated to compute the error
  // set default A range
  Amin_ = 1.; Amax_ = 3.; nA_ =200; // with k in the decay function
  // Amin_ = 0.1; Amax_ = 1.6; nA_ = 100; // with sqrt(k) in the decay function
  
  bestfitS_=-10; // uninitialised
  bestfitA_=-10; // uninitialised
  //FitBaseline(150.6,deg);
  double dS = (Smax_-Smin_)/(nS_-1);
  double dA = (Amax_ - Amin_)/(nA_-1);
  double sVals, aVals;   		
  
  double vecSize = nS_ * nA_;
  DecaySineFuncbaseline decs(h_);
  
  int ndim = 2;
  sa_size_t mydim[ndim];
  mydim[0] = nA_;
  mydim[1] = nS_;
  Chisq_.SetSize(ndim, mydim);
  Chisq1D_.SetSize(nS_);
  Avals_list_.SetSize(nA_);
  Svals_list_.SetSize(nS_);

  //cout << "    Computing chisq from s = "<< Smin_ <<" to s = ";
  //cout << Smax_ <<" in steps of "<< dS << " (ns="<< nS_ <<")" <<endl;
  //cout << "    A is variable : min = " << Amin_ << " to max = " << Amax_ << " in step of "<< dA << endl;
  //cout << "    Total :  compute "<< vecSize << " values of chisq "<< endl;  
  
  //Put the data in a specific format to fit
  int N=kobs_.Size();
  int n=0;
  for(int i=0;i<N;i++)
    {
      if(kobs_(i)>=mink_ && kobs_(i)<=maxk_){ n=n+1;}
    }
  double kdata,Pdata,sigdata;
  GeneralFitData mydata(1,n);  

  double PSFunc(double const* x, double const* p);
  double sloop=0.0,chisq=0.0;
  int Nddl,deg=0,degloop=0;
  bestChisq_=10000000.0; 

  for (int iS=0; iS<nS_; iS++) {
      sVals = Smin_ + iS*dS;
      Svals_list_(iS) = sVals;
      //degmax=2*(maxk_-mink_)/(2.*PI/sVals);
      for (int i=4;i<degmax;i++){
	deg=i;
	
	int kindx = 0;
	
	FitBaseline(sVals,deg);
	mydata.SetDataPtr(0);
	chisq=0.0;
	
	//Create data structure for GeneralFit
	for(int i=0;i<N;i++)
	  { kdata=kobs_(i);Pdata=Pobs_(i)/Pref_(i);sigdata=(sig_(i)*1.2)/Pref_(i);
	    if(kdata>=mink_ && kdata<=maxk_){
	      mydata.AddData(&kdata,Pdata,sigdata);
	      chisq+=pow((Pdata-decs(kdata,2.05127,149.846))/sigdata,2);
	    }
	  }
	
	//Fit the oscillations with "myfunc"
	GeneralFunc myfunc(1,2,PSFunc);
	GeneralFit myfit(&myfunc);
	myfit.SetParam(0,2,0.01,1,3);
	myfit.SetParam(1,150,0.01,120,180);
	myfit.SetData(&mydata);
	myfit.SetEps(1.e-9);
	myfit.Fit();
	Nddl=myfit.GetNddl();
	//myfit.ReCalChi2(Nddl);
	//myfit.PrintFit();//Parm();
	//cout<<"Param ="<<myfit.GetParm()<<endl;//"          s =  parm[1] ="<<myfit.GetParm(1)<<endl;
	//cout<<" Chi2 = "<<myfit.GetChi2()<<endl;
	if(myfit.GetChi2Red()<bestChisq_ && myfit.GetChi2()<1000 && myfit.GetChi2()>0.01)
	  {
	    sloop=sVals; degloop=deg;
	    bestfitA_=myfit.GetParm(0); 
	    bestfitS_=myfit.GetParm(1);
	    errup_=myfit.GetParmErr(1);   
	    errdown_=myfit.GetParmErr(1);
	    errA_=myfit.GetParmErr(0);
	    mybestfitpoly_.Realloc(degree_);
	    mybestfitpoly_=myfitpoly_;
	    bestChisq_=myfit.GetChi2Red();
	  }
      }
  }
  cout<<"Fit res  A = "<<bestfitA_<<" +/- "<<errA_ << "(so detection at " << bestfitA_/errA_ << " sigma ;     s = "<<bestfitS_<<" +/- "<<errup_<<"       Chi2 = "<<bestChisq_<<endl; 
  cout<<" sVals = "<<sloop<<endl;

  // Recompute the baseline with the best-fit parameters
  FitBaseline(sloop,degloop);
  
};


// calculate best-fit and errors
void FitBAObaselineScale::BestfitStdDev(double& bestfit_S, double& bestfit_A, double& bestfit_Chi, double& siglow, double& sighigh, double& sigbest, int nsig)
{
   nsig_ = nsig;
   if (nsig_>3)
	throw ParmError("Can't compute above 3-sigma");

   vector<double> clevels;
   clevels.push_back(0.683); clevels.push_back(0.954); clevels.push_back(0.997);

   double clevel = clevels[nsig_-1];
   ChisqStats chisqstat(Svals_list_, Avals_list_, Chisq_,  2); 
   bestfit_Chi = chisqstat.BestFit(bestfit_A, bestfit_S);

   bestfitA_  = bestfit_A;
   bestfitS_  = bestfit_S;
   bestChisq_ = bestfit_Chi;
   for (int j=0; j< Chisq_.SizeY(); j++) Chisq1D_(j) = Chisq_(bestfit_A,j);
   cout << "bestfitA_ : " << bestfitA_ << endl;
   cout << "raw bestfitS_ : " << bestfitS_ << endl;
 
   //cout << endl << " marginalize chisq "<< endl;
   double dA = (Amax_ - Amin_)/nA_;
   double ds =  (Smax_ - Smin_)/nA_;
   double step[] = {ds, dA};
   //cout <<" step "<< ds << "  "<< dA << "  "<< step[0] << "  "<< step[1] << endl;
   ChisqStats chisqstat_forMarg(Svals_list_, Chisq1D_);

   //chisqstat.GetMarg(step, margChisq_);
   //Computing Errors
   //cout <<" step "<< ds << "  "<< dA << "  "<< step[0] << "  "<< step[1] << endl;

   chisqstat_forMarg.ErrSig(siglow,sighigh,sigbest,clevel,1000);
   bestfitS_ = sigbest;
   errup_   = sighigh - bestfitS_;
   errdown_ = bestfitS_ - siglow;
   //cout << "Fit results: s = "<< bestfitS_ << " + "<< errup_  << " ; -  "<< errdown_ << endl;
};


// write chi-square to a file
void FitBAObaselineScale::WriteChisq(string outfile)
{
   if (bestfitS_<0)
	throw ParmError("Have not calculated best fit s yet!");

   ifstream inp;
   ofstream outp;
   inp.open(outfile.c_str(), ifstream::in);
   inp.close();
   if(inp.fail()) {
     inp.clear(ios::failbit);
     cout << "    Writing chisq to file ..." << outfile.c_str() << endl;
     outp.open(outfile.c_str(), ofstream::out);
     
     outp <<" Best fit s = ";
     outp << bestfitS_ <<", "<< nsig_ <<"-sigma errors: + "<< errup_;
     outp <<" - "<< errdown_ << " with h = " << h_ << endl;
     
     outp << "amplitude is variable, best fit A = "<< bestfitA_ <<endl;
     outp << "first column : Avals(nA = "<<Chisq_.SizeX()<<"), second column : Svals(nS = "<<Chisq_.SizeY()<<"), third column : Chisq" <<endl;

     //wirte first line (A values)
     for(int i=0; i<Avals_list_.Size(); i++) {// loop over A values  
       for(int j=0; j<Svals_list_.SizeY(); j++){// loop over S values
	 outp << Avals_list_(i) << " \t" ;
	 outp << Svals_list_(j) << " \t";
	 outp << Chisq_(i,j) << "\t ";   
	 outp << endl;
       }
     }	
     outp.close();
   } // end of write
   else
     cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

};

// write results
void FitBAObaselineScale::WriteResults(string outfile)
{

  ifstream inp;
  ofstream outp;
  inp.open(outfile.c_str(), ifstream::in);
  inp.close();
  if(inp.fail()) {
    
    inp.clear(ios::failbit);
    //cout << "    Writing results to file ..." << outfile.c_str() << endl;
    outp.open(outfile.c_str(), ofstream::out);

    outp << "redshift : best fit s : n-sigma of errors : +error : -error :  best fit A : hubble parameter h"<<endl;
    
    outp << bestfitS_ <<"   "<<"   "<< errup_;
    outp <<"   "<< errdown_ <<"   "<< bestfitA_ << "   " << bestChisq_ << "   " << h_ << endl<<endl;
    for(int j=0;j<=degree_;j++){
      outp<<mybestfitpoly_[j]<<"*X^"<<j;  //"  ± "<<errCoeff(j)<<endl;
    }

    outp.close();
  } // end of write
  else
    cout << "Error...file """ << outfile.c_str() << """ exists" << endl;
  
};




// RAJOUT PERSO

// write SpecTheo to a file
void FitBAObaselineScale::WriteTable(string outfile)
{

  ifstream inp;
  ofstream outp;
  inp.open(outfile.c_str(), ifstream::in);
  inp.close();
  if(inp.fail()) {
    
    inp.clear(ios::failbit);
    outp.open(outfile.c_str(), ofstream::out);
    
    sa_size_t nc=SpecTheo_.SizeX();
    sa_size_t nr=SpecTheo_.SizeY();
    int i,j;
    
    for(i=0;i<nr;i++){
      for(j=0;j<nc;j++){
	outp <<"  "<<SpecTheo_(j,i)<<"  ;";
      }
      outp<<endl;
    }
    outp.close();
  } // end of write
  else
    cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

};

void FitBAObaselineScale::WriteTableAvrg(string outfile)
{

  ifstream inp;
  ofstream outp;
  inp.open(outfile.c_str(), ifstream::in);
  inp.close();
  if(inp.fail()) {
    
    inp.clear(ios::failbit);
    outp.open(outfile.c_str(), ofstream::out);
    
    sa_size_t nc=avrg_.SizeX();
    sa_size_t nr=avrg_.SizeY();
    int i,j;
    
    for(i=0;i<nr;i++){
      for(j=0;j<nc;j++){
	outp <<avrg_(j,i)<<"  ";
      }
      outp<<endl;
    }
    outp.close();
  } // end of write
  else
    cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

};

//Construct baseline from power spectrum
void FitBAObaselineScale::PSpectSmoothing(double s, int deg)//, TArray<r_8>& SpecTheo)
{
  sa_size_t nc=pspectrum_.SizeX();
  sa_size_t nr=pspectrum_.SizeY();
  int j=0;
  double kloop=0.,dk,ka=2.*PI/s;
  sa_size_t N= int (2*(maxk_-mink_)/ka+2);
  Vector kaverage(N);

  // selecting the range for the k vector with ka/2 step
  while(j*(ka/2)<mink_){ kaverage(0)=j*(ka/2); j=j+1;}
  for(int i=1;i<N;i++){kaverage(i)=kaverage(0)+i*ka/2;}
  
  // interpolation of the power spectrum and error data
  vector<double> kvals,poscvals,perrvals;
  for(int kk=0; kk<nr; kk++) {
    double kv=kobs_(kk);
    kvals.push_back(kv);
    double poscv=Pobs_(kk)*kobs_(kk)*kobs_(kk);
    poscvals.push_back(poscv);
    double perrv=sig_(kk);
    perrvals.push_back(perrv);

  }
  //cout<<endl<<endl;
  size_t npt=0;
  SLinInterp1D interpPS(kvals,poscvals,mink_,maxk_,npt);
  SLinInterp1D interpERR(kvals,perrvals,mink_,maxk_,npt); 
  // average for each j*ka/2, on an interval : [-ka/4;+ka/4]

  int Np=100;
  double pavtot=0.,erravtot=0., deltak=ka/(2*Np);
  Vector paverage(N), erraverage(N), err2average(N);

  for(int n=0;n<N;n++){
    pavtot=0.0; erravtot=0.0;
    for(int i=0;i<Np;i++){
      pavtot=pavtot+interpPS(kaverage(n)+i*deltak-ka/4);
      erravtot=erravtot+interpERR(kaverage(n)+i*deltak-ka/4);
      
    }
    paverage(n)=pavtot/Np;
    err2average(n)=(erravtot/Np)*(erravtot/Np);
  }
  
  sa_size_t avrg_nc=2;
  sa_size_t avrg_nr=N;
  sa_size_t avrg_sz[2];
  avrg_sz[0]=avrg_nc;avrg_sz[1]=avrg_nr;
  avrg_.SetSize(2,avrg_sz);
  
  for(int i=0;i<N;i++){ avrg_(0,i)=kaverage(i); avrg_(1,i)=paverage(i);}
  
  //--- Create a ploynomial objet
  degree_=deg; //degree is a protected variable defined in .h
  //"Poly myfitpoly" protected variable defined in .h
  Vector errCoeff(degree_+1);

  //  fit the polynomial on the data
  double xi2 = myfitpoly_.Fit(kaverage, paverage, err2average, degree_, errCoeff);

  //--- print the result
  for(int i=0; i<=degree_; i++) {
    // cout<<"FitRes: i="<<i<<"  --> fitpoly[i]="<<myfitpoly_[i]<<"  ± "<<errCoeff(i)<<endl;
  }

  if(simu_mode_){
 //On créer de nouveaux vecteurs k et Pwosc(k) avec pas en k constant (Pour spectres CLASS)
  double kmaxCLASS=0.3, kminCLASS=0.01;
  double dkCLASS = (kmaxCLASS-kminCLASS)/nr;
  double X=0.0;
  for(int i=0;i<nr;i++){                                     
    X=myfitpoly_[0];
    SpecTheo_(0,i)=i*dkCLASS+kminCLASS;
    SpecTheo_(1,i)=interpPS(i*dkCLASS)/SpecTheo_(0,i)/SpecTheo_(0,i);
    for(int ii=1;ii<=degree_;ii++){ X=X+myfitpoly_[ii]*pow(SpecTheo_(0,i),ii); }
    SpecTheo_(2,i)= X/(SpecTheo_(0,i)*SpecTheo_(0,i));
    //cout<<SpecTheo_(0,i)<<"   "<<SpecTheo_(1,i)<<"   "<<SpecTheo_(2,i)<<endl;
  }
  }
  else{
  double X=0.0;
  for(int i=0;i<nr;i++){                                     
    X=myfitpoly_[0];
    SpecTheo_(0,i)=kobs_(i);
    SpecTheo_(1,i)=Pobs_(i);
    for(int ii=1;ii<=degree_;ii++){ X=X+myfitpoly_[ii]*pow(SpecTheo_(0,i),ii); }
    SpecTheo_(2,i)= X/(SpecTheo_(0,i)*SpecTheo_(0,i));
  }
  }

  double Xmin=0.0,Xmax=0.0;
  Errpoly_.SetSize(nr);
  for(int i=0;i<nr;i++){                                     
    Xmin=myfitpoly_[0]-errCoeff[0]/myfitpoly_[0];
    Xmax=myfitpoly_[0]+errCoeff[0]/myfitpoly_[0];
    for(int ii=1;ii<=degree_;ii++){ 
      Xmin=Xmin+(myfitpoly_[ii]-errCoeff[ii]/myfitpoly_[ii])*pow(SpecTheo_(0,i),ii); 
      Xmax=Xmax+(myfitpoly_[ii]+errCoeff[ii]/myfitpoly_[ii])*pow(SpecTheo_(0,i),ii);
    }
    Errpoly_(i)= sqrt(pow((Xmax-Xmin),2))/2;
  }
};

// Call PSpectSmoothing to get the baseline and prepare the data for ComputeChisq
void FitBAObaselineScale::FitBaseline(double s, int deg)
{
  sa_size_t nr, nc;    
  PSpectSmoothing(s,deg); 
  nc=SpecTheo_.SizeX();
  nr=SpecTheo_.SizeY();
  vector<double> kvals,pvals,poscvals;
  int icolk=0,icolposc=1,icolp=2;
  for(int kk=0; kk<nr; kk++) {
    double kv=SpecTheo_(icolk,kk);
    kvals.push_back(kv);
    double poscv=SpecTheo_(icolposc,kk);
    poscvals.push_back(poscv);
    double pv=SpecTheo_(icolp,kk);
    pvals.push_back(pv);
  }

  for (int i=0; i<kobs_.Size(); i++) { Pref_(i)=SpecTheo_(i);}

  // interpolate this spectrum at kobs values
  size_t npt=0;
  SLinInterp1D interpYR(kvals,pvals,mink_,maxk_,npt);//interpYR(kvals[0], kvals[nr-1],pvals);
  //cout << "------ kobs ------ Pwosc ----- Psmooth ----"<< endl;
  for (int i=0; i<kobs_.Size(); i++) {
    Pref_(i) = interpYR(kobs_(i));// ORIGINAL
    //cout<< "   "<<kvals[i]<<"   "<<poscvals[i]<<"   "<<pvals[i]<< endl;
  }
  

  vector<double> perr;
  for(int kk=0; kk<nr; kk++) {
    double pe=Errpoly_(kk);
    perr.push_back(pe);
  }

  SLinInterp1D interpERR(kvals,perr,mink_,maxk_,npt);
  for (int i=0; i<kobs_.Size(); i++) {
    Errpoly_(i) = interpERR(kobs_(i))/(Pref_(i)*kobs_(i)*kobs_(i));
  }

  
};

//Function to fit the oscillations
  double PSFunc(double const* x, double const* p)
  {

    DecaySineFuncbaseline decs(0.679);
    //decs.PrintParas();
    double osc = decs(x[0], p[0], p[1]);
    
    return osc;
    }

// FIN RAJOUT PERSO
