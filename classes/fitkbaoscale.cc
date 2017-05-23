#include "fitkbaoscale.h"

FitBAOScale::FitBAOScale(TArray<r_8> pspectrum, SimpleUniverse& su, double zref, string ref_file, bool simu_mode)
  :su_(su),  zref_(zref), simu_mode_(simu_mode)
{
	// fill power spectrum vectors
	InitVect(pspectrum);

	// Cosmological parameters
	h_=su_.h();
	cout << "     h = "<< h_ <<endl;
	// read computed power spectrum without oscillations
	ifstream ifs(ref_file.c_str());
	TArray<r_8> SpecTheo;
	sa_size_t nr, nc;
	SpecTheo.ReadASCII(ifs,nr,nc);
	vector<double> kvals,pvals;
	int icolk=0,icolp=2;
	for(int kk=0; kk<nr; kk++) {
	  double kv=SpecTheo(icolk,kk);
	  kvals.push_back(kv);
	  double pv=SpecTheo(icolp,kk);
	  pvals.push_back(pv);
	}

	// interpolate this spectrum at kobs values
	SLinInterp1D interpYR(kvals[0], kvals[nr-1],pvals);
	for (int i=0; i<kobs_.Size(); i++) Pref_(i) = interpYR(kobs_(i));

	// Compute Pobs/Pref
	cout << "compute Pratio" << endl;
	Pratio();
	
	// set default s range
	Smin_=120.; Smax_=180.; nS_=200; // no need to very fine grid: interpolated to compute the error
	// set default A range
	Amin_ = 1.; Amax_ = 3.; nA_ = 200; // with k in the decay function
	// Amin_ = 0.1; Amax_ = 1.6; nA_ = 100; // with sqrt(k) in the decay function

	bestfitS_=-10; // uninitialised
	bestfitA_=-10; // uninitialised
	cout << " end" << endl;
};


void FitBAOScale::InitVect(TArray<r_8> pspectrum)
{
	int nk = pspectrum.SizeY(); // 2nd dim is k value direction
	cout << "     Power spectrum has "<< nk <<" k values"<<endl;
	
	// set sizes of vectors
	kobs_.SetSize(nk); Pobs_.SetSize(nk); sig_.SetSize(nk);
	Pref_.SetSize(nk);
	Pratio_.SetSize(nk);

	// col 1 = total spectrum, col 6 = shot noise spectrum, col 7 = sigma (computed from (total - shot noise) spectrum (Cecile)
	for (int i=0; i<nk; i++) {
	  kobs_(i) = pspectrum(0,i);
	  if (simu_mode_) {
	    Pobs_(i) = pspectrum(1,i);
	    sig_(i)  = 1.;
	  } else {
	    Pobs_(i) = pspectrum(1,i) - pspectrum(6,i); //pspectrum(1,i);
	    sig_(i)  = pspectrum(7,i);//pspectrum(2,i)
	  }
	  	  //if (i<100) cout <<	kobs_(i) << "    "<<Pobs_(i) << "    "<< sig_(i) <<endl;
	}
};


void FitBAOScale::ComputeChisq(double maxk)
{
   double dS = (Smax_-Smin_)/(nS_-1);
   double dA = (Amax_ - Amin_)/(nA_-1);
   double sVals, aVals;   		
	
   double vecSize = nS_ * nA_;
   DecaySineFunc decs(h_);
   decs.PrintParas();

   int ndim = 2;
   sa_size_t mydim[ndim];
   mydim[0] = nA_;
   mydim[1] = nS_;
   Chisq_.SetSize(ndim, mydim);
   Chisq1D_.SetSize(nS_);
   Avals_list_.SetSize(nA_);
   Svals_list_.SetSize(nS_);

   cout << "    Computing chisq from s = "<< Smin_ <<" to s = ";
   cout << Smax_ <<" in steps of "<< dS << " (ns="<< nS_ <<")" <<endl;
   cout << "    A is variable : min = " << Amin_ << " to max = " << Amax_ << " in step of "<< dA << endl;
   cout << "    Total :  compute "<< vecSize << " values of chisq "<< endl;  

   for(int iA=0; iA<nA_; iA ++){
     aVals = Amin_ + iA*dA;
     Avals_list_(iA) = aVals;

     for (int iS=0; iS<nS_; iS++) {
       sVals = Smin_ + iS*dS;
       Svals_list_(iS) = sVals;

       double chisq=0;
       int kindx = 0;
       
       for (int ik=0; ik<kobs_.Size(); ik++) {
	 if (kobs_(ik) >= maxk || kobs_(ik) <= 0.02) continue;

	 // only use power spectra values below some max k
	 if (sig_(ik)==0)
	   throw ParmError("Error is zero!");
	 
	 // decaying sinusoid at s 
	 double pred = decs(kobs_(ik), aVals, sVals);
	 double diff= Pratio_(ik)-pred;
	 
	 chisq += diff*diff/(sig_(ik)*sig_(ik));	 
	 kindx ++;
       }
       
       Chisq_(iA, iS) = chisq/(kindx-2); // fit with 2 parameters
     }
   }
   cout << "   check compute chisq :  " << endl;
   cout << "   Chisq_: x length="<<Chisq_.SizeX() << " y length="<< Chisq_.SizeY() << " shoulb be equal to (";
   cout << nA_ << ", "<< nS_ <<")"<<endl;
 
};


// calculate best-fit and errors
void FitBAOScale::BestfitStdDev(double& bestfit_S, double& bestfit_A, double& bestfit_Chi, double& siglow, double& sighigh, double& sigbest, int nsig)
{
   nsig_ = nsig;
   if (nsig_>3)
	throw ParmError("Can't compute above 3-sigma");

   vector<double> clevels;
   clevels.push_back(0.683); clevels.push_back(0.954); clevels.push_back(0.997);

   double clevel = clevels[nsig_-1];
   cout <<"    Confidence level = "<< clevel <<endl;
   ChisqStats chisqstat(Svals_list_, Avals_list_, Chisq_,  2); 
   cout <<"    Confidence level = "<< clevel <<endl;
   bestfit_Chi = chisqstat.BestFit(bestfit_A, bestfit_S);

   bestfitA_  = bestfit_A;
   bestfitS_  = bestfit_S;
   bestChisq_ = bestfit_Chi;
   for (int j=0; j< Chisq_.SizeY(); j++) Chisq1D_(j) = Chisq_(bestfit_A,j);
   cout << "bestfitA_ : " << bestfitA_ << endl;
   cout << "raw bestfitS_ : " << bestfitS_ << endl;
 
   cout << endl << " marginalize chisq "<< endl;
   double dA = (Amax_ - Amin_)/nA_;
   double ds =  (Smax_ - Smin_)/nA_;
   double step[] = {ds, dA};
   cout <<" step "<< ds << "  "<< dA << "  "<< step[0] << "  "<< step[1] << endl;
   ChisqStats chisqstat_forMarg(Svals_list_, Chisq1D_);

   //chisqstat.GetMarg(step, margChisq_);
   //Computing Errors
   cout <<" step "<< ds << "  "<< dA << "  "<< step[0] << "  "<< step[1] << endl;

   chisqstat_forMarg.ErrSig(siglow,sighigh,sigbest,clevel,1000);
   bestfitS_ = sigbest;
   errup_   = sighigh - bestfitS_;
   errdown_ = bestfitS_ - siglow;
   cout << "Fit resultst: s = "<< bestfitS_ << " + "<< errup_  << " ; -  "<< errdown_ << endl;
};


// write chi-square to a file
void FitBAOScale::WriteChisq(string outfile)
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
     
     outp << "Redshift of power spectrum = "<< zref_ <<", best fit s = ";
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
void FitBAOScale::WriteResults(string outfile)
{

  ifstream inp;
  ofstream outp;
  inp.open(outfile.c_str(), ifstream::in);
  inp.close();
  if(inp.fail()) {
    
    inp.clear(ios::failbit);
    cout << "    Writing results to file ..." << outfile.c_str() << endl;
    outp.open(outfile.c_str(), ofstream::out);
    
    cout << "    Redshift of power spectrum = "<< zref_ <<", best fit S = ";
    cout << bestfitS_ <<", best fit Amp = ";
    cout << bestfitA_ <<", "<< nsig_ <<"-sigma errors: + "<< errup_ <<" - ";
    cout << errdown_ <<endl;
    
    outp << "redshift : best fit s : n-sigma of errors : +error : -error :  best fit A : hubble parameter h"<<endl;
    
    outp << zref_ <<"   "<< bestfitS_ <<"   "<< nsig_ <<"   "<< errup_;
    outp <<"   "<< errdown_ <<"   "<< bestfitA_ << "   " << bestChisq_ << "   " << h_ << endl;
    
    outp.close();
  } // end of write
  else
    cout << "Error...file """ << outfile.c_str() << """ exists" << endl;
  
};




