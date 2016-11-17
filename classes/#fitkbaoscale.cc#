#include "fitkbaoscale.h"

FitBAOScale::FitBAOScale(TArray<r_8> pspectrum, SimpleUniverse& su, double zref, string ref_file, bool simu_mode)
  :su_(su),  zref_(zref), simu_mode_(simu_mode)
{
	// fill power spectrum vectors
	InitVect(pspectrum);

	// Cosmological parameters
 	double OmegaM=su_.OmegaMatter();
	double OmegaB=su_.OmegaBaryon();
	double OmegaL=su_.OmegaLambda();
	double sig8=su_.Sigma8();
	double n=su_.Ns();
	double R=8; // sigma8 definition
	h_=su_.h();
	cout << "     h = "<< h_ <<endl;
	//	ComputeSmoothPS(OmegaM, OmegaL, OmegaB, h_, sig8, n, R);
	 
	// Cecile - replace new computation by reading the theoretical spectrum without oscillations
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
	
	// set default ka range
	//minka_=0.03; maxka_=0.07; nka_=1000;
	
	minka_=0.035; maxka_=0.05; 
	nka_=1000;
	
	bestfit_=-10; // uninitialised
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

/*
// Calculate SMOOTH fiducial power spectrum
void FitBAOScale::ComputeSmoothPS(double OmegaM, double OmegaL, double OmegaB, 
                                      double h, double sig8, double n, double R)
{
	
	cout << "     Computing smooth fiducial power spectrum ..."<<endl;
	InitialPowerLaw Pkinit(n);
	TransferEH tf(h,OmegaM-OmegaB,OmegaB,T_CMB_K,false);
	//	tf.SetNoOscEnv(1); // want smooth power spectrum
	tf.SetNoOscEnv(2); // want smooth power spectrum Cecile, to do like cmvginit3d
	GrowthEH growth(OmegaM, OmegaL);
 	double growth_at_z = growth(zref_);
 	cout << "     Growth factor at z="<< zref_ <<" = "<< growth_at_z <<endl;
 	PkSpecCalc pkz(Pkinit,tf,growth,zref_);
 	//PkSpectrum0 pk0(Pkinit,tf);
 	//PkSpectrumZ pkz(pk0,growth,zref_);
	pkz.SetZ(0.);
 	cout <<endl<<"     Compute variance for top-hat R="<< R <<" (sigma"<< R;
 	cout <<") at z="<< pkz.GetZ() <<endl;
 	VarianceSpectrum varpk_th(pkz,R,VarianceSpectrum::TOPHAT);
	double kmin=1e-5,kmax=1000.;
	int npt = 10000;
	//double lkmin=log10(kmin), lkmax=log10(kmax);
	double eps=1.e-3;
 	double kfind_th = varpk_th.FindMaximum(kmin,kmax,eps);
 	double pkmax_th = varpk_th(kfind_th);
 	cout <<"     kfind_th = "<< kfind_th <<" ("<< log10(kfind_th);
 	cout <<"), integrand="<< pkmax_th <<endl;
 	double k1=kmin, k2=kmax;
 	int rc = varpk_th.FindLimits(pkmax_th/1.e4,k1,k2,eps);
 	cout <<"     limit_th: rc="<< rc <<" : "<< k1 <<" ("<< log10(k1);
 	cout <<") , "<< k2 <<" ("<< log10(k2) <<")"<<endl;
	double ldlk = (log10(k2)-log10(k1))/npt;
	varpk_th.SetInteg(0.01,ldlk,-1.,4);
	double sr2 = varpk_th.Variance(k1,k2);
	cout <<"     varpk_th="<< sr2 <<"  ->  sigma="<< sqrt(sr2) <<endl;
	double normpkz = sig8*sig8/sr2;
	pkz.SetScale(normpkz);
	cout <<"     Spectrum normalisation = "<< pkz.GetScale() <<endl;
	pkz.SetZ(zref_);
	cout <<endl<<"     Compute variance for Pk at z="<< pkz.GetZ() <<endl;
 	VarianceSpectrum varpk_int(pkz,R,VarianceSpectrum::NOFILTER);
	double kfind_int = varpk_int.FindMaximum(kmin,kmax,eps);
 	double pkmax_int = varpk_int(kfind_int);
 	cout <<"     kfind_int = "<< kfind_int <<" ("<< log10(kfind_int);
 	cout <<"), integrand="<< pkmax_int <<endl;
 	double k1int=kmin, k2int=kmax;
 	int rcint = varpk_int.FindLimits(pkmax_int/1.e4,k1int,k2int,eps);
 	cout <<"     limit_int: rc="<< rcint <<" : "<< k1int <<" ("<< log10(k1int);
 	cout <<") , "<< k2int <<" ("<< log10(k2int) <<")"<<endl;
	double ldlkint = (log10(k2int)-log10(k1int))/npt;
	varpk_int.SetInteg(0.01,ldlkint,-1.,4);
	double sr2int = varpk_int.Variance(k1int,k2int);
	cout <<"     varpk_int="<< sr2int <<"  ->  sigma="<< sqrt(sr2int) <<endl;
	cout << "     ... Finished computing smooth fiducial power spectrum "<<endl;
	cout << endl;

	for (int i=0; i<kobs_.Size(); i++)
		Pref_(i) = pkz(kobs_(i));
};
*/

// Compute the chi-square as a function of ka
void FitBAOScale::ComputeChisq(double maxk)
{

	double dka = (maxka_-minka_)/(nka_-1);
	double amp = 2.5;	// amplitude is fixed
	kavals_.SetSize(nka_);
	Chisq_.SetSize(nka_);
	cout << "    Computing chisq from ka = "<< minka_ <<" to ka = ";
	cout << maxka_ <<" in steps of "<< dka <<endl;
	

	//DecaySineFunc decs(amp, h_);
	DecaySineFunc decs(amp, 1);
	cout << "    Check decaying sine function "<<h_ <<endl;
	//double ka =0.04;
	//for (int ik=0; ik<kobs_.Size(); ik++)
	//		{
	//		double pred=decs(kobs_(ik),ka);
	//		cout << kobs_(ik) <<"   "<< pred<<endl;
	//		}
	decs.PrintParas();

	for (int ika=0; ika<nka_; ika++) {
	  double chisq=0;
	  kavals_(ika) = minka_ + ika*dka;
	  //kavals_(ika) /= h_;  	//k unit should be h/Mpc
	  //cout << " ka = "<< kavals_(ika) << endl;
	  
	  int indx = 0;
	  for (int ik=0; ik<kobs_.Size(); ik++) {
	    
	    // decaying sinusoid at ka 
	    double pred = decs(kobs_(ik),kavals_(ika));
	    double diff= Pratio_(ik)-pred;
	    
	    if (sig_(ik)==0)
	      throw ParmError("Error is zero!");
	    
            // only use power spectra values below some max k
	    
	    if (kobs_(ik)<maxk && kobs_(ik)>0.02){
	      chisq += diff*diff/(sig_(ik)*sig_(ik));
	    
	      //if (isnan(pow(diff,2.)/(sig_(ik)*sig_(ik))))
	      /*if( kavals_(ika)>= 0.0350 &&  kavals_(ika)< 0.0351){
	       	  // cout << kavals_(ika) << " k " << kobs_(ik) << " sig_(ik)="<<sig_(ik)<<", Pratio="<<Pratio_(ik)<<", Pred "<<pred;
		  // cout << ", diff "<< diff << endl;
		   cout << ik << " " << kobs_(ik) << " " << pred << " " << Pratio_(ik) << "  "<< sig_(ik) << " " << chisq << endl;
	      }*/
	      indx ++;
	    }
	  }
	  
	  Chisq_(ika) = chisq/(indx-1);
	  //if( kavals_(ika)> 0.0350 &&  kavals_(ika)< 0.0351)
	   //cout << "ka = " <<  kavals_(ika) << " chisq = "<< chisq  << endl;
	}
};

// Compute the chi-square as a function of ka and A (Adeline modif)
void FitBAOScale::ComputeChisq_Avar(double maxk)
{
   double dka = (maxka_-minka_)/(nka_-1);
   double kaVals_temp = minka_;
   		
   Amin_ = 0;
   Amax_ = 3;
   nA_ = 1000;
   //nA_ = 100; //modif Marion
   double dA = (Amax_ - Amin_)/nA_;
   double amp=Amin_;	
   
   double vecSize = nka_ * nA_;
   
   kavals_.SetSize(vecSize);
   kavals_short_.SetSize(nka_);
   ampvals_.SetSize(vecSize);
   ampvals_short_.SetSize(nA_);
   Chisq_.SetSize(vecSize);
   cout << "    Computing chisq from ka = "<< minka_ <<" to ka = ";
   cout << maxka_ <<" in steps of "<< dka << " (nka="<< nka_ <<")" <<endl;
   cout << "    A is variable : min = " << Amin_ << " to max = " << Amax_ << " in step of "<< dA << endl;
   cout << "    Total :  compute "<< vecSize << " values of chisq "<< endl;  

   int ndim = 2;
   sa_size_t dim[ndim];
   dim[0] = nA_;
   dim[1] = nka_;
   chisq_avar_.SetSize(ndim, dim);

   int indx = 0;
   for(int iA=0; iA<nA_; iA ++){
        amp = Amin_ + iA*dA;
	ampvals_short_(iA) = amp;
	
	//DecaySineFunc decs(amp, h_);
	DecaySineFunc decs(amp, 1);
	//cout << indx << "    Check decaying sine function "<<h_ << " amp = "<< amp <<endl;
	
	//double ka =0.04;
	//for (int ik=0; ik<kobs_.Size(); ik++)
	//		{
	//		double pred=decs(kobs_(ik),ka);
	//		cout << kobs_(ik) <<"   "<< pred<<endl;
	//		}
	

//decs.PrintParas();

	for (int ika=0; ika<nka_; ika++) {
	  kaVals_temp = minka_ + ika*dka;
	 // if(iA==0)
	     kavals_short_(ika) = kaVals_temp;
	  kavals_(indx) = kaVals_temp;
	  ampvals_(indx) = amp;
	  
	  double chisq=0;
	  int kindx = 0;

	  for (int ik=0; ik<kobs_.Size(); ik++) {
	    
	    // decaying sinusoid at ka 
	    double pred = decs(kobs_(ik),kaVals_temp);
	    double diff= Pratio_(ik)-pred;
	    
	    if (sig_(ik)==0)
	      throw ParmError("Error is zero!");
	    
            // only use power spectra values below some max k
	    if (kobs_(ik)<maxk && kobs_(ik)>0.02){
	      chisq += diff*diff/(sig_(ik)*sig_(ik));
	    
	      //if (isnan(pow(diff,2.)/(sig_(ik)*sig_(ik))))
	      /*if( kavals_(ika)> 0.0449 &&  kavals_(ika)< 0.0451 && kobs_(ik)>0.06 && kobs_(ik)<0.07){
	       	   cout << kavals_(ika) << " k " << kobs_(ik) << " sig_(ik)="<<sig_(ik)<<", Pratio="<<Pratio_(ik)<<", Pred "<<pred;
		   cout << ", diff "<< diff << endl;
	      }*/
	      kindx ++;
	    }
	  }
	  
	  Chisq_(indx) = chisq/(kindx-2);
	  chisq_avar_(iA, ika) = chisq/(kindx-2);
	  //cout << " dans computechisq_avar : " << chisq_avar_(iA, ika) << endl;

	  //if( kavals_(ika)> 0.041 &&  kavals_(ika)< 0.042)
	 //    cout << "ka = " <<  kavals_(ika) << " chisq = "<< chisq  << endl;
	  indx ++;
	}
   }
   cout << endl << " test amp "<< endl;
   indx = 0;
   for(int i=0; i<nA_; i++){
     	for(int j=0; j<nka_; j++) {
   	    //cout <<"  A "<< ampvals_short_(i) << "  ka "<< kavals_short_(j) <<" chi "<< chisq_avar_(i,j) <<endl;
	    //cout <<"  A "<< ampvals_(indx) << "  ka "<< kavals_(indx) <<" chi "<< chisq_avar_(i,j) <<endl;
	    indx ++;
	}
   }
   cout << "   check compute chisq :   indx = "<< indx << " Avals.size = "<< ampvals_.Size() << " ka.size = "<< kavals_.Size();
   cout << " chisq.size = " << Chisq_.Size() << endl;
   cout << "   chisq_avar_: x length="<<chisq_avar_.SizeX() << " y length="<< chisq_avar_.SizeY() << " shoulb be equal to (";
   cout << nA_ << ", "<< nka_ <<")"<<endl;
 
};


// write chi-square to a file
void FitBAOScale::WriteChisq(string outfile)
{
   if (bestfit_<0)
	throw ParmError("Have not calculated best fit ka yet!");

   ifstream inp;
   ofstream outp;
   inp.open(outfile.c_str(), ifstream::in);
   inp.close();
   if(inp.fail()) {
	inp.clear(ios::failbit);
	cout << "    Writing chisq to file ..." << outfile.c_str() << endl;
	outp.open(outfile.c_str(), ofstream::out);
		
	outp << "Redshift of power spectrum = "<< zref_ <<", best fit ka = ";
	outp << bestfit_ <<", "<< nsig_ <<"-sigma errors: + "<< errup_;
	outp <<" - "<< errdown_ << " with h = " << h_ << endl;
		
	if(ampvals_.Size() != 0){ //Adeline modif
		outp << "amplitude is variable, best fit A = "<< bestfit_amp_ <<endl;
		outp << "first column : ka(nka = "<<chisq_avar_.SizeY()<<"); first line : amp(nA = "<<chisq_avar_.SizeX()<<")" <<endl;
			
		//write first line (A values)
		outp << "0 \t ";
		for(int j=0; j<ampvals_short_.Size(); j++)
		 	outp << ampvals_short_(j) << " \t" ;
		outp << endl;
			
		for(int j=0; j<chisq_avar_.SizeY(); j++){// loop over ka values
			outp << kavals_short_(j) << " \t";
			for (int i=0;i<chisq_avar_.SizeX();i++){// loop over A values
				outp << chisq_avar_(i,j) << "\t ";
			}
			outp << endl;
		}	
	}
	else{	
		for (int i=0;i<Chisq_.Size();i++)// loop over ka values
			outp << kavals_(i) << "   "<< Chisq_(i)<< endl;
	}
	outp.close();
   } // end of write
   else
	cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

};

// calculate best-fit and errors
void FitBAOScale::BestfitStdDev(double& bestfit, double& siglow, double& sighigh, int nsig)
{
	nsig_ = nsig;
	if (nsig_>3)
		throw ParmError("Can't compute above 3-sigma");

	vector<double> clevels;
	clevels.push_back(0.683); clevels.push_back(0.954); clevels.push_back(0.997);

	double clevel = clevels[nsig_-1];
	cout <<"    Confidence level = "<< clevel <<endl;
	ChisqStats chisqstat(kavals_, Chisq_);
	bestfit = chisqstat.BestFit();
	chisqstat.ErrSig(siglow,sighigh,clevel,100);
	
	bestfit_ = bestfit;
	errup_   = sighigh - bestfit_;
	errdown_ = bestfit_ - siglow;
	
	bestfit_amp_ = bestfit;
};

// calculate best-fit and errors
void FitBAOScale::BestfitStdDev_Avar(double& bestfit_ka, double& bestfit_A, double& siglow, double& sighigh, int nsig)
{
   nsig_ = nsig;
   if (nsig_>3)
	throw ParmError("Can't compute above 3-sigma");

   vector<double> clevels;
   clevels.push_back(0.683); clevels.push_back(0.954); clevels.push_back(0.997);

   double clevel = clevels[nsig_-1];
   cout <<"    Confidence level = "<< clevel <<endl;
   ChisqStats chisqstat(kavals_, ampvals_, Chisq_, chisq_avar_, 2);
   bestfit_ka = chisqstat.BestFit_Avar(bestfit_A);
   //chisqstat.ErrSig(siglow,sighigh,clevel,100);

   bestfit_ = bestfit_ka;
   errup_   = sighigh - bestfit_;
   errdown_ = bestfit_ - siglow;
	
   cout << endl << " marginalize chisq "<< endl;
   double dA = (Amax_ - Amin_)/nA_;
   double dka =  (maxka_ - minka_)/nA_;
   
   ChisqStats chisqstat_forMarg(kavals_short_, ampvals_short_, Chisq_, chisq_avar_, 2);

   double step[] = {dka, dA};
   cout <<" step "<< dka << "  "<< dA << "  "<< step[0] << "  "<< step[1] << endl;
   //TArray<r_8> margChisq_;

    int ndim = 2;
   sa_size_t dim[ndim];
   dim[0] = 2;
   if(nka_>=nA_)
   	dim[1] = nka_;
   else
   	dim[1] = nA_;

   //dim[0] = dim[1] = nA_*nka_; //modif Marion

   margChisq_.SetSize(ndim, dim);
   cout << " margChisq_ dim:  X "<< margChisq_.SizeX() << " Y: "<< margChisq_.SizeY() << endl;
   
   chisqstat.GetMarg(step, margChisq_);
   
   //chisqstat.ErrSig(siglow,sighigh,clevel,100);

   bestfit_ = bestfit_ka;
   errup_   = sighigh - bestfit_;
   errdown_ = bestfit_ - siglow;

   //added by Marion
   bestfit_amp_ = bestfit_A;

   cout << "bestfit_amp_ : " << bestfit_amp_ << endl;
};

void FitBAOScale::WriteChisqMarg(string outfile)
{ 
   if (bestfit_<0)
	throw ParmError("Have not calculated best fit ka yet!");
   if(margChisq_.SizeX()<=0 || margChisq_.SizeY()<=0 )
   	throw ParmError("Have not marginalized chisq yet!");
   
   ifstream inp;
   ofstream outp;
   inp.open(outfile.c_str(), ifstream::in);
   inp.close();
   double temp_ka, temp_amp;
   if(inp.fail()) {
	inp.clear(ios::failbit);
	cout << "    Writing chisq to file ..." << outfile.c_str() << endl;
	outp.open(outfile.c_str(), ofstream::out);
		
	outp << "Redshift of power spectrum = "<< zref_ <<", best fit ka = ";
	outp << bestfit_ <<", "<< nsig_ <<"-sigma errors: + "<< errup_;
	outp <<" - "<< errdown_ << " with h = " << h_ << endl;
		
	if(ampvals_.Size() != 0){ //Adeline modif
	  outp << "amplitude is variable, best fit A = "<< bestfit_amp_ <<endl;
	 
		outp << "ka, A, chisq(ka), chisq(A)" <<endl;
			
		for(int j=0; j<margChisq_.SizeY(); j++){
			temp_ka = kavals_short_(j);
			temp_amp = ampvals_short_(j);
			if(j>kavals_short_.Size() )
				temp_ka = 0;
			if(j>ampvals_short_.Size() )
				temp_amp = 0;
	
			outp <<temp_ka << "\t"<< temp_amp <<" \t";
			for(int i=0; i<margChisq_.SizeX(); i++){
				outp << margChisq_(i,j) << " \t";
			}
			outp << endl;
		}	
	}
	else
		cout << " ERRORS "<< endl;
	
	outp.close();
   } // end of write
   else
	cout << "Error...file """ << outfile.c_str() << """ exists" << endl;   
};


// analytical approx for sample variance error
void FitBAOScale::CalcSigSampVar(double VolCat)
// Analytical approximation of the error on the power spectrum
// due to the sample variance
// Won't be accurate for photo-z power spectrum I expect
{

    sig_.SetSize(kobs_.Size());
    double dk = kobs_(2)-kobs_(1);

    for (int i=0; i<kobs_.Size(); i++) 
        sig_(i) = sqrt(2*pow(2*PI,3)/VolCat*(1./(4*PI*pow(kobs_(i),2)*dk)))*Pratio_(i);

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
		
		cout << "    Redshift of power spectrum = "<< zref_ <<", best fit ka = ";
		cout << bestfit_ <<", "<< nsig_ <<"-sigma errors: + "<< errup_ <<" - ";
		cout << errdown_ <<endl;
		
		outp << "z_ps : best fit ka : n-sigma of errors : +error : -error : hubble parameter h"<<endl;
		
		outp << zref_ <<"   "<< bestfit_ <<"   "<< nsig_ <<"   "<< errup_;
		outp <<"   "<< errdown_ << "   " << h_ << endl;
		
		outp.close();
		} // end of write
    else
	    cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

};


// write reference power spectrum to a file
void FitBAOScale::WriteRefPS(string outfile)
{

	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
		{
		inp.clear(ios::failbit);
		cout << "    Writing reference power spectrum to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
		for (int i=0;i<Pref_.Size();i++)// loop over k values
			outp << kobs_(i) << "   "<< Pref_(i)<< endl;
		outp.close();
		} // end of write
		else
			cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

}

void FitBAOScale::WriteAncillaryInfo(string outfile)
{

	if (bestfit_<0)
		throw ParmError("Have not calculated best fit ka yet!");

	double amp=2.4;
	DecaySineFunc decs(amp,h_);

	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(ios::failbit);
		cout << "    Writing reference power spectrum and decaying sinosoid";
		cout << " to file ..." << outfile.c_str() << endl;
		
		outp.open(outfile.c_str(), ofstream::out);
		
		for (int i=0; i<Pref_.Size(); i++) {// loop over k values
			double pred = decs(kobs_(i),bestfit_);
			outp << kobs_(i) << "   "<< Pref_(i)<< "   "<< pred <<endl;
			}
		outp.close();
		} // end of write
    else
        cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

};

void FitBAOScale::WriteWiggle(string outfile)
{

	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(ios::failbit);
		cout << "    Writing wiggle power spectrum p(k)/psmooth(k) and sigma(Pwiggle)";
		cout << " to file ..." << outfile.c_str() << endl;
		
		outp.open(outfile.c_str(), ofstream::out);
		
		for (int i=0; i<kobs_.Size(); i++){
		    //Pratio_(i) = Pobs_(i)/Pref_(i); 
		    //sig_(i) = sig_(i) / Pref_(i); 
		    outp << kobs_(i) << "   "<< Pratio_(i)<< "   "<< sig_(i) <<endl;
		}
		
		outp.close();
		} // end of write
    else
        cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

};



