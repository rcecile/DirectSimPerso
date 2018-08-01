/* ----
   Project   LSST/BAO/PhotoZ
   Classes to generate galaxie Absolute Magnitude and types starting  
   from Luminosity Functions (LF's)  
   F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud 
   April - July 2017                                          -------  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "multitypzlf.h"
#include "pexceptions.h"
#include "array.h"   // -- for reading the ASCII input file as a TMatrix 
#include "randfmt.h"
#include "ctimer.h"

using namespace std;
using namespace SOPHYA;

//------------------------------------------------------------------------
//-----------  class MultiType_Z_LF --------------------------------------
//------------------------------------------------------------------------

//--- static variables and methods
//-- redshift value at which number of galaxies reaches zero 
static double def_zero_gal_redshift = 20.;
//-- redshift range and number of number of redshifts for which we compute the Integrale[LF](mag) for random drawing of magnitudes
static double def_zmin = 0.;
static double def_zmax = 4.;
static int def_nz = 20;
//-- Schechter integration Gauss-Legendre order, and incremental integration order
static int def_glorder = 100;
static int def_inc_glorder = 5;
//--  do CPU timing or NOT
static bool fgdo_CPUTiming = false;

// define the redshift at which number of galaxies reaches zero
void MultiType_Z_LF::set_def_zero_gal_redshift(double z)
{
  if ((z>0.) && (z<1100.))  def_zero_gal_redshift=z;
  return;
}
// define the redshift range and number of bins for magnitude distribution for all galaxies
void MultiType_Z_LF::set_def_redshift_range(double zmin, double zmax, double nz)
{
  if ((zmin>=0.) && (zmax>zmin) && (nz>0)) {
    def_zmin=zmin;  def_zmax=zmax;  def_nz=nz;
  }
  return;  
}
// define the Gauss-Legendre quadrature order Schechter LF integration from magMin to magMax and for incremental integration
void  MultiType_Z_LF::set_def_schechter_integ_glorder(int glorder, int incglorder)
{
  if (glorder>0)  def_glorder=glorder;
  if (incglorder>0)  def_inc_glorder=incglorder;
  return;
}

void  MultiType_Z_LF::activateCPUTiming(bool fgtim)
{
  fgdo_CPUTiming=fgtim;
  return;
}

/* --Methode-- */
MultiType_Z_LF::MultiType_Z_LF(string const & lfparamsfilename, double magMin, double magMax, size_t nbmag, bool useall)
  : magMin_(magMin), magMax_(magMax), nb_mag_pts_(nbmag), useSchechall_(useall)
{
  Timer* ptm;
  if (fgdo_CPUTiming)  ptm = new Timer("MultiType_Z_LF");
  
  TMatrix<double> lfp;
  ifstream is(lfparamsfilename.c_str());
  if (!is) {cerr << " Err: multitypzlf.cc: " << lfparamsfilename
		 << " does not exist. Exiting ..." << endl; exit(0); }
  sa_size_t nr,nc;
  lfp.ReadASCII (is, nr, nc); //  char clm='#', const char *sep=" \t")
  cout << " MultiType_Z_LF[1]: read from file"<<lfparamsfilename<<" -> nb of lines="<<nr<<" number of columns=" << nc << endl;
  //DBG  lfp.Print(cout);
  double z2last=0.;
  for(sa_size_t j=0; j<lfp.SizeY(); j++)  {
    double z1 = lfp(j,0);
    double z2 = lfp(j,1);
    cout << z1 << "  "  << z2 << endl;
    if ((j>0)&&(fabs(z1-z2last)>1.e-9))
      throw ParmError("MultiType_Z_LF::MultiType_Z_LF: redshift interval error in file");
    pair<double,double> z12;   z12.first=z1;  z12.second=z2;
    zrange.push_back(z12);
    z2last=z2;
    
    sa_size_t off=2;
    MySchechter sfell(lfp(j,off), lfp(j,off+1), lfp(j,off+2));
    cout << " LFP " << lfp(j,off) << "  " << lfp(j,off+1)  << "  " << lfp(j,off+2) << endl;
    vshEll.push_back(sfell);
    off+=3;
    MySchechter sfsp(lfp(j,off), lfp(j,off+1), lfp(j,off+2));
    cout << " LFP " << lfp(j,off) << "  " << lfp(j,off+1)  << "  " << lfp(j,off+2) << endl;
    vshSp.push_back(sfsp);
    off+=3;
    MySchechter sfsb(lfp(j,off), lfp(j,off+1), lfp(j,off+2));
    cout << " LFP " << lfp(j,off) << "  " << lfp(j,off+1)  << "  " << lfp(j,off+2) << endl;
    vshSB.push_back(sfsb);
    off+=3;
    MySchechter sfall(lfp(j,off), lfp(j,off+1), lfp(j,off+2));
    cout << " LFP " << lfp(j,off) << "  " << lfp(j,off+1)  << "  " << lfp(j,off+2) << endl;
    vshAll.push_back(sfall);
  }

  // checking if the first range starts at zero 
  first_z_is_zero=false;
  if (fabs(zrange[0].first)<1.e-12)  {
    zrange[0].first=0.;
    first_z_is_zero=true;
  }

  zero_gal_redshift_=def_zero_gal_redshift;
  cout << " MultiType_Z_LF[2]: Computing Integral of LF's from magMin="<<magMin_<<" to magMax="<<magMax_<<" \n"
       << "  with linear interpolation as a function of redshift , with ngal=0 for redshift z="<<zero_gal_redshift_<<endl; 

  size_t vsz = (first_z_is_zero ? zrange.size()+3 : zrange.size()+4);
  size_t voff = (first_z_is_zero ? 1 : 2);
  
  vector<double> leszm(vsz); // z at the middle of the interval
  vector<double> lesilfEll(vsz); // integral of LF at a given z interval
  vector<double> lesilfSp(vsz);
  vector<double> lesilfSB(vsz);
  vector<double> lesilfAll(vsz);

  if ( first_z_is_zero ) 
    leszm[0]=zrange[0].first;
  else {
    leszm[0]=0.;
    leszm[1]=zrange[0].first;
  }
  for(size_t k=0; k<zrange.size(); k++)  {
    leszm[k+voff]=(zrange[k].first+zrange[k].second)*0.5;
    MySchechter sf=vshEll[k];
    lesilfEll[k+voff]=sf.getIntegral(magMin_,magMax_,def_glorder);
    sf=vshSp[k];
    lesilfSp[k+voff]=sf.getIntegral(magMin_,magMax_,def_glorder);
    sf=vshSB[k];
    lesilfSB[k+voff]=sf.getIntegral(magMin_,magMax_,def_glorder);
    sf=vshAll[k];
    lesilfAll[k+voff]=sf.getIntegral(magMin_,magMax_,def_glorder);
  }
  
  leszm[zrange.size()+voff] = zrange[zrange.size()-1].second;
  leszm[zrange.size()+voff+1] = zero_gal_redshift_;   // maaximum z at which number of galaxies goes to zero

  //DBG  for(size_t ii=0; ii<leszm.size(); ii++) cout<<" ---DBG-- ii="<<ii<<" leszm[ii]="<<leszm[ii]<<endl;

  if ( first_z_is_zero ) {
    lesilfEll[0]=lesilfEll[voff];
    lesilfSp[0]=lesilfSp[voff];
    lesilfSB[0]=lesilfSB[voff];
    lesilfAll[0]=lesilfAll[voff];
  }
  else {
    lesilfEll[0]=lesilfEll[1]=lesilfEll[voff];
    lesilfSp[0]=lesilfSp[1]=lesilfSp[voff];
    lesilfSB[0]=lesilfSB[1]=lesilfSB[voff];
    lesilfAll[0]=lesilfAll[1]=lesilfAll[voff];
  }
  
  lesilfEll[zrange.size()+voff]=lesilfEll[zrange.size()+voff-1];
  lesilfSp[zrange.size()+voff]=lesilfSp[zrange.size()+voff-1];
  lesilfSB[zrange.size()+voff]=lesilfSB[zrange.size()+voff-1];
  lesilfAll[zrange.size()+voff]=lesilfAll[zrange.size()+voff-1];

  lesilfEll[zrange.size()+voff+1]=0.;
  lesilfSp[zrange.size()+voff+1]=0.;
  lesilfSB[zrange.size()+voff+1]=0.;
  lesilfAll[zrange.size()+voff+1]=0.;
  
  ILF_of_z_Ell.DefinePoints(leszm, lesilfEll);
  ILF_of_z_Sp.DefinePoints(leszm, lesilfSp);
  ILF_of_z_SB.DefinePoints(leszm, lesilfSB);
  ILF_of_z_All.DefinePoints(leszm, lesilfAll);

  zmin_=def_zmin;
  zmax_=def_zmax;
  nz_=def_nz;
  dz_=(zmax_-zmin_)/(double)nz_;

  //  --- CPU timing 
  if (fgdo_CPUTiming)  ptm->Split("Done [1]+[2]");
			   
  cout << " MultiType_Z_LF[3]: Computing Sum_[Ell,Spiral,StarBurst]  of Integral_LF[magMin="
       <<magMin_<<","<<"magMax="<<magMax_<<"]  \n for "<<nz_+1<<" redshifts "<<zmin_<<" <=z<= "<<zmax_<<endl;
  vector<double> lesz((size_t)(nz_+1));
  vector<double> lessum((size_t)(nz_+1));
  for(size_t k=0; k<lesz.size(); k++) {
    double redshift=zmin_+(double)k*dz_;
    double sumILF=ILF_of_z_Ell(redshift)+ILF_of_z_Sp(redshift)+ILF_of_z_SB(redshift);
    lesz[k]=redshift;
    lessum[k]=sumILF;
  }
  ILF_of_z_Sum_EllSpSB.DefinePoints(lesz, lessum);
  
  //  --- CPU timing 
  if (fgdo_CPUTiming)  ptm->Split("Done [3]");
			   
  vector<double> lesmag(nb_mag_pts_+1);
  vector<double> lesILF(nb_mag_pts_+1);
  double dmag=(magMax_-magMin_)/(double)nb_mag_pts_;
  lesmag[0]=magMin_;
  double curmag=magMin_;
  for(size_t k=0; k<lesmag.size()-1; k++) {
    lesmag[k]=curmag;  curmag+=dmag;
  }
  lesmag[lesmag.size()-1]=magMax_;

  if (useSchechall_) {  // We use Schechter parameters for All galaxies
    cout << " MultiType_Z_LF[4]: Computing Integral_LF as a function of magnitude from All galaxies Schechter \n"
	 << " parameters for the corresponding " << zrange.size() << " redshift intervals ..."<<endl;
    v_mag_f_ILF.resize(zrange.size(),NULL);
    for(size_t k=0; k<zrange.size(); k++)  {
      SLinInterp1D * mag_f_ILF = new SLinInterp1D;
      double curILF = 0.;
      lesILF[0]=curILF;
      MySchechter sf=vshAll[k];
      for(size_t i=1; i<lesmag.size(); i++) {
	curILF += sf.getIntegral(lesmag[i-1],lesmag[i],def_inc_glorder);
	lesILF[i]=curILF;
      }
      for(size_t i=1; i<lesILF.size(); i++)  lesILF[i] /= curILF;
      lesILF[lesILF.size()-1]=1.;
      mag_f_ILF->DefinePoints(lesILF, lesmag);
      v_mag_f_ILF[k]=mag_f_ILF;
    }
  }
  else {  // We use Sum[Ell+Spiral+StarBurst] for computing magnitude distribution for  All galaxies
    cout << " MultiType_Z_LF[4]: Computing Integral_LF as a function of magnitude from Sum[Ell+Spiral+StarBurst] Schechters \n"
	 << " for "<<nz_+1<<" redshifts values " <<zmin_<<" <=z<= "<<zmax_<<endl;
    v_mag_f_ILF.resize((size_t)nz_,NULL);

    lesz.resize((size_t)nz_);
    for(size_t k=0; k<lesz.size(); k++) {
      double redshift=zmin_+((double)k+0.5)*dz_;
      SLinInterp1D * mag_f_ILF = new SLinInterp1D;
      double curILF = 0.;
      lesILF[0]=curILF;
      MySchechter sfE=getLF_Elliptical(redshift);
      MySchechter sfSp=getLF_Spiral(redshift);
      MySchechter sfSB=getLF_StarBurst(redshift);
      for(size_t i=1; i<lesmag.size(); i++) {
	curILF += sfE.getIntegral(lesmag[i-1],lesmag[i],def_inc_glorder);
	curILF += sfSp.getIntegral(lesmag[i-1],lesmag[i],def_inc_glorder);
	curILF += sfSB.getIntegral(lesmag[i-1],lesmag[i],def_inc_glorder);
	lesILF[i]=curILF;
      }
      for(size_t i=1; i<lesILF.size(); i++)  lesILF[i] /= curILF;
      lesILF[lesILF.size()-1]=1.;
      mag_f_ILF->DefinePoints(lesILF, lesmag);
      v_mag_f_ILF[k]=mag_f_ILF;
    }
  }
  //  --- CPU timing 
  if (fgdo_CPUTiming)  ptm->Split("Done [4]");

  //  Creating / defining default random generator 
  def_rgp_ = new FMTRandGen();
  rgp_ = def_rgp_;
  
  //  --- CPU timing 
  if (fgdo_CPUTiming)  delete ptm;
}

/* --Methode-- */
MultiType_Z_LF::~MultiType_Z_LF()
{
  delete def_rgp_ ;
  for (size_t i=0; i<v_mag_f_ILF.size(); i++)  delete v_mag_f_ILF[i];
}

/* --Methode-- */
size_t MultiType_Z_LF::find_redshift_bin(double z)  const 
{
  if (z < 0.) 
    throw ParmError("MultiType_Z_LF::find_redshift_bin()/ERROR negative redshift !");

  size_t rzb=0;
  for(size_t k=0; k<zrange.size(); k++) {
    if (z < zrange[k].first)  break;
    rzb++;
  }
  if (rzb>0)  rzb--;
  if (rzb >= zrange.size())  rzb=zrange.size()-1;
  return rzb;
}


/* --Methode-- */
MySchechter MultiType_Z_LF::getLF_Elliptical(double z)  const 
{
  size_t kbin = find_redshift_bin(z);
  MySchechter rsf = vshEll[kbin];
  // scale phistar
  rsf.scalePhistar( ILF_of_z_Ell(z)/ILF_of_z_Ell(0.5*(zrange[kbin].first+zrange[kbin].second)) );
  return rsf;
}

/* --Methode-- */
MySchechter MultiType_Z_LF::getLF_Spiral(double z)  const 
{
  size_t kbin = find_redshift_bin(z);
  MySchechter rsf = vshSp[kbin];
  // scale phistar
  rsf.scalePhistar( ILF_of_z_Sp(z)/ILF_of_z_Sp(0.5*(zrange[kbin].first+zrange[kbin].second)) );
  return rsf;
}

/* --Methode-- */
MySchechter MultiType_Z_LF::getLF_StarBurst(double z)   const 
{
  size_t kbin = find_redshift_bin(z);
  MySchechter rsf = vshSB[kbin];
  // scale phistar
  rsf.scalePhistar( ILF_of_z_SB(z)/ILF_of_z_SB(0.5*(zrange[kbin].first+zrange[kbin].second)) );
  return rsf;
}

/* --Methode-- */
MySchechter MultiType_Z_LF::getLF_All(double z)  const 
{
  size_t kbin = find_redshift_bin(z);
  MySchechter rsf = vshAll[kbin];
  // scale phistar
  rsf.scalePhistar( ILF_of_z_All(z)/ILF_of_z_All(0.5*(zrange[kbin].first+zrange[kbin].second)) );
  return rsf;
}

/* --Methode-- */
double MultiType_Z_LF::getGalaxyNumberDensity(double z)  const 
{
  if (useSchechall_)  return ILF_of_z_All.YInterp(z);
  return ILF_of_z_Sum_EllSpSB.YInterp(z);
}

/* --Methode-- */
void MultiType_Z_LF::P_getTypeMagnitude(double z, double& mag, bool fgmag, double& type, bool fgtype)  const 
{
  // DEFINE:
  // type 0  = El_cww_fix2.txt
  // type 10 = Sbc_cww_fix2.txt
  // type 20 = Scd_cww_fix2.txt
  // type 30 = Im_cww_fix2.txt
  // type 40 = SB3_kin_fix2.txt
  // type 50 = SB2_kin_fix2.txt
  // Choice is same as Dahlen et al arXiv:0710.5532

  //USED FOR planck_BAO_LF 
  // int a1=0,  b1=5;    //early (1)
  //  int a2=5,  b2=35;   //late  (2)
  //  int a3=35, b3=51;   //starburst (3)
  // to be used with Dahlen LF
  // separation at 35 : irregular are considered as late galaxies.  JSR/CR 

  int a1=0,  b1=5;    //early (1)
  int a2=5,  b2=15;   //late  (2)
  int a3=15, b3=51;   //starburst (3)
  // to be used with LF Zucca :https://www.aanda.org/articles/aa/pdf/2009/48/aa12665-09.pdf
  // "Galaxies were then classified into four types, 
  // corresponding to colors of E/S0 (type 1), early spirals (type 2), 
  // late spirals (type 3), and irregular and starburst galaxies (type 4)."
  // LF table for type 1, 2, 3+4

  size_t kbin = 0;
  size_t kbinm = 0;
  if (fgtype)  kbin = find_redshift_bin(z);
  if (fgmag) {  // draw magnitude at redshift z 
    if (useSchechall_)  {
      kbinm = (fgtype ? kbin : find_redshift_bin(z));
    }
    else {
      if (z < 0.) 
	throw ParmError("MultiType_Z_LF::getMag(z)/ERROR negative redshift !");
      kbinm = int((z-zmin_)/dz_);
      if (kbinm >= v_mag_f_ILF.size())  kbinm=v_mag_f_ILF.size()-1;
      //      cout << " size " <<v_mag_f_ILF.size()<<"  "  << kbinm <<endl;
    }
    double rnd=rgp_->Flat01();
    mag = v_mag_f_ILF[kbinm]->YInterp(rnd);
  }
  if (fgtype) {  // draw type for magnitude mag at redshift z 
    int typ123=0;
    // Get rescaled Schechter functions for the three types 
    MySchechter rsfE = vshEll[kbin];
    MySchechter rsfSp = vshSp[kbin];
    MySchechter rsfSB = vshSB[kbin];
    rsfE.scalePhistar( ILF_of_z_Ell(z)/ILF_of_z_Ell(0.5*(zrange[kbin].first+zrange[kbin].second)) );
    rsfSp.scalePhistar(ILF_of_z_Sp(z) /ILF_of_z_Sp(0.5*(zrange[kbin].first+zrange[kbin].second)) );
    rsfSB.scalePhistar(ILF_of_z_SB(z) /ILF_of_z_SB(0.5*(zrange[kbin].first+zrange[kbin].second)) );
    double rE = rsfE(mag);
    double rSp = rsfSp(mag);
    double rSB = rsfSB(mag);
    double rSum = rE+rSp+rSB;
    rE /= rSum;  rSp /= rSum;  rSB /= rSum;
    double rnd=rgp_->Flat01();
    if (rnd < rE)  typ123 = 1;  // Elliptical
    else if (rnd < rSp+rE)  typ123 = 2;  // Spiral
    else typ123 = 3;  // StarBurst  

    // simulate more precise type
    switch (typ123) {
    case 1 : {
      double typ1=floor(a1+(b1-a1)*rgp_->Flat01()); 
      type=1.0+typ1/100.0;}
      break;
    case 2 :
      { double typ2=floor(a2+(b2-a2)*rgp_->Flat01()); 
	type=2.0+typ2/100.0; }
      break;
    case 3 :
      { double typ3=floor(a3+(b3-a3)*rgp_->Flat01()); 
	type=3.0+typ3/100.0; }
      break;
    default:
      throw RangeCheckError("Mass2Gal::DrawMagType()/ERROR bad type < 1 or > 3 !");
      break;
    }
 
    /*
    switch (typ123) {
    case 1 : 
      { type=1.0;}
      break;
    case 2 :
      {  
	// Dahlen
	// double typ2=rgp_->Flat01();
	//	if (typ2 < 0.4)  type=2.1 ;
	//if (typ2 < 0.8 && typ2 >= 0.4) type=2.2 ;
	//   if (typ2>= 0.8) type=2.3;

	// Zucca
	type=2.1;
      }
	break;
    case 3 :
      { double typ3=rgp_->Flat01(); 
	// Dahlen
	type= (typ3 < 0.5 ? 3.4 : 3.5);

	// Zucca
	if (typ3 < 0.25)  type=3.2 ;
	if (typ3 < 0.50 && typ3 >= 0.25) type=3.3 ;
	if (typ3 < 0.75 && typ3 >= 0.50) type=3.4 ;
	if (typ3 >= 0.75)                type=3.5 ;
      }
      break;
    default:
      throw RangeCheckError("Mass2Gal::DrawMagType()/ERROR bad type < 1 or > 3 !");
      break;
    }
     */

  //Ell : 1.,0.,0.,0.,0.,0.   (100% en SED El_cww ) 
  //Spiral: 0.,0.4,0.4,0.2,0.,0.  (40% Scd_cww , 40%  Sbc_cww , 20% Im_cww)
  //StarBurst: 0,0,0,0,0.5,0.5  (50% SB2_kin , 50% SB3_kin  ),
  }
  return ;
}


ostream& MultiType_Z_LF::print(ostream& os)  const
{
  os << "# LF_Params: phistar  MStar  alpha ]   magMin="<<magMin_<<" magMax="<<magMax_<<endl;
  os << "# z-range  LF_Params_Ell  LF_Params_Spiral  LF_Params_StarBurst  LF_Params_All "<<endl;
  for(size_t k=0; k<zrange.size(); k++) {
    os << setw(8) << zrange[k].first << "  " << zrange[k].second << " \t " << vshEll[k] << " \t " << vshSp[k]
       << " \t " << vshSB[k] << " \t " << vshAll[k] << endl;
  }
  return os; 
}



