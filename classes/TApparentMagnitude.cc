/*
 * TApparentMagnitude.C
 *
 *  Created on: Oct 28, 2010
 *      Author: gorecki
 */

#include "TApparentMagnitude.h"
#include "TDataCard.h"
#include "TBFilter.h"

TApparentMagnitude::TApparentMagnitude(bool p) :
  MAGNULL(99), ERRMAGNULL(2), celerite(3.e8) {
  noebv = 0;
  gamma = new double[6];
  int indx = 0;
  gamma[indx++] = 0.037;
  gamma[indx++] = 0.038;
  gamma[indx++] = 0.039;
  gamma[indx++] = 0.039;
  gamma[indx++] = 0.040;
  gamma[indx++] = 0.040;
  err_syst = 0.005;
  /*
    cout << endl
    << "###########################################" << endl
    << " Syst err = " << err_syst << endl
    << "###########################################"
    << endl
    << endl;
  */
  zmin=0;
  zmax=10;
  ebvmin=0;
  ebvmax=10;
  typemin=0;
  typemax=1000;

  
  distance = new TDistance();
  distance -> SetCosmologicalParameter();
  print = p;
  
  rand = new TRandom3();
  rand->SetSeed(0);
}


void TApparentMagnitude::Write()
{}



TApparentMagnitude::~TApparentMagnitude() {
  delete gamma;
  delete F0;
  delete fnvis;
  delete fm5;
}

void TApparentMagnitude::LoadDataCard(TDataCard *card){
  if (card->ebvtreatment=="noebv")
    noebv = 1;
  SetSurvey(card);
  
  SetF0(card);
  
  if (survey=="CFHTLS")
    LoadCFHTLSErrorMag(card->CFHTLS_ErrMagFile);
  if (survey=="LSST")
    SetNvisit(card->nVisit, card->nfilter);
  SetFlim();
  
  //SetF0(card);
  
  sed = new TSed(print);
  sed -> SetTypes(card);
  sed -> LoadType(card);
}

TSed *TApparentMagnitude::GetSed(){
  return sed;
}

double TApparentMagnitude::GetEBV(int type){
  if (noebv)
    return 0;
  else
    return sed -> GetEBV(type);
}

void TApparentMagnitude::SetFlim(){
  Flim = new double[nband];
  maglim = new double[nband];
  
  for (int i=0; i<nband; i++){
    // flux limits = F5sigma/5

    Flim[i] = 0;
    maglim[i] = 50;
    
    if (survey=="LSST"){
      Flim[i] = 1./5*pow(10,-0.4*fm5[i])*F0[i]; 
      //cout << "        " << i << " " << Flim[i] << " " << F0[i] << " " << fm5[i] << endl;
      maglim[i] = -2.5 * log(Flim[i]/F0[i])/log(10);
      double x = pow(10, 0.4 * (maglim[i] - fm5[i]));
      double sig = sqrt((0.04 - gamma[i]) * x + gamma[i] * x * x);
      // cout <<  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! m5 = " << fm5[i] << "  MAG LIMITS " << maglim[i] << " +- " << sig <<  " F0 = " << F0[i] << endl;
    }
    FLUXMAGNULL[i] = Flim[i];
  }
}

void TApparentMagnitude::SetSurvey(TDataCard *card){
  survey = card->survey;
  if (survey=="lsst")
    survey="LSST";
  if (survey=="cfhtls")
    survey="CFHTLS";
  nband = card->nfilter;
}

void TApparentMagnitude::LoadCFHTLSErrorMag(string filename){
  fcfh = new TFile(filename.c_str());
  for (int i = 0; i < 5; i++) {
    hsigma_cfh[i] = (TH2F*) fcfh -> Get(Form("hsigmamag%i", i));
  }
}

void TApparentMagnitude::SetParRange(double *vec){
  int i=0;
  zmin = vec[i++];
  zmax = vec[i++];
  typemin = int(vec[i++]);
  typemax = int(vec[i++]);
  ebvmin = vec[i++];
  ebvmax = vec[i++];
}

bool TApparentMagnitude::ParametersInRange(float z, int type, float ebv){
  bool ok=1;
  bool warning = 0;
  if (z<zmin || z>zmax){
    if (warning)
      cout << "WARNING z= " << z << " not in range [" << zmin << ", " << zmax << "] ... skipping" << endl;
    ok=0;
  }
  if (type<typemin || type>typemax){
    if (warning)
      cout << "WARNING type= " << type << " not in range [" << typemin << ", " << typemax << "] ... skipping" << endl;
    ok=0;
  }
  if (ebv<ebvmin || ebv>ebvmax){
    if (warning)
      cout << "WARNING ebv= " << ebv << " not in range [" << ebvmin << ", " << ebvmax << "] ... skipping" << endl;
    ok=0;
  }
  return ok;
}


int TApparentMagnitude::GetNBands(){
  return nband;
}

void TApparentMagnitude::SetKcorrectionTable(TKcorrection *table) {
  KcorrectionTable = table;
}

void TApparentMagnitude::SetFluxTable(TFlux *table) {
  FluxTable = table;
}

void TApparentMagnitude::SetF0(TFilter *filter){
  if (F0==NULL)
    F0 = new double[nband];
  
  for (int ifilter = 0; ifilter < nband; ifilter++) {
    double integral_filter = filter->Integral_Eval_FSurlambdaDlambda(ifilter);
    F0[ifilter] = integral_filter * celerite * pow(10., -0.4 * 46.1);
  }
  // cout << "integral B : " << integral_Bfilter << endl;
}

void TApparentMagnitude::SetF0(TDataCard *card){
  TFilter *filter = new TFilter();
  TBFilter *bfilter = new TBFilter(card);
  
  filter->SetConversionUnitLambda(card->filterlambdaunit);
  filter->SetFileType(card->filterformat);
  filter->LoadFile(card->nfilter, card->filterlist);
  
  F0 = new double[nband];
  
  for (int ifilter = 0; ifilter < nband; ifilter++) {
    double integral_filter = filter->Integral_Eval_FSurlambdaDlambda(ifilter);
    F0[ifilter] = integral_filter * celerite * pow(10., -0.4 * 46.1);
    if (print)
      cout << "F0[ " << ifilter << "] = " << F0[ifilter] << endl;
  }
  integral_Bfilter = bfilter->Integral_Eval_FSurlambdaDlambda(0);
  delete filter;
  delete bfilter;
}

void TApparentMagnitude::SetNewFilter(){
  GetFilters() -> NewFilters();
  SetF0(KcorrectionTable -> GetFilters());
}

TFilter* TApparentMagnitude::GetFilters(){
  return KcorrectionTable -> GetFilters();
}

void TApparentMagnitude::SetDistance(TDistance *dist) {
  distance = dist;
}

void TApparentMagnitude::SetNvisit(double *nvis, int n) {
  /*LSST Science Book uncertainties on apparent magnitude parameters*/
  double msky[] = { 21.8, 22.0, 21.3, 20.0, 19.1, 17.5 };
  double theta[] = { 0.77, 0.73, 0.70, 0.67, 0.65, 0.63 };
  double C[] = { 23.6, 24.57, 24.57, 24.47, 24.19, 23.74 };
  double k[] = { 0.48, 0.21, 0.10, 0.07, 0.06, 0.06 };
  
  double airmass = 1.2;
  
  fm5 = new double[n];

  for (int i = 0; i < n; i++) {
    fm5[i] = C[i] + 0.5 * (msky[i] - 21.) + 2.5 * log(0.7 / theta[i])
      / log(10.) + 1.25 * log(30 * nvis[i] / 30.) / log(10.) - k[i]
      * (airmass - 1);
    
    if (print)
      cout << "Nvisit : " << i << " " << nvis[i] << " fm5 = " << fm5[i] << endl;
  }
}


double TApparentMagnitude::MD(double z) {
  double return_value;
  if (z <= 0) {
    return 0;
  }
  else {
    double distance_lum = distance-> LuminosityDistance(z)
      * distance -> ConversionToParsec() * 1.e-6;
    return_value = 5. * log(distance_lum * 1.e5) / log(10.);
    return return_value;
  }
}

int TApparentMagnitude::GetNTypes(){
  return typemax-typemin;
}


void TApparentMagnitude::ComputeApparentMagnitudes(double *mag, double *err_mag, double mag_abs,
						   int type_s,double z_s, double ebv_s) {

  ComputeApparentMagnitudeTheorique(mag, mag_abs, type_s, z_s, ebv_s);
  ComputeApparentMagnitudeErrors(mag, err_mag);
  ComputeApparentMagnitudeWithErrors_FluxGaus(mag, err_mag);
  ComputeApparentMagnitudeErrors(mag, err_mag);
}


void TApparentMagnitude::EvalApparentMagnitudes(double *mag, double *err_mag, double mag_abs,
						int type_s,double z_s, double ebv_s, int ib) {
  EvalApparentMagnitudeTheorique(mag, mag_abs, type_s, z_s, ebv_s, ib);
  ComputeMagnitudeWithErrors(mag, err_mag,ib);
}

void TApparentMagnitude::ComputeMagnitudeWithErrors(double *mag, double *err_mag, int ib){
  ComputeApparentMagnitudeErrors(mag, err_mag, ib);
  ComputeApparentMagnitudeWithErrors_FluxGaus(mag, err_mag, ib);
  ComputeApparentMagnitudeErrors(mag, err_mag, ib);
}


void TApparentMagnitude::ComputeApparentMagnitudeErrors(double *mag,
							double *error,
							int ib) {
  if (ib<0)
    for (int i = 0; i < nband; i++) {
      if (survey=="LSST"){
	if (mag[i]>maglim[i]){
	  error[i] = ERRMAGNULL;
	}
	else{
	  double x = pow(10, 0.4 * (mag[i] - fm5[i]));
	  double sig = sqrt((0.04 - gamma[i]) * x + gamma[i] * x * x);
	  error[i] = sqrt(sig * sig + err_syst * err_syst);
	}
      }
      if (survey=="CFHTLS")
	error[i] = CFHTLS_errormag(mag[i], i);
    }
  else{
    if (mag[ib]>maglim[ib])
      error[ib] = ERRMAGNULL;
    else{
      double x = pow(10, 0.4 * (mag[ib] - fm5[ib]));
      double sig = sqrt((0.04 - gamma[ib]) * x + gamma[ib] * x * x);
      error[ib] = sqrt(sig * sig + err_syst * err_syst);
    }
  }
}

void TApparentMagnitude::ComputeApparentMagnitudeWithErrors_FluxGaus(double *magnitude,
								     double *error,
								     int ib) {

  if (ib<0)
    for (int i = 0; i < nband; i++) {
      if (magnitude[i]>maglim[i]){
	magnitude[i]=MAGNULL;
      }
      else{
	double flux = F0[i] * pow(10., -0.4 * magnitude[i]);
	double err_flux = 0.4 * error[i] * flux * log(10.);
	double err = rand->Gaus(0, err_flux);
	flux += err; //flux+=r.Gaus(0, err_flux);
	if (flux < Flim[i])
	  magnitude[i] = MAGNULL;
	else {
	  magnitude[i] = -2.5 * log(flux / F0[i]) / log(10.);
	}
      }
    }
  else{
    if (magnitude[ib]>maglim[ib]){
	magnitude[ib]=MAGNULL;
      }
      else{
	double flux = F0[ib] * pow(10., -0.4 * magnitude[ib]);
	double err_flux = 0.4 * error[ib] * flux * log(10.);
	double err = rand->Gaus(0, err_flux);
	flux += err; //flux+=r.Gaus(0, err_flux);
	if (flux < Flim[ib])
	  magnitude[ib] = MAGNULL;
	else {
	  magnitude[ib] = -2.5 * log(flux / F0[ib]) / log(10.);
	}
      }
  }
}


double TApparentMagnitude::ComputeAbsMagFromMag(int ifilter,
						double mag,
						int type,
						double z,
						double ebv) {
  double k = KcorrectionTable -> EvalFromTable(type, z, ebv, ifilter);
  double m_d =  MD(z);
  double AM = mag - m_d - k;
  return AM;
}

void TApparentMagnitude::EvalFluxTheorique(double *flux,
					   double MA,
					   int type,
					   double z,
					   double ebv) {
  double *mag = new double[nband];
  EvalApparentMagnitudeTheorique(mag, MA, type, z, ebv);
  ComputeFlux(mag, flux);
  delete mag;
}

void TApparentMagnitude::EvalApparentMagnitudeTheorique(double *mag,
							double MA,
							int type,
							double z,
							double ebv,
							int ib) {
  double m_d =  MD(z);
  if (ib<0)
    for (int ifilter = 0; ifilter < nband; ifilter++) {
      double k = KcorrectionTable -> EvalFromTable(type, z, ebv, ifilter);
      mag[ifilter] = MA + m_d + k;
      //  cout << ifilter << " mag = " << mag[ifilter] << " = "<< MA << " + " << m_d << " + " << k << endl;
    }
  else{
    double k = KcorrectionTable -> EvalFromTable(type, z, ebv, ib);
    mag[ib] = MA + m_d + k;
  }
  // cout << " mag[3] = " << mag[3] << endl;
}

void TApparentMagnitude::ComputeApparentMagnitudeTheorique(double *mag,
							   double MA,
							   int type,
							   double z,
							   double ebv) {
  KcorrectionTable -> GetSed() -> SetEBV(ebv);
  KcorrectionTable -> GetSed() -> SetEmittedRedshift(z);
  KcorrectionTable -> SetSedNumber(type);
    
  for (int ifilter = 0; ifilter < nband; ifilter++) {
    double k = KcorrectionTable -> Eval(z, ifilter);
    //cout<<" k = "<<k<<endl;
    double m_d =  MD(z);
    mag[ifilter] = MA + m_d + k;
  }
}

void TApparentMagnitude::ObservedBand_FluxGaus(double *magnitude, double *error,
					       int *band_observed) {
  for (int i = 0; i < nband; i++) {
    //cout << i << " " << magnitude[i] << "<?" << maglim[i] << " " << error[i] << endl;
    if (magnitude[i] >= maglim[i] || error[i] > 0.2)
      band_observed[i] = 0;
    else
      band_observed[i] = 1;
  }
}

double TApparentMagnitude::CFHTLS_errormag(double mag, int imag){
  double sig = ERRMAGNULL;
  int ibin = hsigma_cfh[imag] -> GetXaxis() -> FindBin(mag);
  int nbin = hsigma_cfh[imag] -> GetXaxis() -> GetNbins();

  // cout << "bin = " << ibin << " / " << nbin << endl;

  if (ibin <= nbin) {
    double integral = hsigma_cfh[imag]->Integral(ibin, ibin, 0, nbin);
    
    // cout << "integral : " << integral << endl;
    
    if (integral > 0) {
      TH1D *Pdelta = hsigma_cfh[imag]->ProjectionY(Form("hsigmamag%i_py", imag), ibin, ibin, "");
      if (Pdelta -> GetEntries() != 0){
        sig = TMath::Exp(Pdelta -> GetRandom());
      }
      Pdelta -> Delete("");
    }
  }
  
  if (sig>1)
    sig = ERRMAGNULL;
  // cout << "error = " << sig << endl;

  return sig;
}


double TApparentMagnitude::AbsoluteMagnitude(double redshift, double ratio,
					     int type, double ebv, int ifilter) {
  double return_value;
  double integral_sed = FluxTable -> EvalFromTable(type, redshift, ebv, ifilter);
  
  return_value = -2.5 * log((1 + redshift) / celerite * integral_sed
			    / FluxTable -> EvalFromTable(type, 0, ebv, ifilter) * ratio) / log(10.) - MD(redshift) - 46.1;
  return return_value;
}

void TApparentMagnitude::ComputeFlux(float *mag, float *sigma_mag,
				     double *flux, double *sigma_flux, int i) {

  if (mag[i]>=maglim[i]){
    flux[i] = FLUXMAGNULL[i]; 
    sigma_flux[i] = FLUXMAGNULL[i];
    return;
  }
  
  flux[i] = F0[i] * pow(10., -0.4 * mag[i]);
  sigma_flux[i] = 0.4 * sigma_mag[i] * flux[i] * log(10.);
  
  if (flux[i] < Flim[i]){
    flux[i] = FLUXMAGNULL[i];
    sigma_flux[i] = FLUXMAGNULL[i];
  }
  return; 
}
void TApparentMagnitude::ComputeFlux(float *mag, float *sigma_mag,
				     double *flux, double *sigma_flux) {
  for (int i=0; i<nband; i++){
    if (mag[i]>=maglim[i]){
      flux[i] = FLUXMAGNULL[i]; 
      sigma_flux[i] = FLUXMAGNULL[i];
      continue;
    }
    
    flux[i] = F0[i] * pow(10., -0.4 * mag[i]);
    sigma_flux[i] = 0.4 * sigma_mag[i] * flux[i] * log(10.);
    
    if (flux[i] < Flim[i]){
      flux[i] = FLUXMAGNULL[i];
      sigma_flux[i] = FLUXMAGNULL[i];
    }
  }
}

void TApparentMagnitude::ComputeFlux(float *mag, double *flux) {
  for (int i=0; i<nband; i++){
    if (mag[i]<maglim[i]){
      flux[i] = F0[i] * pow(10., -0.4 * mag[i]);
    }
    else{
      flux[i]=Flim[i];
    }
  }
  return; 
}
void TApparentMagnitude::ComputeFlux(double *mag, double *flux) {
  for (int i=0; i<nband; i++){
    if (mag[i]<maglim[i]){
      flux[i] = F0[i] * pow(10., -0.4 * mag[i]);
    }
    else{
      flux[i]=Flim[i];
    }
  }
  return; 
}

float TApparentMagnitude::ComputeMag(double flux, int i) {
  if (flux<=0)
    return 99;
  else
    return -2.5*log10(flux/F0[i]);
}


double TApparentMagnitude::GetFLUXMAGNULL(int i){
  return FLUXMAGNULL[i];
}


double TApparentMagnitude::ComputeFluxTemplate(double redshift, int type,
					       double ebv, int i) {
  double flux = FluxTable -> EvalFromTable(type, redshift, ebv, i);
  return flux;
}
