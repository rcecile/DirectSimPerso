#include "powerspecfromgrids.h"

/*
  * herited from Alex Abate powerspec
  * but simplified with better shot-noise computation
  */
PowerSpecFromGrids::PowerSpecFromGrids(TArray<r_8> galdens, double Dx, double ratio_AngDiam)
: drho_(galdens)
// Reads in galaxy fluctuation grid and computes Fourier transform
// The Dkx, Dky, Dkz values are computed from the array size and Dx
{
	cout << "    PowerSpecFromGrids::PowerSpecFromGrids():"<<endl;
	
	Dx_=Dx, Dy_=Dx, Dz_=Dx;
	ratio_AngDiam_ =ratio_AngDiam;
   	ComputeFourier(drho_,four_);
	SetDKDR();

	zc_=-1;
	
	cout << "    Printing Fourier grid info ...."<<endl;
	cout << "    Nx,Ny,Nz = "<< Nx_ <<","<< Ny_ <<","<< Nz_ <<endl;
	cout << "    Dkx,Dky,Dkz = "<< Dkx_ <<","<< Dky_ <<","<< Dkz_ <<endl;
	cout << "    cell size = dx,dy,dx = "<< Dx_ <<","<< Dy_ <<","<< Dz_ <<endl;

	cout << "    Printing drho_ array info (should match Nx,Ny,Nz above) ...."<<endl;
	cout << "    SizeX,SizeY,SizeZ = "<< drho_.SizeX() <<","<< drho_.SizeY();
	cout <<","<< drho_.SizeZ() <<endl;

	cout << "    Printing four_ array info ...."<<endl;
	cout << "    SizeX,SizeY,SizeZ = "<< four_.SizeX() <<","<< four_.SizeY();
	cout << ","<< four_.SizeZ() <<endl;
	cout << "    EXIT PowerSpecFromGrids::PowerSpecFromGrids():"<<endl<<endl;
	
};

void PowerSpecFromGrids::ComputeFourier(TArray<r_8>& real_array, TArray< complex< r_8 > >& four_array)
{
    // FFT galaxy flunction grid
	FFTWServer ffts;
	ffts.FFTForward(real_array, four_array);

	cout << "    PowerSpecFromGrids::ComputeFourier(): "<<endl;
	cout << "    check size of arrays ..."<<endl;
	cout << "    size of drho = [sizex sizey sizez] = ["<< real_array.SizeX();
	cout << "  "<< real_array.SizeY() <<"  "<< real_array.SizeZ() <<"]"<<endl; 
	cout << "    size of four = [sizex sizey sizez] = ["<< four_array.SizeX();
	cout << "  "<< four_array.SizeY() <<"  "<< four_array.SizeZ() <<"]"<<endl;
	cout << "    EXIT PowerSpecFromGrids::ComputeFourier(): "<<endl<<endl;
	return;

};


void PowerSpecFromGrids::ComputeFourierBack(TArray< complex< r_8 > >& four_array, TArray<r_8>& real_array)
{
	// FFT convolved galaxy grid
	FFTWServer ffts;
	ffts.FFTBackward(four_array, real_array);

	cout << "    PowerSpecFromGrids::ComputeFourierBack(): "<<endl;
	cout << "    check size of arrays ..."<<endl;
	cout << "    size of drho = [sizex sizey sizez] = ["<< real_array.SizeX();
	cout << "  "<< real_array.SizeY() <<"  "<< real_array.SizeZ() <<"]"<<endl; 
	cout << "    size of four = [sizex sizey sizez] = ["<< four_array.SizeX();
	cout << "  "<< four_array.SizeY() <<"  "<< four_array.SizeZ() <<"]"<<endl;
	cout << "    EXIT PowerSpecFromGrids::ComputeFourierBack(): "<<endl<<endl;
	return;
};


void PowerSpecFromGrids::SetDKDR()
//Read in info on FT grid
{
	Nx_ = drho_.SizeX();// will be short dimension
	Ny_ = drho_.SizeY();
	Nz_ = drho_.SizeZ();
	Dkx_ = (2*PI)/(Nx_*Dx_*ratio_AngDiam_);
	Dky_ = (2*PI)/(Ny_*Dy_*ratio_AngDiam_);
	Dkz_ = (2*PI)/(Nz_*Dz_*ratio_AngDiam_);

	NRtot_ = Nx_*Ny_*Nz_; // number of pixels in the survey
	NCx_ = Nx_/2+1; // think this should be the length of the short dimension in the four_ array
	// Kny = 2/Dx // Nyquist frequency of grid, should be sufficiently greater than klin
	// klin~0.1h/Mpc at z=0, ~0.2h/Mpc at z=1
	
};


double PowerSpecFromGrids::AccumulatePowerSpectra(HProf& hp)
// Compute spectrum from "four_" and fill profile histogram "hp"
// Power spectrum in hp is NOT yet normalised
// Radial direction is assumed to be the third dimenion (z-coord)
// four_ : in the format: four_(nkx,nky,nkz)
// FIRST dimension is short dimension
 { 
   cout << "    PowerSpecFromGrids::AccumulatePowerSpectra():"<<endl;
   
   if(hp.NBins()<0) 
     throw ParmError("ERROR! HProf bins undefined");
   hp.Zero();
  
   double sum=0;
   // wavenumber k = n * 2pi/L where n=wave index and L = length of grid side
   
   for(sa_size_t iz=0; iz<four_.SizeZ(); iz++) { // assumed RADIAL direction
     
     // get wave number 3RD/RADIAL dim
     double kz = iz;// wave index for +ve freq
     if (iz > four_.SizeZ()/2) 
       kz = four_.SizeZ()-iz;// wave index for -ve freq
     kz *= Dkz_; // wave number
     
     
     for(sa_size_t iy=0; iy<four_.SizeY(); iy++) { 
       
       // get wave number 2ND dim
       double ky = iy; // wave index for +ve freq
       if (iy > four_.SizeY()/2) 
	 ky = four_.SizeY()-iy; // wave index for -ve freq
       ky *= Dky_; // wave number
       
       for(sa_size_t ix=0; ix<four_.SizeX(); ix++) { // THIS IS THE SHORT DIMENSION, no neg freq
	 
	 // get wave number 1ST dim (straightforward: no neg freq)
	 double kx = ix*Dkx_;
	 	 
	 // k modulus (wavevector length)
	 double kmod = sqrt((kx*kx+ky*ky+kz*kz));
	 
	 // Fourier component
	 complex< r_8 > za = four_(ix, iy, iz);
	 
	 // Fourier component * its complex conjugate
	 double pk    = za.real()*za.real()+za.imag()*za.imag();
	  
	 hp.Add(kmod, pk);  
	 sum+=pk; 
	 
       }// end loop over X dim
     }// end loop over Y dim
   }// end loop over Z dim

   cout <<"    Sum of Fourier coefficients sq ="<< sum <<endl;
  
   return sum;
   cout << "    EXIT PowerSpecFromGrids::AccumulatePowerSpectra():"<<endl<<endl;
   
 };


void PowerSpecFromGrids::WritePS(string fname, HProf& Pdata, r_4 Voldata, HProf& Pdata_noise, double h, double nGalGrid, int nGridData,bool comp_SN)
// Writes power spectra to a text file
// the format is:
// [k values] 
// [normalised galaxy PS] 
// [normalized shot-noise PS] 
// [error on normalised galaxy PS] 
{
	cout << "    PowerSpec::WritePS()" <<endl;
	Histo Pdatah=Pdata.GetHisto();	
	Histo Pnoiseh=Pdata_noise.GetHisto();
	
	Histo fracmodok;

	cout <<"check Nbins : " << Pdatah.NBins() <<endl;
	
	// number of k bins	
	int_4 nbk=Pdatah.NBins();
	
	cout << "    Number of bins="<<nbk<<endl;
	
	//read histogram values into a file
	ifstream inp;
	ofstream outp;
	inp.open(fname.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
	  inp.clear(ios::failbit);
	  cout << "Writing to file ..." << fname.c_str() << endl<< endl<< endl;
	  outp.open(fname.c_str(), ofstream::out);
	  
	  outp << "@HUBBLE " << h << " Volume of data per grid = "<< Voldata ;
	  outp << ", total number of galaxies in weight grids = " << nGalGrid <<endl;
	  outp << "@Nx,Ny,Nz,R(Mpc),zc = "<< Nx_ <<","<< Ny_ <<","<< Nz_ <<",";
	  outp << Dx_ <<","<< zc_ << ", number of data grids = " << nGridData << endl;
	  
	  double sigma_factor = 2. * PI * sqrt(1. / Nx_ / Ny_ / Nz_ / Dx_ / Dy_ / Dz_ / nGridData / Pdatah.BinWidth());
	  r_8 kvals,Pnorm,Pnoise,sigmaP,SNmean=0.;
	  
	  if (comp_SN) {  // SN is the constant part of the PS of the Poisson grid
	    int nSNmean = 0;
	    for(int_4 i=0;i<nbk;i++) 
	      if (Pdatah.BinCenter(i) > 0.2 &&  Pdatah.BinCenter(i) <  0.4) {
		SNmean += Pnoiseh.operator()(i) * Voldata * nGridData;
		nSNmean ++;
	      }
	    SNmean /= (double) nSNmean;
	    cout << "Shot-noise level computed from the Poisson grid = "<< SNmean << endl;
	  }

	  for(int_4 i=0;i<nbk;i++) {
	    
	    kvals = Pdatah.BinCenter(i);			
	    Pnorm  = Pdatah.operator()(i)  * Voldata * nGridData;
	    sigmaP = Pnorm / kvals * sigma_factor;

	    if (comp_SN) Pnoise = SNmean;
	    else Pnoise = Pnoiseh.operator()(i) * Voldata * nGridData;		
	    
	    outp << kvals <<"   \t"<< Pnorm <<"   \t" << Pnoise <<"   \t" << sigmaP << endl;
	  }
	  outp.close();
	}
	else
	  cout << "Error...file """ << fname.c_str() << """ exists" << endl;
	
	cout << "    EXIT PowerSpecFromGrids::WritePS()"<<endl;
};
