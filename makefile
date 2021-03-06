################################################################################
##  makefile                                                                  ##
################################################################################

ifndef SOPHYABASE
	SOPHYABASE := /sps/lsst/users/rcecile/BAOProgs/LUC/SrcSophyaPersoPourTests/SObjs/
	echo "sophyabase = ${SOPHYABASE}"
endif

include $(SOPHYABASE)/include/sophyamake.inc

#### Automatically set up GALSIM variable if it doesn't exist
ifndef GALSIM
	GALSIM := ${PWD}
#	echo "GALSIM = ${GALSIM}"
endif

#define GEN3D_FLOAT

OBJ = ${GALSIM}/objs

EXE = ${GALSIM}/exe

BAOPROGS = ${GALSIM}/baoprogs

ROOTOUT = ${GALSIM}/root

TESTS =${GALSIM}/testfiles

MYCL = ${GALSIM}/classes

PROGS = ${GALSIM}/progs

KCORR = ${GALSIM}/kCorrections

LIBH := $(MYCL)/cosmocalcs.h $(MYCL)/geneutils.h $(MYCL)/gftdist.h \
$(MYCL)/schechter.h $(MYCL)/sinterp.h $(MYCL)/simdata.h $(MYCL)/reddening.h \
$(MYCL)/sedfilter.h $(MYCL)/sedpca.h $(MYCL)/genefluct3d.h  $(MYCL)/pkspectrum.h \
$(MYCL)/mass2gal.h $(MYCL)/powerspec.h $(MYCL)/matrix.h $(MYCL)/igm.h \
$(MYCL)/myschechter.h  $(MYCL)/multitypzlf.h   $(MYCL)/powerspecfromgrids.h \
$(MYCL)/hpoly.h $(MYCL)/shapelets.h $(MYCL)/em.h $(MYCL)/cat2grid.h  $(MYCL)/projgrid.h \
$(MYCL)/fitkbaoscale.h $(MYCL)/fitkbaobaselinescale.h $(MYCL)/chisqstats.h \
$(MYCL)/TAM.h $(MYCL)/TApparentMagnitude.h $(MYCL)/TDataCard.h $(MYCL)/TFilter.h $(MYCL)/TIGM.h \
$(MYCL)/TReddening.h $(MYCL)/TBFilter.h $(MYCL)/TDistance.h $(MYCL)/TFlux.h \
$(MYCL)/TKcorrection.h  $(MYCL)/TSed.h
#$(MYCL)/constcosmo.h
#$(MYCL)/root_plots.h

LIBO := $(OBJ)/cosmocalcs.o $(OBJ)/geneutils.o $(OBJ)/gftdist.o \
$(OBJ)/schechter.o $(OBJ)/sinterp.o $(OBJ)/simdata.o $(OBJ)/reddening.o \
$(OBJ)/sedfilter.o $(OBJ)/sedpca.o $(OBJ)/genefluct3d.o  $(OBJ)/pkspectrum.o $(OBJ)/mass2gal.o \
$(OBJ)/myschechter.o  $(OBJ)/multitypzlf.o  $(OBJ)/powerspecfromgrids.o \
$(OBJ)/powerspec.o $(OBJ)/matrix.o $(OBJ)/igm.o $(OBJ)/hpoly.o $(OBJ)/shapelets.o\
$(OBJ)/em.o $(OBJ)/cat2grid.o $(OBJ)/fitkbaoscale.o $(OBJ)/fitkbaobaselinescale.o $(OBJ)/chisqstats.o \
$(OBJ)/TAM.o $(OBJ)/TApparentMagnitude.o $(OBJ)/TDataCard.o $(OBJ)/TFilter.o $(OBJ)/TIGM.o \
$(OBJ)/TReddening.o $(OBJ)/TBFilter.o $(OBJ)/TDistance.o $(OBJ)/TFlux.o \
$(OBJ)/TKcorrection.o  $(OBJ)/TSed.o
#$(OBJ)/root_plots.o

# root libraries
ROOTLIB = $(shell root-config --libs) -lMathMore
ROOTINC = $(shell root-config --incdir)
MYLIB = -lfftw3_threads  
MINUIT = -lMinuit


################################################################################

#all : progs tests bao
all : progs bao

progs : addIGMToSED analyzeBPZ baseSimulation calculateKcorrections cfhtColors \
colorDistributions convertSEDS fitLSSTspectra lineOfSightLymanAlpha lineOfSightMagnitude \
lsstPicklesLibrary lymanAlphaToDensity pcaTemplates photoZdist priorFitter \
make_cat sdssElColors sdssPicklesLibrary simdensity  \
simulateAbsorberLinesOfSight simulateLSSTobs simulateLSSTobsFromTruth 

tests : test2Dinterp testbasesim testEMalgorithm testErrors testgoodsmagsim \
testKcorrColors testKcorrMethod testLF testLymanAlphaAbs testMadau testMeiksin \
testSimReadKcorr testsimulateIGM testSimulation testTemplateFitting 

bao : addGausszerr computepsfromgrids fitkbao fitkbaobaseline getpzconvf getsf grid_data \
 subfromfull 

clean : 
	rm  $(OBJ)/* $(EXE)/*

# MAIN PROGS

addIGMToSED : $(EXE)/addIGMToSED
	@echo 'makefile : addIGMToSED made'

analyzeBPZ : $(EXE)/analyzeBPZ
	@echo 'makefile : analyzeBPZ made'

baseSimulation : $(EXE)/baseSimulation
	@echo 'makefile : baseSimulation made'

calculateKcorrections : $(EXE)/calculateKcorrections
	@echo 'makefile : calculateKcorrections made'

cfhtColors : $(EXE)/cfhtColors
	@echo 'makefile : cfhtColors made'

colorDistributions	: $(EXE)/colorDistributions
	@echo 'makefile : colorDistributions made'

convertSEDS	: $(EXE)/convertSEDS
	@echo 'makefile : convertSEDS made'

fitLSSTspectra : $(EXE)/fitLSSTspectra
	@echo 'makefile : fitLSSTspectra made'

lineOfSightLymanAlpha : $(EXE)/lineOfSightLymanAlpha
	@echo 'makefile : lineOfSightLymanAlpha made'

lineOfSightMagnitude : $(EXE)/lineOfSightMagnitude
	@echo 'makefile : lineOfSightMagnitude made'

lsstPicklesLibrary : $(EXE)/lsstPicklesLibrary
	@echo 'makefile : lsstPicklesLibrary made'

lymanAlphaToDensity : $(EXE)/lymanAlphaToDensity
	@echo 'makefile : lymanAlphaToDensity made'

pcaTemplates : $(EXE)/pcaTemplates
	@echo 'makefile : pcaTemplates made'

photoZdist : $(EXE)/photoZdist
	@echo 'makefile : photoZdist made'

priorFitter : $(EXE)/priorFitter
	@echo 'makefile : priorFitter made'

#projectTemplates : $(EXE)/projectTemplates
#	@echo 'makefile : projectTemplates made'

make_cat : $(EXE)/make_cat
	@echo 'makefile : make_cat made'

sdssElColors : $(EXE)/sdssElColors
	@echo 'makefile : sdssElColors made'

sdssPicklesLibrary : $(EXE)/sdssPicklesLibrary
	@echo 'makefile : sdssPicklesLibrary made'

simdensity : $(EXE)/simdensity
	@echo 'makefile : simdensity made'

simulateAbsorberLinesOfSight : $(EXE)/simulateAbsorberLinesOfSight
	@echo 'makefile : simulateAbsorberLinesOfSight made'

simulateCFHTobs : $(EXE)/simulateCFHTobs
	@echo 'makefile : simulateCFHTobs made'

simulateLSSTobs : $(EXE)/simulateLSSTobs
	@echo 'makefile : simulateLSSTobs made'	

simulateLSSTobsFromTruth : $(EXE)/simulateLSSTobsFromTruth
	@echo 'makefile : simulateLSSTobsFromTruth made'

# BAO PROGS

addGausszerr : $(EXE)/addGausszerr
	@echo 'makefile : addGausszerr made'

computepsfromgrids : $(EXE)/computepsfromgrids
	@echo 'makefile : computepsfromgrids made'

fitkbao : $(EXE)/fitkbao
	@echo 'makefile : fitkbao made'

fitkbaobaseline : $(EXE)/fitkbaobaseline
	@echo 'makefile : fitkbaobaseline made'

getpzconvf : $(EXE)/getpzconvf
	@echo 'makefile : getpzconvf made'

getsf : $(EXE)/getsf
	@echo 'makefile : getsf made'

grid_data : $(EXE)/grid_data
	@echo 'makefile : grid_data made'

subfromfull : $(EXE)/subfromfull
	@echo 'makefile : subfromfull made'

# TESTING PROGS

test2Dinterp : $(EXE)/test2Dinterp 
	@echo 'makefile :test2Dinterp made'

testbasesim : $(EXE)/testbasesim 
	@echo 'makefile :testbasesim made'

test : $(EXE)/test
	@echo 'makefile : test made'

testEMalgorithm : $(EXE)/testEMalgorithm
	@echo 'makefile : testEMalgorithm made'

testErrors : $(EXE)/testErrors
	@echo 'makefile : testErrors made'

testgoodsmagsim : $(EXE)/testgoodsmagsim
	@echo 'makefile : testgoodsmagsim made'

testKcorrColors : $(EXE)/testKcorrColors 
	@echo 'makefile : testKcorrColors made'

testKcorrMethod : $(EXE)/testKcorrMethod 
	@echo 'makefile : testKcorrMethod made'

testLF : $(EXE)/testLF 
	@echo 'makefile :testLF made'

testLymanAlphaAbs : $(EXE)/testLymanAlphaAbs 
	@echo 'makefile :testLymanAlphaAbs made'

testMadau : $(EXE)/testMadau
	@echo 'makefile : testMadau made'

testMeiksin : $(EXE)/testMeiksin
	@echo 'makefile : testMeiksin made'

testSimReadKcorr : $(EXE)/testSimReadKcorr
	@echo 'makefile : testSimReadKcorr made'

testsimulateIGM : $(EXE)/testsimulateIGM
	@echo 'makefile : testsimulateIGM made'

testSimulation : $(EXE)/testSimulation
	@echo 'makefile : testSimulation made'

testTemplateFitting : $(EXE)/testTemplateFitting
	@echo 'makefile : testTemplateFitting made'

## programs below here have not been CHECKED or maybe even finished...

testpsdenscube : $(EXE)/testpsdenscube
	@echo 'makefile :  testpsdenscube made'

testsimdensity : $(EXE)/testsimdensity
	@echo 'makefile :  testsimdensity made'

###################### MAIN PROGRAMS ###########################################

# ADD LINE OF SIGHT TRANSMISSON TO SEDS IN A LIBRARY
$(EXE)/addIGMToSED :	$(OBJ)/addIGMToSED.o $(LIBO) 
	mkdir -p $(EXE)
	@echo 'SOPHYA LIB'
	@echo 	$(SOPHYAEXTSLBLIST)
	@echo 'SOPHYA LIB'
	$(CXXLINK) -o $(EXE)/addIGMToSED $(OBJ)/addIGMToSED.o $(LIBO) \
	$(SOPHYAEXTSLBLIST)	$(MYLIB) $(ROOTLIB)

$(OBJ)/addIGMToSED.o : $(PROGS)/addIGMToSED.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/addIGMToSED.o $(PROGS)/addIGMToSED.cc

# ANALYZE A BPZ CATALOG
$(EXE)/analyzeBPZ : $(OBJ)/analyzeBPZ.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/analyzeBPZ $(OBJ)/analyzeBPZ.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/analyzeBPZ.o : $(PROGS)/analyzeBPZ.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/analyzeBPZ.o $(PROGS)/analyzeBPZ.cc 

# BASE SIMULATION
$(EXE)/baseSimulation : $(OBJ)/baseSimulation.o $(LIBO) 
	mkdir -p $(EXE)
	$(CXXLINK) -o $(EXE)/baseSimulation $(OBJ)/baseSimulation.o $(LIBO) \
	$(SOPHYAEXTSLBLIST)	$(MYLIB) $(ROOTLIB)

$(OBJ)/baseSimulation.o : $(PROGS)/baseSimulation.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/baseSimulation.o $(PROGS)/baseSimulation.cc

# CALCULATE K CORRECTION TABLES
$(EXE)/calculateKcorrections : $(OBJ)/calculateKcorrections.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(KCORR)
	$(CXXLINK) -o $(EXE)/calculateKcorrections $(OBJ)/calculateKcorrections.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/calculateKcorrections.o : $(PROGS)/calculateKcorrections.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/calculateKcorrections.o $(PROGS)/calculateKcorrections.cc

# CALCULATE U-G, I-Z CFHT colors w/ and wo/ host galaxy reddening
$(EXE)/cfhtColors : $(OBJ)/cfhtColors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(KCORR)
	$(CXXLINK) -o $(EXE)/cfhtColors $(OBJ)/cfhtColors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/cfhtColors.o : $(PROGS)/cfhtColors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/cfhtColors.o $(PROGS)/cfhtColors.cc

# COLOR DISTRIBUTIONS
$(EXE)/colorDistributions : $(OBJ)/colorDistributions.o $(LIBO) 
	mkdir -p $(EXE)
	$(CXXLINK) -o $(EXE)/colorDistributions $(OBJ)/colorDistributions.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/colorDistributions.o : $(PROGS)/colorDistributions.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/colorDistributions.o $(PROGS)/colorDistributions.cc

# CONVERT SED UNITS
$(EXE)/convertSEDS : $(OBJ)/convertSEDS.o
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/convertSEDS $(OBJ)/convertSEDS.o $(SOPHYAEXTSLBLIST)

$(OBJ)/convertSEDS.o : $(PROGS)/convertSEDS.cc
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -o $(OBJ)/convertSEDS.o $(PROGS)/convertSEDS.cc 

# FIT LSST SPECTRA TO CWWK
$(EXE)/fitLSSTspectra : $(OBJ)/fitLSSTspectra.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/fitLSSTspectra $(OBJ)/fitLSSTspectra.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/fitLSSTspectra.o : $(PROGS)/fitLSSTspectra.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/fitLSSTspectra.o $(PROGS)/fitLSSTspectra.cc

# SIMULATE LINE OF SIGHT LYMAN ALPHA TRANSMISSION
$(EXE)/lineOfSightLymanAlpha : $(OBJ)/lineOfSightLymanAlpha.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/lineOfSightLymanAlpha $(OBJ)/lineOfSightLymanAlpha.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lineOfSightLymanAlpha.o : $(PROGS)/lineOfSightLymanAlpha.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/lineOfSightLymanAlpha.o $(PROGS)/lineOfSightLymanAlpha.cc

# U,G MAGNITUDES WITH IGM LINE OF SIGHT
$(EXE)/lineOfSightMagnitude : $(OBJ)/lineOfSightMagnitude.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/lineOfSightMagnitude $(OBJ)/lineOfSightMagnitude.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lineOfSightMagnitude.o : $(PROGS)/lineOfSightMagnitude.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/lineOfSightMagnitude.o $(PROGS)/lineOfSightMagnitude.cc

# CALCULATE LSST COLORS FOR PICKLES' LIBRARY OF STARS 
$(EXE)/lsstPicklesLibrary : $(OBJ)/lsstPicklesLibrary.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/lsstPicklesLibrary $(OBJ)/lsstPicklesLibrary.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lsstPicklesLibrary.o : $(PROGS)/lsstPicklesLibrary.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/lsstPicklesLibrary.o $(PROGS)/lsstPicklesLibrary.cc

# LYMAN ALPHA ALONG LINE OF SIGHT CONVERTED TO DENSITY
$(EXE)/lymanAlphaToDensity : $(OBJ)/lymanAlphaToDensity.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/lymanAlphaToDensity $(OBJ)/lymanAlphaToDensity.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lymanAlphaToDensity.o : $(PROGS)/lymanAlphaToDensity.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/lymanAlphaToDensity.o $(PROGS)/lymanAlphaToDensity.cc

# TEMPLATE PCA
$(EXE)/pcaTemplates : $(OBJ)/pcaTemplates.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/pcaTemplates $(OBJ)/pcaTemplates.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/pcaTemplates.o : $(PROGS)/pcaTemplates.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/pcaTemplates.o $(PROGS)/pcaTemplates.cc 

# PHOTO-Z DISTRIBUTION
$(EXE)/photoZdist : $(OBJ)/photoZdist.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/photoZdist $(OBJ)/photoZdist.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/photoZdist.o : $(PROGS)/photoZdist.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/photoZdist.o $(PROGS)/photoZdist.cc 

# PRIOR FITTER
$(EXE)/priorFitter : $(OBJ)/priorFitter.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/priorFitter $(OBJ)/priorFitter.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) $(MINUIT)

$(OBJ)/priorFitter.o : $(PROGS)/priorFitter.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/priorFitter.o $(PROGS)/priorFitter.cc 

# PROJECT TEMPLATES
#$(EXE)/projectTemplates : $(OBJ)/projectTemplates.o $(LIBO) 
#	mkdir -p $(EXE)
#	mkdir -p $(ROOTOUT)
#	$(CXXLINK) -o $(EXE)/projectTemplates $(OBJ)/projectTemplates.o $(LIBO) \
#	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

#$(OBJ)/projectTemplates.o : $(PROGS)/projectTemplates.cc $(LIBH)  
#	mkdir -p $(OBJ)
#	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/projectTemplates.o $(PROGS)/projectTemplates.cc 


# SIMULATE CATALOG OF BASIC GALAXY PROPERTIES FROM OVER-DENSITY GRID - NEW VERSION
$(EXE)/make_cat : $(OBJ)/make_cat.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/make_cat $(OBJ)/make_cat.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/make_cat.o : $(PROGS)/make_cat.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/make_cat.o $(PROGS)/make_cat.cc

# CALCULATE SDSS COLORS OF ELLIPTICAL GALAXY
$(EXE)/sdssElColors : $(OBJ)/sdssElColors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/sdssElColors $(OBJ)/sdssElColors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/sdssElColors.o : $(PROGS)/sdssElColors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/sdssElColors.o $(PROGS)/sdssElColors.cc 

# CALCULATE SDSS COLORS FOR PICKLES' LIBRARY OF STARS 
$(EXE)/sdssPicklesLibrary : $(OBJ)/sdssPicklesLibrary.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/sdssPicklesLibrary $(OBJ)/sdssPicklesLibrary.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/sdssPicklesLibrary.o : $(PROGS)/sdssPicklesLibrary.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/sdssPicklesLibrary.o $(PROGS)/sdssPicklesLibrary.cc

# SIMULATE OVERDENSITY GRID
$(EXE)/simdensity : $(OBJ)/simdensity.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/simdensity $(OBJ)/simdensity.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simdensity.o : $(PROGS)/simdensity.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/simdensity.o $(PROGS)/simdensity.cc 

# SIMULATE LINE OF SIGHT ABSORBER DISTRIBUTIONS 
$(EXE)/simulateAbsorberLinesOfSight : $(OBJ)/simulateAbsorberLinesOfSight.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/simulateAbsorberLinesOfSight $(OBJ)/simulateAbsorberLinesOfSight.o \
	$(LIBO) $(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateAbsorberLinesOfSight.o : $(PROGS)/simulateAbsorberLinesOfSight.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/simulateAbsorberLinesOfSight.o \
	$(PROGS)/simulateAbsorberLinesOfSight.cc 

# SIMULATE CFHT OBSERVATIONS FROM BASE SIMULATION
$(EXE)/simulateCFHTobs : $(OBJ)/simulateCFHTobs.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/simulateCFHTobs $(OBJ)/simulateCFHTobs.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateCFHTobs.o : $(PROGS)/simulateCFHTobs.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/simulateCFHTobs.o \
	$(PROGS)/simulateCFHTobs.cc 

# SIMULATE LSST OBSERVATIONS FROM BASE SIMULATION
$(EXE)/simulateLSSTobs : $(OBJ)/simulateLSSTobs.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/simulateLSSTobs $(OBJ)/simulateLSSTobs.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateLSSTobs.o : $(PROGS)/simulateLSSTobs.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/simulateLSSTobs.o \
	$(PROGS)/simulateLSSTobs.cc 

# SIMULATE LSST OBSERVATIONS FROM INPUT IMSIM CATALOG OF TRUE MAGS AND REDSHIFTS
$(EXE)/simulateLSSTobsFromTruth : $(OBJ)/simulateLSSTobsFromTruth.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/simulateLSSTobsFromTruth $(OBJ)/simulateLSSTobsFromTruth.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateLSSTobsFromTruth.o : $(PROGS)/simulateLSSTobsFromTruth.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/simulateLSSTobsFromTruth.o \
	$(PROGS)/simulateLSSTobsFromTruth.cc 

###################### BAO PROGRAMS ############################################

# ADD GAUSSIAN Z ERROR TO CATALOG
$(EXE)/addGausszerr : $(OBJ)/addGausszerr.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -std=c++11  -o $(EXE)/addGausszerr $(OBJ)/addGausszerr.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)


$(OBJ)/addGausszerr.o : $(BAOPROGS)/addGausszerr.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/addGausszerr.o \
	$(BAOPROGS)/addGausszerr.cc 

# COMPUTE POWER SPECTRA

$(EXE)/computepsfromgrids : $(OBJ)/computepsfromgrids.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/computepsfromgrids $(OBJ)/computepsfromgrids.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/computepsfromgrids.o : $(BAOPROGS)/computepsfromgrids.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11  -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/computepsfromgrids.o \
	$(BAOPROGS)/computepsfromgrids.cc 

# FIT K BAO
$(EXE)/fitkbao : $(OBJ)/fitkbao.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/fitkbao $(OBJ)/fitkbao.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/fitkbao.o : $(BAOPROGS)/fitkbao.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11  -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/fitkbao.o \
	$(BAOPROGS)/fitkbao.cc

# FIT K BAO WITH BASELINE
$(EXE)/fitkbaobaseline : $(OBJ)/fitkbaobaseline.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/fitkbaobaseline $(OBJ)/fitkbaobaseline.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/fitkbaobaseline.o : $(BAOPROGS)/fitkbaobaseline.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/fitkbaobaseline.o \
	$(BAOPROGS)/fitkbaobaseline.cc

# CALCULATE PHOTO-Z CONVOLUTION FUNCTION
$(EXE)/getpzconvf : $(OBJ)/getpzconvf.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -std=c++11  -o $(EXE)/getpzconvf $(OBJ)/getpzconvf.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/getpzconvf.o : $(BAOPROGS)/getpzconvf.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/getpzconvf.o \
	$(BAOPROGS)/getpzconvf.cc 

# CALCULATE SELECTION FUNCTION OF OBSERVED CATALOG
$(EXE)/getsf : $(OBJ)/getsf.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/getsf $(OBJ)/getsf.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/getsf.o : $(BAOPROGS)/getsf.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/getsf.o \
	$(BAOPROGS)/getsf.cc 

# GRID GALAXY DATA
$(EXE)/grid_data : $(OBJ)/grid_data.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/grid_data $(OBJ)/grid_data.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/grid_data.o : $(BAOPROGS)/grid_data.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/grid_data.o \
	$(BAOPROGS)/grid_data.cc

# SIMULATE RANDOM CATALOGS
# $(EXE)/sim_mcgrids : $(OBJ)/sim_mcgrids.o $(LIBO)
#	mkdir -p $(EXE)
#	mkdir -p $(ROOTOUT)
#	$(CXXLINK) -o $(EXE)/sim_mcgrids $(OBJ)/sim_mcgrids.o $(LIBO) \
#	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

# $(OBJ)/sim_mcgrids.o : $(BAOPROGS)/sim_mcgrids.cc $(LIBH)
#	mkdir -p $(OBJ)
#	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/sim_mcgrids.o \
#	$(BAOPROGS)/sim_mcgrids.cc

# GET DATA SUB GRID
$(EXE)/subfromfull : $(OBJ)/subfromfull.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/subfromfull $(OBJ)/subfromfull.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/subfromfull.o : $(BAOPROGS)/subfromfull.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/subfromfull.o \
	$(BAOPROGS)/subfromfull.cc


###################### TESTING PROGRAMS ########################################


# TEST GOLDENCUT 2D INTERPOLATION
#$(EXE)/testGoldenInterp2D : $(OBJ)/testGoldenInterp2D.o $(LIBO) 
#	mkdir -p $(EXE)
#	mkdir -p $(TESTS)
#	$(CXXLINK) -o $(EXE)/testGoldenInterp2D $(OBJ)/testGoldenInterp2D.o $(LIBO) \
#	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

#$(OBJ)/testGoldenInterp2D.o : $(PROGS)/testGoldenInterp2D.cc $(LIBH)  
#	mkdir -p $(OBJ)
#	$(CXXCOMPILE) -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testGoldenInterp2D.o $(PROGS)/testGoldenInterp2D.cc


# TEST 2D INTERPOLATION
$(EXE)/test2Dinterp : $(OBJ)/test2Dinterp.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/test2Dinterp $(OBJ)/test2Dinterp.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/test2Dinterp.o : $(PROGS)/test2Dinterp.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/test2Dinterp.o $(PROGS)/test2Dinterp.cc

# TEST BASESIM
$(EXE)/testbasesim : $(OBJ)/testbasesim.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testbasesim $(OBJ)/testbasesim.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testbasesim.o : $(PROGS)/testbasesim.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testbasesim.o $(PROGS)/testbasesim.cc

# TEST EXPECTATION-MAXIMIZATION ALGORITHM
$(EXE)/testEMalgorithm : $(OBJ)/testEMalgorithm.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testEMalgorithm $(OBJ)/testEMalgorithm.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testEMalgorithm.o : $(PROGS)/testEMalgorithm.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testEMalgorithm.o $(PROGS)/testEMalgorithm.cc

# TEST LSST ERRORS
$(EXE)/testErrors : $(OBJ)/testErrors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testErrors $(OBJ)/testErrors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testErrors.o : $(PROGS)/testErrors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testErrors.o $(PROGS)/testErrors.cc

# TEST GOODS SIM
$(EXE)/testgoodsmagsim : $(OBJ)/testgoodsmagsim.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testgoodsmagsim $(OBJ)/testgoodsmagsim.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testgoodsmagsim.o : $(PROGS)/testgoodsmagsim.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testgoodsmagsim.o $(PROGS)/testgoodsmagsim.cc

# TEST K CORRECTION 
$(EXE)/testKcorrColors : $(OBJ)/testKcorrColors.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testKcorrColors $(OBJ)/testKcorrColors.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testKcorrColors.o : $(PROGS)/testKcorrColors.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testKcorrColors.o $(PROGS)/testKcorrColors.cc

# TEST K CORRECTION INTERPOLATION
$(EXE)/testKcorrMethod : $(OBJ)/testKcorrMethod.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testKcorrMethod $(OBJ)/testKcorrMethod.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testKcorrMethod.o : $(PROGS)/testKcorrMethod.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testKcorrMethod.o $(PROGS)/testKcorrMethod.cc

# TEST LF 
$(EXE)/testLF : $(OBJ)/testLF.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testLF $(OBJ)/testLF.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testLF.o : $(PROGS)/testLF.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -o $(OBJ)/testLF.o $(PROGS)/testLF.cc 

# TEST Lyman-alpha absorption calculation parts 
$(EXE)/testLymanAlphaAbs : $(OBJ)/testLymanAlphaAbs.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testLymanAlphaAbs $(OBJ)/testLymanAlphaAbs.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testLymanAlphaAbs.o : $(PROGS)/testLymanAlphaAbs.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -o $(OBJ)/testLymanAlphaAbs.o $(PROGS)/testLymanAlphaAbs.cc 

# TEST MADAU
$(EXE)/testMadau : $(OBJ)/testMadau.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testMadau $(OBJ)/testMadau.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testMadau.o : $(PROGS)/testMadau.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testMadau.o $(PROGS)/testMadau.cc

# TEST MEIKSIN
$(EXE)/testMeiksin : $(OBJ)/testMeiksin.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testMeiksin $(OBJ)/testMeiksin.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testMeiksin.o : $(PROGS)/testMeiksin.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testMeiksin.o $(PROGS)/testMeiksin.cc

# TEST SIMULATE DATA USING K CORRECTION TABLES
$(EXE)/testSimReadKcorr : $(OBJ)/testSimReadKcorr.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	mkdir -p $(KCORR)
	$(CXXLINK) -o $(EXE)/testSimReadKcorr $(OBJ)/testSimReadKcorr.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/testSimReadKcorr.o : $(PROGS)/testSimReadKcorr.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testSimReadKcorr.o $(PROGS)/testSimReadKcorr.cc 

# TEST SIMULATE IGM
$(EXE)/testsimulateIGM : $(OBJ)/testsimulateIGM.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testsimulateIGM $(OBJ)/testsimulateIGM.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/testsimulateIGM.o : $(PROGS)/testsimulateIGM.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testsimulateIGM.o $(PROGS)/testsimulateIGM.cc 

# TEST SIMULATION
$(EXE)/testSimulation : $(OBJ)/testSimulation.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(ROOTOUT)
	$(CXXLINK) -o $(EXE)/testSimulation $(OBJ)/testSimulation.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/testSimulation.o : $(PROGS)/testSimulation.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testSimulation.o $(PROGS)/testSimulation.cc

# TEST TEMPLATE FITTING
$(EXE)/testTemplateFitting : $(OBJ)/testTemplateFitting.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testTemplateFitting $(OBJ)/testTemplateFitting.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testTemplateFitting.o : $(PROGS)/testTemplateFitting.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testTemplateFitting.o $(PROGS)/testTemplateFitting.cc

# TEMPORARY TEST CODE
$(EXE)/test : $(OBJ)/test.o $(LIBO)
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/test $(OBJ)/test.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/test.o : $(PROGS)/test.cc $(LIBH)
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/test.o $(PROGS)/test.cc

## programs below here have not been CHECKED or maybe even finished...

# TEST POWER SPECTRUM FROM OVER-DENSITY CUBE
$(EXE)/testpsdenscube : $(OBJ)/testpsdenscube.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testpsdenscube $(OBJ)/testpsdenscube.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testpsdenscube.o : $(PROGS)/testpsdenscube.cc $(LIBH)  
	mkdir -p $(OBJ)
	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testpsdenscube.o $(PROGS)/testpsdenscube.cc 

# TEST DENSITY SIMULATION
$(EXE)/testsimdensity : $(OBJ)/testsimdensity.o $(LIBO) 
	mkdir -p $(EXE)
	mkdir -p $(TESTS)
	$(CXXLINK) -o $(EXE)/testsimdensity $(OBJ)/testsimdensity.o $(LIBO) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testsimdensity.o : $(PROGS)/testsimdensity.cc $(LIBH)  
	mkdir -p $(OBJ)/usr/include/c++/4.6/bits/stl_vector.h:142:9: error: invalid use of incomplete type ???struct SOPHYA::FunRan???

	$(CXXCOMPILE) -std=c++11 -I$(MYCL) -I$(ROOTINC) -o $(OBJ)/testsimdensity.o $(PROGS)/testsimdensity.cc 

############################# UTILITIES ########################################

#$(OBJ)/root_plots.o : $(MYCL)/root_plots.cpp
#	$(CXXCOMPILE) -o $(OBJ)/root_plots.o -I$(ROOTINC) -I$(MYCL) -c $(MYCL)/root_plots.cpp


################################  THE CLASSES ##################################

$(OBJ)/%.o : $(MYCL)/%.cc $(MYCL)/%.h 
	mkdir -p $(OBJ)
	mkdir -p $(TESTS)
	$(CXXCOMPILE) -std=c++11 -I$(ROOTINC) -I$(MYCL) -o $@ $<

#$(OBJ)genefluct3df.o: $(MYCL)genefluct3d.cc $(MYCL)genefluct3d.h
#	$(CXXCOMPILE) -I$(MYCL) -DGEN3D_FLOAT -o $@ $(MYCL)genefluct3d.cc



