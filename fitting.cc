// The following contains the code for performing an amplitude fit to the input data. In this file you will find examples of different
// lineshapes and a model independent description of the D_0^{*0} component 


#include <cstdlib>
#include <iostream>
#include <vector>


#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TTree.h"

#include "LauSimpleFitModel.hh"
#include "LauBkgndDPModel.hh"
#include "LauDaughters.hh"
#include "LauEffModel.hh"
#include "LauIsobarDynamics.hh"
#include "LauMagPhaseCoeffSet.hh"
#include "LauModIndPartWaveMagPhase.hh"
#include <LauModIndPartWaveRealImag.hh>
#include "LauRealImagCoeffSet.hh"
#include "LauResonanceMaker.hh"
#include "LauVetoes.hh"
#include "LauParameter.hh"
#include "LauResonanceInfo.hh"
#include "LauComplex.hh"
#include "LauAbsModIndPartWave.hh"
#include "Lau1DCubicSpline.hh"

using namespace std;



void usage( std::ostream& out, const TString& progName )
{
	out<<"Usage:\n";
	out<<progName<<" gen [nExpt = 1] [firstExpt = 0]\n";
	out<<"or\n";
	out<<progName<<" fit <iFit> [nExpt = 1] [firstExpt = 0]"<<std::endl;
}



int main( int argc, char** argv )
{
	// Process command-line arguments
	// Usage:
	// ./GenFit3pi gen [nExpt = 1] [firstExpt = 0]
	// or
	// ./GenFit3pi fit <iFit> [nExpt = 1] [firstExpt = 0]
	if ( argc < 2 ) {
		usage( std::cerr, argv[0] );
		return EXIT_FAILURE;
	}

	TString command = argv[1];
	command.ToLower();
	Int_t iFit(0);
	Int_t nExpt(1);
	Int_t firstExpt(0);
	if ( command == "gen" ) {
		if ( argc > 2 ) {
			nExpt = atoi( argv[2] );
			if ( argc > 3 ) {
				firstExpt = atoi( argv[3] );
			}
		}
	} else if ( command == "fit" ) {
		if ( argc < 3 ) {
			usage( std::cerr, argv[0] );
			return EXIT_FAILURE;
		}
		iFit = atoi( argv[2] );
		if ( argc > 3 ) {
			nExpt = atoi( argv[3] );
			if ( argc > 4 ) {
				firstExpt = atoi( argv[4] );
			}
		}
	} else {
		usage( std::cerr, argv[0] );
		return EXIT_FAILURE;
	}
	Bool_t squareDP = kTRUE;

	LauDaughters* daughters = new LauDaughters("B-", "pi-", "pi-", "D+", squareDP);
	LauVetoes* vetoes = new LauVetoes();

	LauEffModel* effModel = new LauEffModel(daughters, vetoes);

	TFile *effHistFile = TFile::Open("efficiency.root", "read");
	TH2* effHist = dynamic_cast<TH2*>(effHistFile->Get("h_ratio;1"));
	Bool_t useInterpolation = kTRUE;
	Bool_t fluctuateBins = kFALSE;
	Bool_t useUpperHalf = kTRUE;
	effModel->setEffHisto(effHist, useInterpolation, fluctuateBins, 0.0, 0.0, useUpperHalf, squareDP);

	// Set the values of the Blatt-Weisskopf barrier radii and whether they are fixed or floating
	LauResonanceMaker& resMaker = LauResonanceMaker::get();
	resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Parent,     4.0 );
	resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Charm,      4.0 );
	resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Light,      4.0 );
	resMaker.setDefaultBWRadius( LauBlattWeisskopfFactor::Beauty,     4.0 );
	resMaker.fixBWRadius( LauBlattWeisskopfFactor::Parent,  kTRUE);
	resMaker.fixBWRadius( LauBlattWeisskopfFactor::Charm,   kTRUE);
	resMaker.fixBWRadius( LauBlattWeisskopfFactor::Light,   kTRUE);
	resMaker.fixBWRadius( LauBlattWeisskopfFactor::Beauty,  kTRUE);

	LauIsobarDynamics* sigModel = new LauIsobarDynamics(daughters, effModel);

	LauAbsResonance* reson(0);

	reson = sigModel->addResonance("D*0_2",  2, LauAbsResonance::RelBW);
	reson->changeResonance(2.4607, 0.0475, 2);
    reson->fixMass(kFALSE);
	reson->fixWidth(kFALSE);
	LauParameter* mass = reson->getMassPar();
	LauParameter* width = reson->getWidthPar();
	mass->addGaussianConstraint(2.4607, 0.7e-3); 
	width->addGaussianConstraint(0.0475, 0.8e-3);

	reson = sigModel->addResonance("D*0_0",	2, LauAbsResonance::MIPW_MagPhase);
	LauModIndPartWaveMagPhase* mipw = dynamic_cast<LauModIndPartWaveMagPhase*>(reson);
	if (mipw==nullptr){
		std::cout << "MIPW pointer is null" << std::endl;
		return 0;
	}

	//vector of knot masses - ignore ends as they are dealt with internally
	std::set<Double_t> knot_mass{
	2.10,
	2.20,
	2.30,
	2.40,
	2.50,
	2.60,
	2.70,
	2.80,
	2.90,
	3.10,
	4.10,
	};

	mipw->defineKnots(knot_mass);
	mipw->floatKnotsSecondStage(kFALSE);

	//Set magnitude and phase for the knots (including the end points here) 
	mipw->setKnotAmp(0,   0.12,  -2.82,kFALSE,kFALSE);
	mipw->setKnotAmp(1,   0.58,  -1.56,kFALSE,kFALSE);
	mipw->setKnotAmp(2,   0.73,  -1.00,kFALSE,kFALSE);
	mipw->setKnotAmp(4,   0.68,  -0.42,kFALSE,kFALSE);
	mipw->setKnotAmp(3,   0.5,    0.0,kTRUE,kTRUE);
	mipw->setKnotAmp(5,   0.23,  -0.00,kFALSE,kFALSE);
	mipw->setKnotAmp(6,   0.23,  -0.42,kFALSE,kFALSE);
	mipw->setKnotAmp(7,   0.15,  -0.31,kFALSE,kFALSE);
	mipw->setKnotAmp(8,   0.17,  -0.63,kFALSE,kFALSE);
	mipw->setKnotAmp(9,   0.20,  -0.87,kFALSE,kFALSE);
	mipw->setKnotAmp(10,  0.14,  -1.16,kFALSE,kFALSE);
	mipw->setKnotAmp(11,  0.08,   1.02,kFALSE,kFALSE);
	mipw->setKnotAmp(12,  0.0,  0.0,kTRUE, kTRUE);


    // Rel. BW description of D*0_0
	// reson = sigModel->addResonance("D*0_0",  2, LauAbsResonance::RelBW);
	// reson->changeResonance(2.3, 0.27, 0);
    // reson->fixMass(kTRUE);
	// reson->fixWidth(kTRUE);
	// LauParameter* mass_1 = reson->getMassPar();
	// LauParameter* width_1 = reson->getWidthPar();
	// mass_1->addGaussianConstraint(2.3, 0.019); 
	// width_1->addGaussianConstraint(0.27, 0.004);

	reson = sigModel->addResonance("NonReson", 0, LauAbsResonance::FlatNR);

	sigModel->setASqMaxValue(0.5);


	LauSimpleFitModel* fitModel = new LauSimpleFitModel(sigModel);
	std::vector<LauAbsCoeffSet*> coeffset;
	coeffset.push_back ( new LauRealImagCoeffSet("D*0_0",  -5.68103e-01, -4.04082e-01, kFALSE, kFALSE)) ;
	coeffset.push_back ( new LauRealImagCoeffSet("D*0_2", 1.0, 0.0, kTRUE, kTRUE)) ;
	coeffset.push_back ( new LauRealImagCoeffSet("NonReson",-1.06978e+00 , -1.92616e+00, kFALSE, kFALSE)) ;
	for (std::vector<LauAbsCoeffSet*>::iterator iter=coeffset.begin(); iter!=coeffset.end(); ++iter) {
		fitModel->setAmpCoeffSet(*iter);
	}

	// Set the signal yield and define whether it is fixed or floated
	const Double_t nSigEvents = 12551.0;
	LauParameter* signalEvents = new LauParameter("signalEvents", nSigEvents, -2.0*nSigEvents, 2.0*nSigEvents, kTRUE);
    fitModel->setNSigEvents(signalEvents);

	// Set the number of experiments to generate or fit and which
	// experiment to start with
	fitModel->setNExpts( nExpt, firstExpt );

	// Switch on/off calculation of asymmetric errors.
	fitModel->useAsymmFitErrors(kFALSE);

	fitModel->twoStageFit(kFALSE);

	fitModel->doEMLFit(kTRUE);

	fitModel->doPoissonSmearing(kTRUE);

	fitModel->useRandomInitFitPars(kFALSE);

	fitModel->doSFit("evtWeight", 1);

	TString fitToyFileName("fittoy");
	fitToyFileName += iFit;
	fitToyFileName += ".root";
	fitModel->compareFitData(1, fitToyFileName);

	// TString splotFileName("splot_");
	// splotFileName += iFit;
	// splotFileName += ".root";
	// fitModel->writeSPlotData(splotFileName, "splot", kFALSE);

	TString dataFile("input-root-file.root");
	TString treeName("tree");
	TString rootFileName("");
	TString tableFileName("");
	if (command == "fit") {
		rootFileName = "newtrial_"; rootFileName += iFit;
		rootFileName += "_expt_"; rootFileName += firstExpt;
		rootFileName += "-"; rootFileName += (firstExpt+nExpt-1);
		rootFileName += ".root";
		tableFileName = "newtrial_tab"; tableFileName += iFit;
	} else {
		rootFileName = "dummy.root";
		tableFileName = "gen3piResults";
	}

	// Execute the generation/fit
	fitModel->run( command, dataFile, treeName, rootFileName, tableFileName );

	return EXIT_SUCCESS;
}