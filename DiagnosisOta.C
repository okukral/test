#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
using namespace std;
using namespace RooFit;

// Old version of the macro, not updated
// Use DiagnosisMacro.C instead
//void rf605_profilell()
void DiagnosisOta()
{
	/*
  // C r e a t e   m o d e l   a n d   d a t a s e t
  // -----------------------------------------------

  // Observable
  RooRealVar x("x","x",-20,20) ;

  // Model (intentional strong correlations)
  RooRealVar mean("mean","mean of g1 and g2",-5,-10,10) ;
  RooRealVar sigma_g1("sigma_g1","width of g1",2) ;
  RooGaussian g1("g1","g1",x,mean,sigma_g1) ;

  RooRealVar mean2("mean2", "mean2 of  g2", 0, -10, 10);
  RooRealVar sigma_g2("sigma_g2","width of g2",4,3.0,6.0) ;
  RooGaussian g2("g2","g2",x,mean2,sigma_g2) ;

  RooRealVar frac("frac","frac",0.8,0.0,1.0) ;
  RooAddPdf model("model","model",RooArgList(g1,g2),frac) ;

  // Generate 1000 events
  RooDataSet* data = model.generate(x,1000) ;

  // Plot data and PDF overlaid
  RooPlot* xframe = x.frame(Title("Model and Data"));
  data->plotOn(xframe);
  model.plotOn(xframe);
  model.plotOn(xframe, Components(RooArgSet(g1)), LineColor(kRed), LineStyle(kDashed));
  model.plotOn(xframe, Components(RooArgSet(g2)), LineColor(kRed), LineStyle(kDotted));

  new TCanvas("PlotModel", "PlotModel", 800, 800);
  gPad->SetLeftMargin(0.15);
  xframe->GetYaxis()->SetTitleOffset(1.4);
  xframe->Draw();


  //*/
	//}
	//void nothing(){

	// R e a d   w o r k s p a c e   f r o m   f i l e
	// -----------------------------------------------
	///*
	// Open input file with workspace
	TFile *f = new TFile("FIT_DATA_Psi2SJpsi_PPPrompt_Bkg_SecondOrderChebychev_pt65300_rap016_cent0200_262620_263757.root");

	// Retrieve workspace from file
	RooWorkspace* w = (RooWorkspace*)f->Get("workspace");

	// Retrieve x,model and data from workspace
	RooRealVar* x = w->var("invMass");
	RooAbsPdf* model = w->pdf("pdfMASS_Tot_PP");
	if (model == 0){ model = w->pdf("pdfMASS_Tot_PbPb"); }
	if (model == 0){
		cout << "[ERROR] pdf failed to load from the workspace" << endl;// return false;
	}

	RooAbsData* data = w->data("dOS_DATA_PP");
	if (data == 0){ model = w->pdf("dOS_DATA_PbPb"); }
	if (data == 0){
		cout << "[ERROR] data failed to load from the workspace" << endl; //return false;
	}


	// Print structure of composite p.d.f.
	model->Print("t");

	// F i t   m o d e l   t o   d a t a ,   p l o t   m o d e l 
	// ---------------------------------------------------------

	// Fit model to data
	//model->fitTo(*data);

	// Plot data and PDF overlaid
	RooPlot* xframe = x->frame(Title("J/psi Model and Data"));
	data->plotOn(xframe);
	model->plotOn(xframe);

	// Overlay the background component of model with a dashed line
	//model->plotOn(xframe, Components("bkg"), LineStyle(kDashed));

	// Overlay the background+sig2 components of model with a dotted line
	//model->plotOn(xframe, Components("bkg,sig2"), LineStyle(kDotted));


	RooArgSet* paramSet1 = model->getDependents(data);
	paramSet1->Print("v");  //O: Just check
	RooArgSet* paramSet2 = model->getParameters(data);
	paramSet2->Print("v");
	cout << "Nparams" << paramSet2->getSize()<<endl<<endl;
	TString ParamName;

	TIterator* iter = paramSet2->createIterator();
	TObject* var = iter->Next();
	while (var!=0) {
		ParamName = var->GetName();
		cout << ParamName << endl<<endl;
		var = iter->Next();
	}


	// Draw the frame on the canvas
	new TCanvas("PlotModel", "PlotModel", 1000, 1000);
	gPad->SetLeftMargin(0.15);
	xframe->GetYaxis()->SetTitleOffset(2.0);
	xframe->Draw();
	//*/

//}
//void nothing(){

	// C o n s t r u c t   p l a i n   l i k e l i h o o d
	// ---------------------------------------------------
	/*
	// Construct unbinned likelihood
	RooAbsReal* nll = model.createNLL(*data, NumCPU(2));

	// Minimize likelihood w.r.t all parameters before making plots
	RooMinuit(*nll).migrad();

	// Plot likelihood scan frac 
	RooPlot* frame1 = frac.frame(Bins(100), Range(0.01, 0.95), Title("LL and profileLL in frac"));
	nll->plotOn(frame1, ShiftToZero());

	// Plot likelihood scan in sigma_g2
	RooPlot* frame2 = sigma_g2.frame(Bins(100), Range(3.3, 5.0), Title("LL and profileLL in sigma_g2"));
	nll->plotOn(frame2, ShiftToZero());//*/


	// Construct unbinned likelihood
	RooAbsReal* nll = model->createNLL(*data, NumCPU(4));

	// Minimize likelihood w.r.t all parameters before making plots
	RooMinuit(*nll).migrad();

	// Plot likelihood scan frac 
	RooRealVar* N_Bkg_PP = w->var("N_Bkg_PP");
	RooPlot* frame1 = N_Bkg_PP->frame(Bins(5), Range(39000, 45000), Title("LL and profileLL in N background"));
	nll->plotOn(frame1, ShiftToZero());

	// Plot likelihood scan in sigma_g2
	RooRealVar* N_Jpsi_PP = w->var("N_Jpsi_PP");
	RooPlot* frame2 = N_Jpsi_PP->frame(Bins(5), Range(210000, 220000), Title("LL and profileLL in N Jpsi"));
	nll->plotOn(frame2, ShiftToZero());//*/

	// C o n s t r u c t   p r o f i l e   l i k e l i h o o d
	// -----------------------------------------------------------------------
	/*
	RooAbsReal* pll = nll->createProfile(*N_Bkg_PP);

	cout << "HERE 2" << endl;
	// Plot the profile likelihood in frac
	pll->plotOn(frame1, LineColor(kRed), RooFit::Precision(-1));

	// Adjust frame maximum for visual clarity
	frame1->SetMinimum(0);
	frame1->SetMaximum(20);

	TCanvas *c = new TCanvas("CLikelihoodResult", "CLikelihoodResult", 800, 600);
	c->cd(1); gPad->SetLeftMargin(0.15); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw();

	delete pll;
	delete nll;*/



	// C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   f r a c
	// -----------------------------------------------------------------------

	// The profile likelihood estimator on nll for frac will minimize nll w.r.t
	// all floating parameters except frac for each evaluation
///*	cout << "HERE 1" << endl;

	RooAbsReal* pll_frac = nll->createProfile(*N_Bkg_PP);

	cout << "HERE 2" << endl;
	// Plot the profile likelihood in frac
	pll_frac->plotOn(frame1, LineColor(kRed), RooFit::Precision(-1));

	// Adjust frame maximum for visual clarity
	frame1->SetMinimum(0);
	frame1->SetMaximum(20);
	cout << "HERE 3" << endl;
	//*/
	///*
	
	// C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   s i g m a _ g 2 
	// -------------------------------------------------------------------------------

	// The profile likelihood estimator on nll for sigma_g2 will minimize nll
	// w.r.t all floating parameters except sigma_g2 for each evaluation
	RooAbsReal* pll_sigmag2 = nll->createProfile(*N_Jpsi_PP);

	// Plot the profile likelihood in sigma_g2
	pll_sigmag2->plotOn(frame2, LineColor(kRed), RooFit::Precision(-1));
	//*/

	// Adjust frame maximum for visual clarity
	frame2->SetMinimum(0);
	frame2->SetMaximum(20);



	// Make canvas and draw RooPlots
	TCanvas *c = new TCanvas("rf605_profilell", "rf605_profilell", 1000, 600);
	c->Divide(2);
	c->cd(1); gPad->SetLeftMargin(0.15); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw();
	c->cd(2); gPad->SetLeftMargin(0.15); frame2->GetYaxis()->SetTitleOffset(1.4); frame2->Draw();

	delete pll_frac;
	delete pll_sigmag2;
	delete nll;
}


