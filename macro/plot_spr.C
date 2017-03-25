//
// Generates outfile with TTree for
// SPR, mean, NPH, and PID for
// reconstructed thetaC from
// EIC simulation
//
// Set avr = 1 to use averaged
// LUT reconstruction
//
// Author: Lee Allison, 2017

#include "TCanvas.h"
#include "TCollection.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TList.h"
#include "TROOT.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TTree.h"
#include "TPaveStats.h"

void formatStats(TCanvas *&canv, TPaveStats *&pavestats, TH1D *&hist,
				 double track, double mean, double spr)
{
	// clean up  pavestats
	pavestats->Clear();
	
	// define bounds of stats box
	// numbers are percentage of full pad
	double x1 = 0.73;  pavestats->SetX1NDC(x1); // left
	double x2 = 0.95;  pavestats->SetX2NDC(x2); // right
	double y1 = 0.77;  pavestats->SetY1NDC(y1); // bottom
	double y2 = 0.95;  pavestats->SetY2NDC(y2); // top

	pavestats->AddText(Form("%.2f#circ track",track));
	pavestats->AddText(Form("entries = %d",hist->GetEntries());
	pavestats->AddText(Form("mean = %.3f mrad",mean));
	pavestats->AddText(Form("SPR = %.3f mrad",spr));
	pavestats->SetTextSize(0.03);

	// remove default stats and draw new stats box
	hist->SetStats(0);
	pavestats->Draw();
	canv->Modified();

	return;
}

void plot_spr( bool averaged = 0 )
{
	gStyle->SetOptFit(1); // show fitting parameters
	gErrorIgnoreLevel = kWarning; // ignore 'Info in..' messages
	
	// save path
	TString avr = averaged ? "_avr" : "";
	TString savepath = "./analysis/";

	// want to print canvases?
	bool printCanv  = true;

	// variables for output trees
	int pid(0);
	double track(0); // polar track angle
	double nph(0), spr(0), mean(0);

	const int bins = 120;
	int pids[] = {211, 321};

	// create new file and trees to hold analysis results
	TFile *out_file = new TFile(savepath+Form("analysis%s.root",avr.Data()),"recreate");
	TTree *out_tree = new TTree("analysis","analysis tree");

	// make branches for simulation tree
	out_tree->Branch("track",&track,"track/D");
	out_tree->Branch("pid",&pid,"pid/I");
	out_tree->Branch("nph",&nph,"nph/D");
	out_tree->Branch("spr",&spr,"spr/D");
	out_tree->Branch("mean",&mean,"mean/D");

	// define paths to data and simulation files
	// assumes I'm only looking at 3CS lens
	TString simdir;
	if(averaged)
		simdir = "../simulation/reco/avr/";
	else
		simdir = "../simulation/reco/";
	
	// grab files from data directory
	TSystemDirectory dir(simdir,simdir);
	TList *files = dir.GetListOfFiles();
	TSystemFile *file;
	TString fname;
	TIter next(files);

	// reused objects in loop
	TCanvas *canv = new TCanvas(); // canvas for drawing
	TPaveStats *pavestats = new TPaveStats(); // hist stats
	pavestats->SetBorderSize(1);
	pavestats->SetFillColor(kWhite);

	TFile *sim_file;  // for data and simulation
	TTree *dirc, *reco;

	// loop over each root file in directory
	// and extract nph, spr, and mean thetaC
	// for both data and sim w/ and w/o bgsub
	int n(0); // counter
	while( file=(TSystemFile*)next() )
	{
		fname = file->GetName();
		
		if( file->IsDirectory() || !fname.EndsWith("_spr.root") )
			continue; // don't process directories or non-root files

		// define data file and trees
		sim_file = TFile::Open(simdir+fname);
		dirc = (TTree*)sim_file->Get("dirc");
		reco = (TTree*)sim_file->Get("reco");
		dirc->SetBranchAddress("theta",&track);
		dirc->SetBranchAddress("nph",&nph);
		dirc->GetEntry(0);

		//if(track != 75) continue;

		cout << "Processing " << track << ": " << fname << endl;
		
		// do the same for kaons and pions
		for(int i=0; i<2; i++) //pid=211; pid<=321; pid+=110)
		{
			pid = pids[i];
			TString particle;
			if(pid==211){ particle = "pi+"; } 
			else if(pid==321){ particle = "K+"; }
			else{ particle == "proton"; }

			// define histograms for reconstructed thetaC
			TH1D *theta = new TH1D("theta","reconstructed theta",bins,0.7,0.9);

			// define fit functions
			TString fitfunc = "gaus(0)+pol0(3)";
			TF1 *fit  = new TF1("fit",fitfunc);

			// set time cut to be +- 500ps
			double time(0.5);
			TString cut = Form("PID==%d && abs(diff-%f)<1",pid,time);

			// project thetaC onto histograms
			reco->Project("theta","theta",cut);

			// set initial values for fit and
			// fit reconstructed thetaC
			mean = theta->GetBinCenter(theta->GetMaximumBin());
			fit->SetParameters(theta->GetMaximum(),
							   mean,
							   0.01,
							   theta->GetBinContent(1));
			fit->SetParName(1,"#theta_{C}");
			fit->SetParName(2,"#sigma");
			theta->Fit("fit","Q","",0.78,0.88);

			// set mean and SPR values
			mean     = 1000*fit->GetParameter(1);
			spr      = 1000*fit->GetParameter(2);
		
			theta->SetTitle(Form("sim #theta_{C} %s reco, %.2f#circ",particle.Data(),track));
			theta->GetYaxis()->SetRangeUser(0,1.1*theta->GetMaximum());
			formatStats(canv,pavestats,theta,track,mean,spr);
			canv->Print(savepath+Form("thetaC%s/theta_%s_%.2f.png",avr.Data(),particle.Data(),track));

			// write TTrees
			out_tree->Fill();

			// clear persistant objects and increment
			// counter for next loop
			canv->Clear();
			sim_file->Clear();
			reco->Clear();
			dirc->Clear();

			// clean up pointers
			delete fit;
			delete theta;
		}
		n++;
		cout << "finishing track " << track << endl;
		//if(n>0) break;
	}

	//cout << "out of the loop" << endl;
	// write output file
	out_file->Write();
	cout << "out_file written" << endl;
}

