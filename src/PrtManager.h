// -----------------------------------------
// PrtManager.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtManager_h
#define PrtManager_h

#include "globals.hh"

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "TClonesArray.h"

#include "PrtEvent.h"
#include "PrtHit.h"
#include "PrtTrackInfo.h"

class PrtManager
{
  static PrtManager* fInstance;
  TFile *fRootFile;
  TTree *fTree;
  PrtEvent *fEvent;
  PrtTrackInfo *fTrackInfo;
  PrtHit *fHit;
  TH1F *fHist;

public:
  PrtManager(G4String outfile, G4int runtype);
  ~PrtManager(){};
  static PrtManager* Instance(G4String outfile="hits.root", G4int runtype=0);
  void Save()             { fRootFile->Write(); }
  void Fill();
  void FillLut();
  void AddEvent(PrtEvent event);
  void AddHit(PrtHit hit);
  void AddTrackInfo(PrtTrackInfo trackinfo);
  PrtEvent* Event(){ return fEvent; }
  
  // Mutators
  void SetRunType(int val){ fRunType = val; }
  void SetPhysList(int val){ fPhysList = val; }
  void SetGeometry(int val){ fGeometry = val; }
  void SetBeamDimension(int val){ fBeamDimension = val; }
  void SetRadiator(int val){ fRadiator = val; }
  void SetLens(int val){ fLens = val; }
  void SetMcpLayout(int val){ fMcpLayout = val; }
  void SetAngle(double val){ fAngle = val; }
  void SetParticle(int val){ fParticle = val; }
  void SetMomentum(TVector3 val){ fMomentum = val; }
  void SetCurrentCherenkov(double val){ fCurrentCherenkov = val; }
  void SetShift(double val){ fShift = val; }
  void SetTest(double val){ fTest = val; }

  // Accessors
  int GetRunType(){ return fRunType; }
  int GetPhysList(){ return fPhysList; }
  int GetGeometry(){ return fGeometry; }
  int GetBeamDinsion(){ return fBeamDimension; }
  int GetRadiator(){ return fRadiator; }
  int GetLens(){ return fLens; }
  int GetMcpLayout(){ return fMcpLayout; }
  double GetAngle(){ return fAngle; }
  int GetParticle(){ return fParticle; }
  TVector3 GetMomentum(){ return fMomentum; }
  double GetCurrentCherenkov(){ return fCurrentCherenkov; }
  double GetShift(){ return fShift; }
  double GetTest(){ return fTest; }
  
private: 
  int fRunType;
  int fPhysList;
  int fGeometry;
  int fRadiator;
  int fLens;
  int fMcpLayout;
  double fAngle;
  int fParticle;
  int fBeamDimension;
  TVector3 fMomentum;
  TClonesArray *fLut;
  TClonesArray *fTrackInfoArray;
  double fCurrentCherenkov;
  double fShift;
  double fTest;
};

#endif
