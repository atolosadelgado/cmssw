#include "SimG4Core/Application/interface/RunAction.h"
#include "SimG4Core/Application/interface/SimRunInterface.h"

#include "SimG4Core/Notification/interface/BeginOfRun.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <fstream>

#include "TH1D.h"
#include "TFile.h"

RunAction::RunAction(const edm::ParameterSet& p, SimRunInterface* rm, bool master)
    : m_runInterface(rm), m_stopFile(p.getParameter<std::string>("StopFile")), fIsMaster(master) {
      hHGCal_eprofilez = std::make_unique<TH1D>("hHGCal_eprofilez",
                                                "Energy profile;z (mm);E (MeV)",
                                                nbins_zprofile,
                                                zmin_zprofile,
                                                zmax_zprofile
                                                );
      hHGCal_eprofilez->SetDirectory(0);
    }

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run* aRun) {
  if (!m_stopFile.empty() && std::ifstream(m_stopFile.c_str())) {
    edm::LogWarning("SimG4CoreApplication") << "RunAction::BeginOfRunAction: termination signal received";
    m_runInterface->abortRun(true);
  }
  BeginOfRun r(aRun);
  m_beginOfRunSignal(&r);
}

void RunAction::EndOfRunAction(const G4Run* aRun) {
  EndOfRun r(aRun);
  m_endOfRunSignal(&r);
  // if is master, do nothing else
  if(fIsMaster) return;

  int ThreadIndex = m_runInterface->getThreadIndex();
  std::string ofilename = "test" + std::to_string(ThreadIndex) + ".root";
  TFile * ofile = new TFile( ofilename.c_str() ,"recreate");
  hHGCal_eprofilez->Write();
  ofile->Close();

}

void RunAction::FillHGCalEprofilez(double zabs_mm, double edep_MeV){hHGCal_eprofilez->Fill(zabs_mm, edep_MeV);}
