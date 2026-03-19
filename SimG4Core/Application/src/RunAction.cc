#include "SimG4Core/Application/interface/RunAction.h"
#include "SimG4Core/Application/interface/SimRunInterface.h"

#include "SimG4Core/Notification/interface/BeginOfRun.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <fstream>

#include "G4AnalysisManager.hh"

RunAction::RunAction(const edm::ParameterSet& p, SimRunInterface* rm, bool)
    : m_runInterface(rm), m_stopFile(p.getParameter<std::string>("StopFile")), m_outputFile(p.getParameter<std::string>("outputFile")) {
                  auto analysisManager = G4AnalysisManager::Instance();
    }

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run* aRun) {
  if (!m_stopFile.empty() && std::ifstream(m_stopFile.c_str())) {
    edm::LogWarning("SimG4CoreApplication") << "RunAction::BeginOfRunAction: termination signal received";
    m_runInterface->abortRun(true);
  }
  BeginOfRun r(aRun);
  m_beginOfRunSignal(&r);

  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root"); // set in macrofile
  analysisManager->SetVerboseLevel(10);
  analysisManager->SetNtupleMerging(true);  // important for MT

  analysisManager->CreateNtuple("tree", "tree for HCAL 2006 TB experiment");

  analysisManager->CreateNtupleDColumn("ECAL_eresponse");

  analysisManager->CreateNtupleDColumn("HCAL_eresponse");

  analysisManager->CreateNtupleDColumn("ECAL_eresponse_raw");

  analysisManager->CreateNtupleDColumn("HCAL_eresponse_raw");

  analysisManager->FinishNtuple();
  analysisManager->SetFileName(m_outputFile);
  analysisManager->OpenFile(); // name set in macrofile

}

void RunAction::EndOfRunAction(const G4Run* aRun) {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  EndOfRun r(aRun);
  m_endOfRunSignal(&r);
}
