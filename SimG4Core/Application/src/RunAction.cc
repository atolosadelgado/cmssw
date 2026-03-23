#include "SimG4Core/Application/interface/RunAction.h"
#include "SimG4Core/Application/interface/SimRunInterface.h"

#include "SimG4Core/Notification/interface/BeginOfRun.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <fstream>

#include "G4AnalysisManager.hh"
#include "G4PhysicsModelCatalog.hh"

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
  // create histograms of initial and final energy, and lifetime of particles
  {
      int nmodels = G4PhysicsModelCatalog::Entries();
      std::vector<std::pair<int, G4String>> particles = {
            {11,   "electron"},   // e-
            {22,   "gamma"},
            {2112, "neutron"},
            {211,  "piPlus"},
            {-211, "piMinus"},
            {111,  "pi0"},
            {2212, "proton"},
            {0, "others"}

      };
      for (const auto& [pdg, name] : particles) {
            analysisManager->CreateH2(
                "hE0_" + name, "", 2500, -15, 10, nmodels, 0, nmodels);

            analysisManager->CreateH2(
                "hEf_" + name, "", 2500, -15, 10, nmodels, 0, nmodels);

            analysisManager->CreateH2(
                "hTf_" + name, "", 2500, -15, 10, nmodels, 0, nmodels);

      }
  } // end create histograms
}

void RunAction::EndOfRunAction(const G4Run* aRun) {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  EndOfRun r(aRun);
  m_endOfRunSignal(&r);
}
