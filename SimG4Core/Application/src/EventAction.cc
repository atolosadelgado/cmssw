#include "SimG4Core/Application/interface/EventAction.h"
#include "SimG4Core/Application/interface/SimRunInterface.h"
#include "SimG4Core/Notification/interface/TmpSimEvent.h"
#include "SimG4Core/Notification/interface/TmpSimVertex.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Notification/interface/CMSSteppingVerbose.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Randomize.hh"

#include "G4AnalysisManager.hh"


EventAction::EventAction(const edm::ParameterSet& p,
                         SimRunInterface* rm,
                         SimTrackManager* iManager,
                         CMSSteppingVerbose* sv)
    : m_runInterface(rm),
      m_trackManager(iManager),
      m_SteppingVerbose(sv),
      m_stopFile(p.getParameter<std::string>("StopFile")),
      m_printRandom(p.getParameter<bool>("PrintRandomSeed")),
      m_debug(p.getUntrackedParameter<bool>("debug", false)) {}

void EventAction::BeginOfEventAction(const G4Event* anEvent) {
  atd_ecal_energy = 0;
  atd_ecal_energy_raw = 0;
  atd_hcal_energy = 0;
  atd_hcal_energy_raw = 0;
  BeginOfEvent e(anEvent);
  m_beginOfEventSignal(&e);

  if (m_printRandom) {
    edm::LogVerbatim("SimG4CoreApplication")
        << "BeginOfEvent " << anEvent->GetEventID() << " Random number: " << G4UniformRand();
  }

  if (nullptr != m_SteppingVerbose) {
    m_SteppingVerbose->beginOfEvent(anEvent);
  }
}

void EventAction::EndOfEventAction(const G4Event* anEvent) {
  // std::cout << "EventACtion::EndOfEventAction: " << anEvent->GetEventID() << std::endl;
  // std::cout << "\t atd_ecal_energy : " <<  atd_ecal_energy << std::endl;
  // std::cout << "\t atd_hcal_energy : " <<  atd_hcal_energy << std::endl;
  // std::cout << "\t atd_ecal_energy_raw : " <<  atd_ecal_energy_raw << std::endl;
  // std::cout << "\t atd_hcal_energy_raw : " <<  atd_hcal_energy_raw << std::endl;
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(0, atd_ecal_energy);
  analysisManager->FillNtupleDColumn(1, atd_hcal_energy);
  analysisManager->FillNtupleDColumn(2, atd_ecal_energy_raw );
  analysisManager->FillNtupleDColumn(3, atd_hcal_energy_raw );
  analysisManager->AddNtupleRow();
  if (m_printRandom) {
    edm::LogVerbatim("SimG4CoreApplication")
        << "EventACtion::EndOfEventAction: " << anEvent->GetEventID() << " Random number: " << G4UniformRand();
  }
  if (!m_stopFile.empty() && std::ifstream(m_stopFile.c_str())) {
    edm::LogWarning("SimG4CoreApplication")
        << "EventACtion::EndOfEventAction: termination signal received at event " << anEvent->GetEventID();
    // soft abort run
    m_runInterface->abortRun(true);
  }
  if (anEvent->GetNumberOfPrimaryVertex() == 0) {
    edm::LogWarning("SimG4CoreApplication") << "EventACtion::EndOfEventAction: event " << anEvent->GetEventID()
                                            << " must have failed (no G4PrimaryVertices found) and will be skipped";
    return;
  }

  m_trackManager->storeTracks();

  // dispatch now end of event
  EndOfEvent e(anEvent);
  m_endOfEventSignal(&e);

  // delete transient objects
  m_trackManager->reset();
}

void EventAction::abortEvent() { m_runInterface->abortEvent(); }
