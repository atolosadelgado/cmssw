#ifndef SimG4Core_RunAction_H
#define SimG4Core_RunAction_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimG4Core/Notification/interface/SimActivityRegistry.h"

#include "G4UserRunAction.hh"

#include <string>
#include <memory>

#include "SimG4Core/Application/interface/SecondaryEscapeCounter.h"

class SimRunInterface;
class BeginOfRun;
class EndOfRun;
class TH1D;

class RunAction : public G4UserRunAction {
public:
  explicit RunAction(const edm::ParameterSet& ps, SimRunInterface*, bool master);
  ~RunAction() override;

  void BeginOfRunAction(const G4Run* aRun) override;
  void EndOfRunAction(const G4Run* aRun) override;

  void FillHGCalEprofilez(double zabs_mm, double edep_MeV);

  SimActivityRegistry::BeginOfRunSignal m_beginOfRunSignal;
  SimActivityRegistry::EndOfRunSignal m_endOfRunSignal;

  void RegisterCreationOfParticle(const G4Track* track)
  {
    totalSecondaryCounter.RegisterCreation(track);
    gammaSecondaryCounter.RegisterCreation(track);
    electronSecondaryCounter.RegisterCreation(track);
  }
  void RegisterEndOfParticle(const G4Track* track)
  {
    totalSecondaryCounter.RegisterEnd(track);
    gammaSecondaryCounter.RegisterEnd(track);
    electronSecondaryCounter.RegisterEnd(track);
  }

  void SecondaryCounterFillHistogramAndReset(){
    totalSecondaryCounter.FillHistogramsAndResetCounters();
    gammaSecondaryCounter.FillHistogramsAndResetCounters();
    electronSecondaryCounter.FillHistogramsAndResetCounters();
  }

private:
  SimRunInterface* m_runInterface;
  std::string m_stopFile;
  int nbins_zprofile = 180000;
  double zmin_zprofile = 3000;
  double zmax_zprofile = 6000;
  std::unique_ptr<TH1D> hHGCal_eprofilez;
  bool fIsMaster = {false};
  SecondaryEscapeCounter totalSecondaryCounter = {""};
  SecondaryEscapeCounter gammaSecondaryCounter = {"gamma"};
  SecondaryEscapeCounter electronSecondaryCounter = {"e-"};

};

#endif
