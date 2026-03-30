#ifndef SimG4Core_TrackingAction_H
#define SimG4Core_TrackingAction_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimG4Core/Notification/interface/SimActivityRegistry.h"

#include "G4UserTrackingAction.hh"
#include "G4Region.hh"
#include "G4Track.hh"

#include <vector>
#include <map>
#include <string>
#include <unordered_map>

class SimTrackManager;
class TrackWithHistory;
class BeginOfTrack;
class EndOfTrack;
class CMSSteppingVerbose;
class TrackInformation;
class G4ParticleDefinition;

class TrackingAction : public G4UserTrackingAction {
public:
  explicit TrackingAction(SimTrackManager*, CMSSteppingVerbose*, const edm::ParameterSet& ps);
  ~TrackingAction() override = default;

  void PreUserTrackingAction(const G4Track* aTrack) override;
  void PostUserTrackingAction(const G4Track* aTrack) override;

  inline TrackWithHistory* currentTrackWithHistory() { return currentTrack_; }
  inline const G4Track* geant4Track() const { return g4Track_; }
  inline G4TrackingManager* getTrackManager() { return fpTrackingManager; }

  SimActivityRegistry::BeginOfTrackSignal m_beginOfTrackSignal;
  SimActivityRegistry::EndOfTrackSignal m_endOfTrackSignal;
    std::map<int, std::pair<const G4ParticleDefinition*,double>> trackIDmap;

private:
  SimTrackManager* trackManager_;
  CMSSteppingVerbose* steppingVerbose_;
  const G4Track* g4Track_ = nullptr;
  TrackInformation* trkInfo_ = nullptr;
  TrackWithHistory* currentTrack_ = nullptr;
  int endPrintTrackID_;
  bool checkTrack_;
  bool doFineCalo_;
  bool saveCaloBoundaryInformation_;
  double ekinMin_;
  std::vector<double> ekinMinRegion_;
  std::vector<G4Region*> ptrRegion_;
  int PDG_OTHERS = 0;
  std::map<int,int> pdgID_to_histoID_map = { {11,0},
                                             {22,5},
                                             {2112,10},
                                             {211,15},
                                             {-211,20},
                                             {111,25},
                                             {2212,30},
                                             {PDG_OTHERS,35}
                                            };
  double GetParticleHistoID(const G4Track * track) const {
      int pdgID = track->GetParticleDefinition()->GetPDGEncoding();
      auto particleInformation = pdgID_to_histoID_map.find(pdgID);
      if( pdgID_to_histoID_map.end() == particleInformation )
        return pdgID_to_histoID_map.at(PDG_OTHERS);
      else
        return particleInformation->second;
    }
    std::unordered_map<std::string,int> fProcNameId = {{"muBrems", 1},
                                                        {"muPairProd", 2},
                                                        {"alphaInelastic", 3},
                                                        {"ionInelastic", 4},
                                                        {"positronNuclear", 5},
                                                        {"anti_lambdaInelastic", 6},
                                                        {"He3Inelastic", 7},
                                                        {"sigma+Inelastic", 8},
                                                        {"tInelastic", 9},
                                                        {"anti_sigma-Inelastic", 10},
                                                        {"electronNuclear", 11},
                                                        {"anti_neutronInelastic", 12},
                                                        {"hIoni", 13},
                                                        {"nCapture", 14},
                                                        {"pi+Inelastic", 15},
                                                        {"sigma-Inelastic", 16},
                                                        {"hBertiniCaptureAtRest", 17},
                                                        {"anti_sigma+Inelastic", 18},
                                                        {"muMinusCaptureAtRest", 19},
                                                        {"eIoni", 20},
                                                        {"eBrem", 21},
                                                        {"Decay", 22},
                                                        {"phot", 23},
                                                        {"compt", 24},
                                                        {"hFritiofCaptureAtRest", 25},
                                                        {"annihil", 26},
                                                        {"muIoni", 27},
                                                        {"neutronInelastic", 28},
                                                        {"hadElastic", 29},
                                                        {"protonInelastic", 30},
                                                        {"photonNuclear", 31},
                                                        {"CoulombScat", 32},
                                                        {"anti_protonInelastic", 33},
                                                        {"kaon0LInelastic", 34},
                                                        {"hPairProd", 35},
                                                        {"kaon+Inelastic", 36},
                                                        {"conv", 37},
                                                        {"dInelastic", 38},
                                                        {"lambdaInelastic", 39},
                                                        {"kaon0SInelastic", 40},
                                                        {"kaon-Inelastic", 41},
                                                        {"pi-Inelastic", 42},
                                                        {"ionIoni", 43},
                                                        {"hBrems", 44},
                                                    };
};
#endif
