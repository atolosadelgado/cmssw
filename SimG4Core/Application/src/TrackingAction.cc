#include "SimG4Core/Application/interface/TrackingAction.h"

#include "SimG4Core/Notification/interface/CurrentG4Track.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfTrack.h"
#include "SimG4Core/Notification/interface/TrackInformation.h"
#include "SimG4Core/Notification/interface/TrackWithHistory.h"
#include "SimG4Core/Notification/interface/SimTrackManager.h"
#include "SimG4Core/Notification/interface/CMSSteppingVerbose.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4UImanager.hh"
#include "G4AnalysisManager.hh"
#include "G4TrackingManager.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4MaterialCutsCouple.hh"
//#define EDM_ML_DEBUG
#include "G4RegionStore.hh"
#include "G4Threading.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"

TrackingAction::TrackingAction(SimTrackManager* stm, CMSSteppingVerbose* sv, const edm::ParameterSet& p)
    : trackManager_(stm),
      steppingVerbose_(sv),
      endPrintTrackID_(p.getParameter<int>("EndPrintTrackID")),
      checkTrack_(p.getUntrackedParameter<bool>("CheckTrack", false)),
      doFineCalo_(p.getParameter<bool>("DoFineCalo")),
      saveCaloBoundaryInformation_(p.getParameter<bool>("SaveCaloBoundaryInformation")),
      ekinMin_(p.getParameter<double>("PersistencyEmin") * CLHEP::GeV),
      ekinMinRegion_(p.getParameter<std::vector<double>>("RegionEmin")) {
  double eth = p.getParameter<double>("EminFineTrack") * CLHEP::MeV;
  if (doFineCalo_ && eth < ekinMin_) {
    ekinMin_ = eth;
  }
  edm::LogVerbatim("SimG4CoreApplication") << "TrackingAction: boundary: " << saveCaloBoundaryInformation_
                                           << "; DoFineCalo: " << doFineCalo_ << "; ekinMin(MeV)=" << ekinMin_;
  if (!ekinMinRegion_.empty()) {
    ptrRegion_.resize(ekinMinRegion_.size(), nullptr);
  }
}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
  g4Track_ = aTrack;
  currentTrack_ = new TrackWithHistory(aTrack, aTrack->GetParentID());

  BeginOfTrack bt(aTrack);
  m_beginOfTrackSignal(&bt);

  trkInfo_ = static_cast<TrackInformation*>(aTrack->GetUserInformation());

  // Always save primaries
  if (trkInfo_->isPrimary()) {
    trackManager_->cleanTracksWithHistory();
    currentTrack_->setToBeSaved();
  }

  if (nullptr != steppingVerbose_) {
    steppingVerbose_->trackStarted(aTrack, false);
    if (aTrack->GetTrackID() == endPrintTrackID_) {
      steppingVerbose_->stopEventPrint();
    }
  }
  double ekin = aTrack->GetKineticEnergy();

#ifdef EDM_ML_DEBUG
  edm::LogVerbatim("DoFineCalo") << "PreUserTrackingAction: Start processing track " << aTrack->GetTrackID()
                                 << " pdgid=" << aTrack->GetDefinition()->GetPDGEncoding()
                                 << " ekin[GeV]=" << ekin / CLHEP::GeV << " vertex[cm]=("
                                 << aTrack->GetVertexPosition().x() / CLHEP::cm << ","
                                 << aTrack->GetVertexPosition().y() / CLHEP::cm << ","
                                 << aTrack->GetVertexPosition().z() / CLHEP::cm << ")"
                                 << " parentid=" << aTrack->GetParentID();
#endif
  if (ekin > ekinMin_) {
    // Each track with energy above the threshold should be saved
    trkInfo_->putInHistory();
  }
}
#include "G4RegionStore.hh"
#include "G4Proton.hh"
void TrackingAction::PostUserTrackingAction(const G4Track* aTrack) {

  auto ff_myscoring = [&](const G4Track* aTrack){
    trackIDmap[aTrack->GetTrackID()] = {aTrack->GetParticleDefinition(), aTrack->GetVertexKineticEnergy()};

    // if no creator process, return early
    if(0 == aTrack->GetParentID() ) return;
    const G4VProcess * track_creator_process = aTrack->GetCreatorProcess();
    if (nullptr == track_creator_process) return;
    G4RegionStore * regionStore = G4RegionStore::GetInstance();
    auto * fRegionEcal = regionStore->FindOrCreateRegion("EcalRegion");
    auto * fRegionHcal = regionStore->FindOrCreateRegion("HcalRegion");
    auto * trackRegion = aTrack->GetLogicalVolumeAtVertex()->GetRegion();
    if(trackRegion != fRegionEcal && trackRegion != fRegionHcal )
        return;

    int hIDe0 = this->GetParticleHistoID(aTrack);
    int hIDe0_n = hIDe0 + 1;
    int hIDe0_pi= hIDe0 + 2;
    int hIDef = hIDe0+3;
    int hIDtf = hIDe0+4;
    double e0 = aTrack->GetVertexKineticEnergy();
    double ef = aTrack->GetKineticEnergy();
    double tf = aTrack->GetLocalTime();

    int pindex = 0;
    auto procIt = fProcNameId.find(track_creator_process->GetProcessName());
    if(fProcNameId.end() == procIt ){
        //std::cerr << "\tAlvaro warning: creator process name <"
        //          << track_creator_process->GetProcessName()
        //          << "> not found in fProcNameId" << std::endl;
        pindex = 0;
    }
    else
        pindex = procIt->second + 1;

    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH2(hIDe0, std::log10(e0) ,pindex);
    if(auto it = trackIDmap.find(aTrack->GetParentID()); it != trackIDmap.end()){
        if(G4Neutron::Neutron() == it->second.first)
        {
            analysisManager->FillH2(hIDe0_n, std::log10(e0) ,pindex);
        }
        else if(G4PionMinus::PionMinus() == it->second.first ||
                G4PionPlus::PionPlus() == it->second.first ||
                G4PionZero::PionZero() == it->second.first
                )
        {
            analysisManager->FillH2(hIDe0_pi, std::log10(e0) ,pindex);
        }
    }
    analysisManager->FillH2(hIDef, std::log10(ef) ,pindex);
    analysisManager->FillH2(hIDtf, std::log10(tf) ,pindex);
  };
  ff_myscoring(aTrack);


  // Tracks in history may be upgraded to stored secondary tracks,
  // which cross the boundary between Tracker and Calo
  int id = aTrack->GetTrackID();
  bool ok = (trkInfo_->storeTrack() || currentTrack_->saved());
  if (trkInfo_->crossedBoundary()) {
    currentTrack_->setCrossedBoundaryPosMom(id, trkInfo_->getPositionAtBoundary(), trkInfo_->getMomentumAtBoundary());
    ok = (ok || saveCaloBoundaryInformation_ || doFineCalo_);
  }
  if (ok) {
    currentTrack_->setToBeSaved();
  }

  bool withAncestor = (trkInfo_->getIDonCaloSurface() == id || trkInfo_->isAncestor());
  bool isInHistory = trkInfo_->isInHistory();

  trackManager_->addTrack(currentTrack_, aTrack, isInHistory, withAncestor);

#ifdef EDM_ML_DEBUG
  edm::LogVerbatim("TrackingAction") << "TrackingAction end track=" << id << "  "
                                     << aTrack->GetDefinition()->GetParticleName() << " proposed to be saved= " << ok
                                     << " end point " << aTrack->GetPosition();
#endif

  if (!isInHistory) {
    delete currentTrack_;
  }

  EndOfTrack et(aTrack);
  m_endOfTrackSignal(&et);
}
