#ifndef SECONDARYESCAPECOUNTER_HH
#define SECONDARYESCAPECOUNTER_HH

#include "G4Track.hh"

#include "TH1D.h"
#include "TFile.h"
#include "G4ParticleTable.hh"

class SecondaryEscapeCounter {
public:
    SecondaryEscapeCounter(const G4String & particleName = "")
        : fParticleNameFilter(particleName),
          fTotalSecondaries(0), fEscapingSecondaries(0) {
              std::string fHistTotalName;
              std::string fHistEscapingName;
              if(fParticleNameFilter.empty())
              {
                  fHistTotalName = "hTotalSecondaries";
                  fHistEscapingName = "hTotalSecondariesEscaping";
              }
              else {
                  std::string s = particleName;
                  // remove non alpha characters
                  s.erase(std::remove_if(s.begin(), s.end(),
                           [](unsigned char c){ return !std::isalpha(c); }),
                    s.end());
                  fHistTotalName = "hTotal" + s;
                  fHistEscapingName = "hTotalEscaping" + s;

                  fParticleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
              }
              fHistTotal = std::make_unique<TH1D>(fHistTotalName.c_str(),fHistTotalName.c_str(), 30000, 0, 3e6);
              fHistEscaping = std::make_unique<TH1D>(fHistEscapingName.c_str(),fHistEscapingName.c_str(), 30000, 0, 3e6);

              fHistTotal->SetDirectory(0);
              fHistEscaping->SetDirectory(0);



        }

    void RegisterCreation(const G4Track* track) {
        // return if particle is not a secondary
        if (track->GetParentID() <= 0)
            return;

        // return if no nullptr (all secondaries) or if particle definition does not match
        if (fParticleDefinition && fParticleDefinition != track->GetDefinition())
            return;

        // register creation of secondary
        fTotalSecondaries++;
        fCreationVolume[track->GetTrackID()] = track->GetVolume();
    }


    void RegisterEnd(const G4Track* track)
    {
        // return if particle is not a secondary
        if (track->GetParentID() <= 0)
            return;

        // return if no nullptr (all secondaries) or if particle definition does not match
        if (fParticleDefinition && fParticleDefinition != track->GetDefinition())
            return;

        // return if secondary was not added to the map (this should not happen)
        auto it = fCreationVolume.find(track->GetTrackID());
        if (it == fCreationVolume.end())
        {
            std::cout << "\twarning in secondary counter, RegisterEnd, secondary not registered at begining!!" << std::endl;
            return;
        }

        const G4VPhysicalVolume* creationVol = it->second;

        const G4Step* step = track->GetStep();
        const G4VPhysicalVolume* endVol =
            step ? step->GetPostStepPoint()->GetPhysicalVolume() : nullptr;

        if (endVol != creationVol)
            fEscapingSecondaries++;

        fCreationVolume.erase(it);
    }


    void ResetCounters() {
        fTotalSecondaries = 0;
        fEscapingSecondaries = 0;
        fCreationVolume.clear();
    }
    void FillHistograms() {
        fHistTotal->Fill(fTotalSecondaries);
        fHistEscaping->Fill(fEscapingSecondaries);
    }
    void FillHistogramsAndResetCounters() {
        FillHistograms();
        ResetCounters();
    }
    void WriteHistogram(TFile * f) {
        f->cd();
        fHistTotal->Write();
        fHistEscaping->Write();
    }


    // Getters
    G4int GetTotalSecondaries() const { return fTotalSecondaries; }
    G4int GetEscapingSecondaries() const { return fEscapingSecondaries; }
    const G4String& GetParticleNameFilter() const { return fParticleNameFilter; }
    const G4ParticleDefinition * GetParticleDefinition() const { return fParticleDefinition; }

private:
    G4String fParticleNameFilter; // empty = all
    G4ParticleDefinition * fParticleDefinition = {nullptr};

    G4int fTotalSecondaries;
    G4int fEscapingSecondaries;
    std::unordered_map<G4int, const G4VPhysicalVolume*> fCreationVolume;

    std::unique_ptr<TH1D> fHistTotal;
    std::unique_ptr<TH1D> fHistEscaping;
};

#endif

