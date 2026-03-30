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
      std::vector<std::pair<std::string,int>> histovector;
      int id = 0;
      for (const auto& [pdg, name] : particles) {
            id = analysisManager->CreateH2(
                "hE0_" + name, "", 1000, -5, 5, nmodels+1, -1, nmodels);
            histovector.push_back({"hE0_" + name, id});
            id = analysisManager->CreateH2(
                "hE0_" + name + "_n", "", 1000, -5, 5, nmodels+1, -1, nmodels);

            histovector.push_back({"hE0_" + name + "_n", id});
            id = analysisManager->CreateH2(
                "hE0_" + name + "_pi", "", 1000, -5, 5, nmodels+1, -1, nmodels);

            histovector.push_back({"hE0_" + name + "_pi", id});

            id = analysisManager->CreateH2(
                "hEf_" + name, "",  1000, -10, 10, nmodels+1, -1, nmodels);

            histovector.push_back({"hEf_" + name, id});

            id = analysisManager->CreateH2(
                "hTf_" + name, "",  1000, -10, 10, nmodels+1, -1, nmodels);

            histovector.push_back({"hTf_" + name, id});
      }
      for(auto [hname, hid] : histovector)
          std::cout << hname.c_str() << ' ' << hid << std::endl;
  } // end create histograms
  int nH2s = analysisManager->GetNofH2s();
  for(int i = 0; i<nH2s; ++i)
    std::cout << "\t" << i << "\t" << analysisManager->GetH2Name(i) << std::endl;

  this->PrintGeant4Configuration();
}

void RunAction::EndOfRunAction(const G4Run* aRun) {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  EndOfRun r(aRun);
  m_endOfRunSignal(&r);
}


#include "G4ProductionCutsTable.hh"
#include "G4RegionStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4UserLimits.hh"
#include "G4ios.hh"
#include "G4HadronicParameters.hh"
#include "G4EmParameters.hh"
#include "G4HadronicProcessStore.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"
#include "G4PhysicsModelCatalog.hh"

#include <fstream>

void RunAction::PrintGeant4Configuration()
{
    static std::ofstream dumpFile("geant4_dump.txt");
    //save buffer to restore it later
    auto cout_buf = G4cout.rdbuf();
    auto cerr_buf = G4cerr.rdbuf();
    // redirect G4cout
    G4cout.rdbuf(dumpFile.rdbuf());
    G4cerr.rdbuf(dumpFile.rdbuf());

    G4cout << "==============================" << G4endl;
    G4cout << " GEANT4 CONFIGURATION DUMP " << G4endl;
    G4cout << "==============================" << G4endl;


    G4cout << "\n=== Production Cuts Table ===\n" << G4endl;

    auto pct = G4ProductionCutsTable::GetProductionCutsTable();
    pct->DumpCouples();


    G4cout << "================================\n";
    G4cout << "\n=== Regions and Production Cuts ===\n" << G4endl;

    auto regionStore = G4RegionStore::GetInstance();

    for (size_t i = 0; i < regionStore->size(); ++i)
    {
        auto region = (*regionStore)[i];
        G4cout << "Region: " << region->GetName() << G4endl;

        auto cuts = region->GetProductionCuts();
        if (cuts)
        {
            G4cout << "  gamma cut: " << cuts->GetProductionCut("gamma") << G4endl;
            G4cout << "  e- cut   : " << cuts->GetProductionCut("e-") << G4endl;
            G4cout << "  e+ cut   : " << cuts->GetProductionCut("e+") << G4endl;
            G4cout << "  proton cut: " << cuts->GetProductionCut("proton") << G4endl;
        }
        else
        {
            G4cout << "  No production cuts attached" << G4endl;
        }
        auto * ul = region->GetUserLimits();
        if (ul)
        {
            G4Track dummyTrack;
            G4cout << "  MaxAllowedStep: "
                << ul->GetMaxAllowedStep(dummyTrack) << G4endl;

            G4cout << "  MaxTrackLength: "
                << ul->GetUserMaxTrackLength(dummyTrack) << G4endl;

            G4cout << "  MaxTime: "
                << ul->GetUserMaxTime(dummyTrack) << G4endl;

            G4cout << "  MinEkine: "
                << ul->GetUserMinEkine(dummyTrack) << G4endl;

            G4cout << "  MinRange: "
                << ul->GetUserMinRange(dummyTrack) << G4endl;
        }
        else
        {
            G4cout << "  No user limits attached" << G4endl;
        }
    }
    {

        G4cout << "================================\n";
        G4cout << "\n=== User Limits per Logical Volume ===\n" << G4endl;

        auto lvStore = G4LogicalVolumeStore::GetInstance();

        G4Track dummyTrack;
        for (auto lv : *lvStore)
        {
            auto ul = lv->GetUserLimits();
            if (ul)
            {
                G4cout << "LogicalVolume: " << lv->GetName() << G4endl;

                G4cout << "  MaxAllowedStep: "
                    << ul->GetMaxAllowedStep(dummyTrack) << G4endl;

                G4cout << "  MaxTrackLength: "
                    << ul->GetUserMaxTrackLength(dummyTrack) << G4endl;

                G4cout << "  MaxTime: "
                    << ul->GetUserMaxTime(dummyTrack) << G4endl;

                G4cout << "  MinEkine: "
                    << ul->GetUserMinEkine(dummyTrack) << G4endl;

                G4cout << "  MinRange: "
                    << ul->GetUserMinRange(dummyTrack) << G4endl;
            }
        }
    }

    std::set<std::string> processnames;
    G4cout << "================================\n";
    G4cout << "\n=== Processes per Particle ===\n" << G4endl;

    auto particleTable = G4ParticleTable::GetParticleTable();
    auto it = particleTable->GetIterator();
    it->reset();

    while ((*it)())
    {
        auto particle = it->value();
        auto pm = particle->GetProcessManager();

        if (!pm) continue;

        G4cout << "\nParticle (pdgID): "
            << particle->GetParticleName()
            << " (" << particle->GetPDGEncoding() << ")"
            << G4endl;

        G4ProcessVector * pv = pm->GetProcessList();

        for (int i = 0; i < pm->GetProcessListLength(); ++i)
        {
            auto proc = (*pv)[i];

            G4cout << "  Process: "
                << proc->GetProcessName()
                << G4endl;

            processnames.emplace(proc->GetProcessName());

            // -------------------------------------------------
            // only for hadronic processes
            // -------------------------------------------------
            auto hadProc = dynamic_cast<G4HadronicProcess*>((*pv)[i]);

            if (!hadProc) continue;

            auto& models = hadProc->GetHadronicInteractionList();

            for (auto model : models)
            {
                G4cout << "    Model: "
                    << model->GetModelName()
                    << " | Emin = " << model->GetMinEnergy()/CLHEP::GeV << " GeV"
                    << " | Emax = " << model->GetMaxEnergy()/CLHEP::GeV << " GeV"
                    << G4endl;
            }

        }
    }


    G4cout << "===== G4EmParameters =====\n";
    G4EmParameters::Instance()->StreamInfo(G4cout);
    G4EmParameters::Instance()->Dump();

    G4cout << "================================\n";
    G4cout << "===== G4HadronicParameters =====\n";
    G4HadronicParameters* hadronic_params = G4HadronicParameters::Instance();
    // Energy limits
    G4cout << "MaxEnergy: " << hadronic_params->GetMaxEnergy() << "\n";

    G4cout << "MinEnergyTransitionFTF_Cascade: "
              << hadronic_params->GetMinEnergyTransitionFTF_Cascade() << "\n";
    G4cout << "MaxEnergyTransitionFTF_Cascade: "
              << hadronic_params->GetMaxEnergyTransitionFTF_Cascade() << "\n";

    G4cout << "MinEnergyTransitionQGS_FTF: "
              << hadronic_params->GetMinEnergyTransitionQGS_FTF() << "\n";
    G4cout << "MaxEnergyTransitionQGS_FTF: "
              << hadronic_params->GetMaxEnergyTransitionQGS_FTF() << "\n";

    G4cout << "MinEnergyINCLXX_Pbar: "
              << hadronic_params->GetMinEnergyINCLXX_Pbar() << "\n";
    G4cout << "MaxEnergyINCLXX_Pbar: "
              << hadronic_params->GetMaxEnergyINCLXX_Pbar() << "\n";

    G4cout << "EnergyThresholdForHeavyHadrons: "
              << hadronic_params->EnergyThresholdForHeavyHadrons() << "\n";

    // Cross section factors
    G4cout << "XSFactorNucleonInelastic: "
              << hadronic_params->XSFactorNucleonInelastic() << "\n";
    G4cout << "XSFactorNucleonElastic: "
              << hadronic_params->XSFactorNucleonElastic() << "\n";

    G4cout << "XSFactorPionInelastic: "
              << hadronic_params->XSFactorPionInelastic() << "\n";
    G4cout << "XSFactorPionElastic: "
              << hadronic_params->XSFactorPionElastic() << "\n";

    G4cout << "XSFactorHadronInelastic: "
              << hadronic_params->XSFactorHadronInelastic() << "\n";
    G4cout << "XSFactorHadronElastic: "
              << hadronic_params->XSFactorHadronElastic() << "\n";

    G4cout << "XSFactorEM: "
              << hadronic_params->XSFactorEM() << "\n";

    // Flags
    G4cout << "EnableBCParticles: "
              << hadronic_params->EnableBCParticles() << "\n";
    G4cout << "EnableHyperNuclei: "
              << hadronic_params->EnableHyperNuclei() << "\n";
    G4cout << "ApplyFactorXS: "
              << hadronic_params->ApplyFactorXS() << "\n";

    G4cout << "EnableCRCoalescence: "
              << hadronic_params->EnableCRCoalescence() << "\n";

    G4cout << "EnableIntegralInelasticXS: "
              << hadronic_params->EnableIntegralInelasticXS() << "\n";
    G4cout << "EnableIntegralElasticXS: "
              << hadronic_params->EnableIntegralElasticXS() << "\n";

    G4cout << "EnableDiffDissociationForBGreater10: "
              << hadronic_params->EnableDiffDissociationForBGreater10() << "\n";

    G4cout << "EnableCoherentChargeExchange: "
              << hadronic_params->EnableCoherentChargeExchange() << "\n";

    G4cout << "EnableNeutronGeneralProcess: "
              << hadronic_params->EnableNeutronGeneralProcess() << "\n";

    // Verbosity
    G4cout << "VerboseLevel: "
              << hadronic_params->GetVerboseLevel() << "\n";

    // Energy-momentum conservation parameters
    G4cout << "EPRelativeLevel: "
              << hadronic_params->GetEPRelativeLevel() << "\n";
    G4cout << "EPAbsoluteLevel: "
              << hadronic_params->GetEPAbsoluteLevel() << "\n";
    G4cout << "EPReportLevel: "
              << hadronic_params->GetEPReportLevel() << "\n";

    G4cout << "BinaryDebug: "
              << hadronic_params->GetBinaryDebug() << "\n";

    // Environment variables
    G4cout << "DirPARTICLEXS: "
              << hadronic_params->GetDirPARTICLEXS() << "\n";
    G4cout << "PhysListDocDir: "
              << hadronic_params->GetPhysListDocDir() << "\n";
    G4cout << "PhysListName: "
              << hadronic_params->GetPhysListName() << "\n";

    // Thresholds
    G4cout << "NeutronKineticEnergyThresholdForSVT: "
              << hadronic_params->GetNeutronKineticEnergyThresholdForSVT() << "\n";

    G4cout << "TimeThresholdForRadioactiveDecay: "
              << hadronic_params->GetTimeThresholdForRadioactiveDecay() << "\n";

    G4cout << "================================\n";
    G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
    deex->StreamInfo(G4cout);

    G4cout << "================================\n";
    G4cout << "===== G4HadronicProcessStore =====\n";
    G4HadronicProcessStore::Instance()->Dump(1);

    G4cout << "================================\n";
    G4PhysicsModelCatalog::PrintAllInformation();


    dumpFile.flush();
    G4cout.rdbuf(cout_buf);
    G4cerr.rdbuf(cerr_buf);
}
