#include "SimG4Core/Application/interface/RunManagerMT.h"
#include "SimG4Core/Application/interface/PrimaryTransformer.h"
#include "SimG4Core/Application/interface/SimRunInterface.h"
#include "SimG4Core/Application/interface/RunAction.h"
#include "SimG4Core/Application/interface/ParametrisedEMPhysics.h"
#include "SimG4Core/Application/interface/ExceptionHandler.h"

#include "SimG4Core/Geometry/interface/DDDWorld.h"
#include "SimG4Core/Geometry/interface/CustomUIsession.h"
#include "SimG4Core/Geometry/interface/CMSG4RegionReporter.h"

#include "SimG4Core/Physics/interface/PhysicsListFactory.h"
#include "SimG4Core/PhysicsLists/interface/CMSMonopolePhysics.h"
#include "SimG4Core/CustomPhysics/interface/CMSExoticaPhysics.h"

#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"

#include "SimG4Core/Notification/interface/SimTrackManager.h"
#include "SimG4Core/Notification/interface/BeginOfJob.h"
#include "SimG4Core/Notification/interface/CurrentG4Track.h"
#include "SimG4Core/Geometry/interface/CMSG4CheckOverlap.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/DDCMS/interface/DDCompactView.h"

#include "SimDataFormats/Forward/interface/LHCTransportLinkContainer.h"

#include "HepPDT/ParticleDataTable.hh"

#include "G4Timer.hh"
#include "G4GeometryManager.hh"
#include "G4ScoringManager.hh"
#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "G4MTRunManagerKernel.hh"
#include "G4UImanager.hh"

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleTable.hh"
#include "G4CascadeInterface.hh"
#include "G4EmParameters.hh"
#include "G4LossTableManager.hh"
#include "G4HadronicParameters.hh"
#include "G4NuclearLevelData.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysListUtil.hh"
#include "G4PhysicsListHelper.hh"

#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>

RunManagerMT::RunManagerMT(edm::ParameterSet const& p)
    : m_PhysicsTablesDir(p.getUntrackedParameter<std::string>("PhysicsTablesDirectory", "")),
      m_StorePhysicsTables(p.getUntrackedParameter<bool>("StorePhysicsTables", false)),
      m_RestorePhysicsTables(p.getUntrackedParameter<bool>("RestorePhysicsTables", false)),
      m_check(p.getUntrackedParameter<bool>("CheckGeometry")),
      m_geoFromDD4hep(p.getParameter<bool>("g4GeometryDD4hepSource")),
      m_score(p.getParameter<bool>("UseCommandBaseScorer")),
      m_stepverb(p.getUntrackedParameter<int>("SteppingVerbosity", 0)),
      m_regionFile(p.getUntrackedParameter<std::string>("FileNameRegions", "")),
      m_pPhysics(p.getParameter<edm::ParameterSet>("Physics")),
      m_pRunAction(p.getParameter<edm::ParameterSet>("RunAction")),
      m_Init(p.getParameter<edm::ParameterSet>("Init")),
      m_G4Commands(p.getParameter<std::vector<std::string> >("G4Commands")),
      m_EnableHGCalSubregions(p.getUntrackedParameter<bool>("EnableHGCalSubregions", false)){
  m_physicsList.reset(nullptr);
  m_world.reset(nullptr);
  m_runInterface.reset(nullptr);

  m_kernel = new G4MTRunManagerKernel();
  m_stateManager = G4StateManager::GetStateManager();
  double th = p.getParameter<double>("ThresholdForGeometryExceptions") * CLHEP::GeV;
  bool tr = p.getParameter<bool>("TraceExceptions");
  m_stateManager->SetExceptionHandler(new ExceptionHandler(th, tr));
  if (m_check) {
    m_CheckOverlap = p.getUntrackedParameter<edm::ParameterSet>("G4CheckOverlap");
  }
  m_UIsession = new CustomUIsession();
  G4UImanager::GetUIpointer()->SetCoutDestination(m_UIsession);
  G4UImanager::GetUIpointer()->SetMasterUIManager(true);
  G4PhysListUtil::InitialiseParameters();
  G4LossTableManager::Instance();
}

RunManagerMT::~RunManagerMT() { delete m_UIsession; }

void RunManagerMT::initG4(const DDCompactView* pDD,
                          const cms::DDCompactView* pDD4hep,
                          const HepPDT::ParticleDataTable* fPDGTable) {
  if (m_managerInitialized) {
    edm::LogWarning("SimG4CoreApplication") << "RunManagerMT::initG4 was already done - exit";
    return;
  }
  bool cuts = m_pPhysics.getParameter<bool>("CutsPerRegion");
  bool protonCut = m_pPhysics.getParameter<bool>("CutsOnProton");
  int verb = m_pPhysics.getUntrackedParameter<int>("Verbosity", 0);
  edm::LogVerbatim("SimG4CoreApplication")
      << "RunManagerMT: start initialising of geometry DD4hep: " << m_geoFromDD4hep << "\n"
      << "              cutsPerRegion: " << cuts << " cutForProton: " << protonCut << "\n"
      << "              G4 verbosity: " << verb;

  G4Timer timer;
  timer.Start();

  m_world = std::make_unique<DDDWorld>(pDD, pDD4hep, m_catalog, verb, cuts, protonCut);
  G4VPhysicalVolume* world = m_world.get()->GetWorldVolume();

  m_kernel->SetVerboseLevel(verb);
  edm::LogVerbatim("SimG4CoreApplication")
      << "RunManagerMT: Define cuts: " << cuts << " Geant4 run manager verbosity: " << verb;

  const G4RegionStore* regStore = G4RegionStore::GetInstance();
  const G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
  const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
  unsigned int numPV = pvs->size();
  unsigned int numLV = lvs->size();
  unsigned int nn = regStore->size();
  edm::LogVerbatim("SimG4CoreApplication")
      << "RunManagerMT: " << numPV << " physical volumes; " << numLV << " logical volumes; " << nn << " regions.";

#if G4VERSION_NUMBER >= 1130
  G4GeometryManager::GetInstance()->RequestParallelOptimisation(false, false);
#endif

  if (m_check) {
    m_kernel->SetVerboseLevel(2);
  }
  m_kernel->DefineWorldVolume(world, true);
  m_registry.dddWorldSignal_(m_world.get());
  m_stateManager->SetNewState(G4State_PreInit);

  if (m_EnableHGCalSubregions) {
    addRegions();
  }
  // Create physics list
  edm::LogVerbatim("SimG4CoreApplication") << "RunManagerMT: create PhysicsList";

  std::unique_ptr<PhysicsListMakerBase> physicsMaker(
      PhysicsListFactory::get()->create(m_pPhysics.getParameter<std::string>("type")));
  if (physicsMaker.get() == nullptr) {
    throw cms::Exception("Configuration") << "Unable to find the Physics list requested";
  }
  m_physicsList = physicsMaker->make(m_pPhysics, m_registry);

  PhysicsList* phys = m_physicsList.get();
  if (phys == nullptr) {
    throw cms::Exception("Configuration") << "Physics list construction failed!";
  }
  if (m_stepverb > 0) {
    verb = std::max(verb, 1);
  }
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
  G4EmParameters::Instance()->SetVerbose(verb);
  G4EmParameters::Instance()->SetWorkerVerbose(std::max(verb - 1, 0));
  G4PhysicsListHelper::GetPhysicsListHelper();

  // exotic particle physics
  double monopoleMass = m_pPhysics.getUntrackedParameter<double>("MonopoleMass", 0);
  if (monopoleMass > 0.0) {
    phys->RegisterPhysics(new CMSMonopolePhysics(fPDGTable, m_pPhysics));
  }
  bool exotica = m_pPhysics.getUntrackedParameter<bool>("ExoticaTransport", false);
  if (exotica) {
    CMSExoticaPhysics exo(phys, m_pPhysics);
  }

  // adding GFlash, Russian Roulette for eletrons and gamma,
  // step limiters on top of any Physics Lists
  phys->RegisterPhysics(new ParametrisedEMPhysics("EMoptions", m_pPhysics));

  if (m_RestorePhysicsTables) {
    m_physicsList->SetPhysicsTableRetrieved(m_PhysicsTablesDir);
  }
  edm::LogVerbatim("SimG4CoreApplication") << "RunManagerMT: start initialisation of PhysicsList for master";

  m_physicsList->SetDefaultCutValue(m_pPhysics.getParameter<double>("DefaultCutValue") * CLHEP::cm);
  m_physicsList->SetCutsWithDefault();
  m_kernel->SetPhysics(phys);

  edm::LogVerbatim("SimG4CoreApplication") << "RunManagerMT: PhysicsList and cuts are defined";

  // Enable couple transportation
  if (m_score) {
    G4ScoringManager* scManager = G4ScoringManager::GetScoringManager();
    scManager->SetVerboseLevel(1);
  }
  // Geant4 UI commands before initialisation of physics
  if (!m_G4Commands.empty()) {
    edm::LogVerbatim("SimG4CoreApplication") << "RunManagerMT: Requested UI commands: ";
    for (const std::string& command : m_G4Commands) {
      edm::LogVerbatim("SimG4CoreApplication") << "    " << command;
      G4UImanager::GetUIpointer()->ApplyCommand(command);
    }
  }

  setupVoxels();

  m_stateManager->SetNewState(G4State_Init);
  edm::LogVerbatim("SimG4CoreApplication") << "RunManagerMT: G4State is Init";
  m_kernel->InitializePhysics();
  if (verb > 0) {
    G4EmParameters::Instance()->Dump();
  }
  m_kernel->SetUpDecayChannels();

  if (m_kernel->RunInitialization()) {
    m_managerInitialized = true;
  } else {
    throw cms::Exception("LogicError") << "G4RunManagerKernel initialization failed!";
  }

  if (m_check) {
    checkVoxels();
  }

  if (m_StorePhysicsTables) {
    std::ostringstream dir;
    dir << m_PhysicsTablesDir << '\0';
    std::string cmd = std::string("/control/shell mkdir -p ") + m_PhysicsTablesDir;
    if (!std::ifstream(dir.str().c_str(), std::ios::in))
      G4UImanager::GetUIpointer()->ApplyCommand(cmd);
    m_physicsList->StorePhysicsTable(m_PhysicsTablesDir);
  }
  // Appload nuclear level data up to Z=84
  G4NuclearLevelData::GetInstance()->UploadNuclearLevelData(84);

  if (verb > 1) {
    m_physicsList->DumpCutValuesTable();
  }
  edm::LogVerbatim("SimG4CoreApplication")
      << "RunManagerMT: Physics is initilized, now initialise user actions, verb=" << verb;

  initializeUserActions();

  // G4Region dump file name
  runForPhase2();

  // Geometry checks
  if (m_check || !m_regionFile.empty()) {
    CMSG4CheckOverlap check(m_CheckOverlap, m_regionFile, m_UIsession, world);
  }

  m_stateManager->SetNewState(G4State_PreInit);
  G4HadronicParameters::Instance()->SetVerboseLevel(std::max(verb - 1, 0));

  // If the Geant4 particle table is needed, decomment the lines below
  //
  //G4ParticleTable::GetParticleTable()->DumpTable("ALL");
  //
  m_stateManager->SetNewState(G4State_GeomClosed);
  m_currentRun = new G4Run();
  m_userRunAction->BeginOfRunAction(m_currentRun);
  timer.Stop();
  G4cout.precision(4);
  G4cout << "RunManagerMT: initG4 done " << timer << G4endl;
}

void RunManagerMT::initializeUserActions() {
  m_runInterface = std::make_unique<SimRunInterface>(this, true);
  m_userRunAction = new RunAction(m_pRunAction, m_runInterface.get(), true);
  Connect(m_userRunAction);
}

void RunManagerMT::Connect(RunAction* runAction) {
  runAction->m_beginOfRunSignal.connect(m_registry.beginOfRunSignal_);
  runAction->m_endOfRunSignal.connect(m_registry.endOfRunSignal_);
}

void RunManagerMT::stopG4() {
  edm::LogVerbatim("SimG4CoreApplication") << "RunManagerMT::stopG4";
  G4GeometryManager::GetInstance()->OpenGeometry();
  m_stateManager->SetNewState(G4State_Quit);
  if (!m_runTerminated) {
    terminateRun();
  }
  edm::LogVerbatim("SimG4CoreApplication") << "RunManagerMT::stopG4 done";
}

void RunManagerMT::terminateRun() {
  edm::LogVerbatim("SimG4CoreApplication") << "RunManagerMT::terminateRun";
  if (nullptr != m_userRunAction) {
    m_userRunAction->EndOfRunAction(m_currentRun);
    delete m_userRunAction;
    m_userRunAction = nullptr;
  }
  if (!m_runTerminated) {
    m_kernel->RunTermination();
  }
  m_runTerminated = true;
  edm::LogVerbatim("SimG4CoreApplication") << "RunManagerMT::terminateRun done";
}

void RunManagerMT::checkVoxels() {
  const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
  int numLV = (int)lvs->size();
  edm::LogVerbatim("SimG4CoreApplication") << "RunManagerMT: nLV=" << numLV;
  int nvox = 0;
  int nslice = 0;
  for (int i = 0; i < numLV; ++i) {
    auto lv = (*lvs)[i];
    auto nd = lv->GetNoDaughters();
    auto vox = lv->GetVoxelHeader();
    auto sma = lv->GetSmartless();
    auto reg = lv->GetRegion();
    std::size_t nsl = (nullptr == vox) ? 0 : vox->GetNoSlices();
    if (0 < nsl) {
      nslice += nsl;
      std::string rname = (nullptr != reg) ? reg->GetName() : "";
      edm::LogVerbatim("Voxels") << " " << i << ". Nd=" << nd << " Nsl=" << nsl << " Smartless=" << sma << " "
                                 << lv->GetName() << " Region: " << rname;
      ++nvox;
    }
  }
  edm::LogVerbatim("SimG4CoreApplication")
      << "RunManagerMT: nLV=" << numLV << " NlvVox=" << nvox << " Nslices=" << nslice;
}

void RunManagerMT::setupVoxels() {
  double density = m_Init.getParameter<double>("DefaultVoxelDensity");
  std::vector<std::string> rnames = m_Init.getParameter<std::vector<std::string> >("VoxelRegions");
  std::vector<double> rdensities = m_Init.getParameter<std::vector<double> >("VoxelDensityPerRegion");
  int nr = 0;
  std::size_t n = rnames.size();
  if (n == rdensities.size()) {
    nr = (int)n;
  }
  const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
  for (auto const& lv : *lvs) {
    double den = density;
    if (0 < nr) {
      std::string nam = lv->GetRegion()->GetName();
      for (int i = 0; i < nr; ++i) {
        if (nam == rnames[i]) {
          den = rdensities[i];
          break;
        }
      }
    }
    lv->SetSmartless(den);
  }
  edm::LogVerbatim("SimG4CoreApplication")
      << "RunManagerMT: default voxel density=" << density << "; number of regions with special density " << nr;
}

void RunManagerMT::runForPhase2() {
  const G4RegionStore* regStore = G4RegionStore::GetInstance();
  for (auto const& r : *regStore) {
    const G4String& name = r->GetName();
    if (name == "HGCalRegion" || name == "FastTimerRegionETL" || name == "FastTimerRegionBTL") {
      m_isPhase2 = true;
      break;
    }
  }
}

// This method should be extended in order to add new regions via names of logical volume
// and to add a new production cuts for these regions
void RunManagerMT::addRegions() {
  struct mat_cut_couple_t{
    // Constructor with 1 cut delegates to ctor with cut per particle
    mat_cut_couple_t(std::string imatname, std::string ilvname_pattern, double icut_mm)
        : mat_cut_couple_t(std::move(imatname), std::move(ilvname_pattern),
                           icut_mm, icut_mm, icut_mm) {}

    // Constructor with cuts per particle
    mat_cut_couple_t(std::string imatname, std::string ilvname_pattern,
                     double icutg_mm, double icute_mm, double icutp_mm)
        : matname(std::move(imatname)),
          lvname_pattern(std::move(ilvname_pattern)),
          cutg_mm(icutg_mm), cute_mm(icute_mm), cutp_mm(icutp_mm) {}
    std::string matname;
    std::string lvname_pattern;
    double cutg_mm;
    double cute_mm;
    double cutp_mm;
  };
  std::vector<mat_cut_couple_t> mat_cut_couple_vector = {
    // list of material, logical volume regex and cut (in mm)
    // global cut, if exact regex, material is ignored
    {"global","HGCal",3},
    // common materials to EE and HE
    {"Epoxy","^HGCal.*",0.1},
    {"HGC_G10-FR4","^HGCal.*",0.1},
    {"Silicon","^HGCal.*",0.03},
    //EE cuts
    {"Kapton","^HGCalEEKapton.*",0.1},
    {"WCu","^HGCal.*",0.3},
    {"Lead","^HGCal.*",2},
    {"StainlessSteel","^HGCalEEAbsorberStainlessSteel.*",0.1},
    {"HGC_EEConnector","^HGCal.*",0.1},
    {"Copper","^HGCalEEAbsorberCopper.*",0.1},
    {"Copper","^HGCalEECoolingPlate.*",3},
    //HE cuts
    {"Polystyrene","^HGCalScintillator.*",0.03},
    {"H_Scintillator","^HGCalHEScintillatorSensitive.*",0.03},
    {"Kapton","^HGCalHEKapton.*",0.03},
    {"Copper","^HGCalHECoolingPlate.*",3},
    {"HGC_HEConnector","^HGCal.*",0.1},
    {"StainlessSteel","^HGCalHESteelCover.*",2},
    {"StainlessSteel","^HGCalHEAbsorber.*",41.5},
    {"StainlessSteel","^HGCalBackPlate.*",100},
  };

  /// this lambda function creates a region a given couple material-regex of lvname
  auto createG4Region = [&](const mat_cut_couple_t & couple)->void {
    // set new cuts for existing HGCal envelope region
    if("global" == couple.matname)
    {
      // retrieve existing region
      // if it does not exist, G4 throws JustWarning G4Exception
	    G4Region * rg = G4RegionStore::GetInstance()->GetRegion("HGCalRegion");
        G4ProductionCuts * HGCalmatcuts = rg->GetProductionCuts();
        // Set cut values (in mm)
        HGCalmatcuts->SetProductionCut(couple.cutg_mm * CLHEP::mm, G4ProductionCuts::GetIndex("gamma"));
        HGCalmatcuts->SetProductionCut(couple.cute_mm * CLHEP::mm, G4ProductionCuts::GetIndex("e-"));
        HGCalmatcuts->SetProductionCut(couple.cute_mm * CLHEP::mm, G4ProductionCuts::GetIndex("e+"));
        HGCalmatcuts->SetProductionCut(couple.cutp_mm * CLHEP::mm, G4ProductionCuts::GetIndex("proton"));
	    return;
    }
    // existing method addG4Region requires vector with ptr to logical volumes and cuts
    // the mat_cut_couple_t object contains material name and a regex of lv name
    // 1. Retrieve pointer to material
    const G4Material * mat_ptr = G4Material::GetMaterial(couple.matname);
    if(! mat_ptr)
    {
      std::string error_message = std::string("Input material <")
                                + couple.matname
                                + "> does not exist\nIn function: "
                                + __PRETTY_FUNCTION__;
      throw std::runtime_error(error_message);
    }
    // 2. Identify lv by regex of its name
    std::regex lvname_pattern_rg(couple.lvname_pattern);

    // 2.1 Accumulate ptr of lv in a vector
    std::vector<G4LogicalVolume*> lv_vector;

    // 2.2 Loop over the lv store and append to vector lv if material & name regex match
    G4LogicalVolumeStore * lv_store = G4LogicalVolumeStore::GetInstance();
    for (const auto& lv : *lv_store)
    {
        if( (lv->GetMaterial() == mat_ptr) && std::regex_match( lv->GetName(), lvname_pattern_rg) )
            lv_vector.push_back(lv);
    }

    // 3. Create unique region name
    std::string rg_name = "HGCalRegion_" + couple.matname + "_" + couple.lvname_pattern;

    // 4. Pass the vector with lv, region name and cuts to the method to create the region
    this->addG4Region(lv_vector,
                      rg_name,
                      couple.cutg_mm * CLHEP::mm,
                      couple.cute_mm * CLHEP::mm,
                      couple.cute_mm * CLHEP::mm,
                      couple.cutp_mm * CLHEP::mm
                      );
  };
  //---------------- end of createG4Region

  // iterate over all the material-lv regex name cut couples
  for( auto & mc : mat_cut_couple_vector)
    createG4Region(mc);

  // Dump each lv and its assigned region and cuts
  CMSG4RegionReporter rep;
  rep.ReportRegions("g4region.txt");
}

// This is a utility method to add a G4Region
void RunManagerMT::addG4Region(const std::vector<G4LogicalVolume*>& v,
                               const std::string& rName,
                               double cutg,
                               double cute,
                               double cutp,
                               double cuti) {
  if (v.empty()) {
    return;
  }
  auto reg = new G4Region((G4String)rName);
  for (auto const& lv : v) {
    reg->AddRootLogicalVolume(lv, true);
  }
  auto cuts = new G4ProductionCuts();
  cuts->SetProductionCut(cutg, 0);
  cuts->SetProductionCut(cute, 1);
  cuts->SetProductionCut(cutp, 2);
  cuts->SetProductionCut(cuti, 3);
  reg->SetProductionCuts(cuts);
}
