##
## my-script-Run4-G4HepEm.py
##
import FWCore.ParameterSet.Config as cms

RUN4=True
G4HEPEM=True

if RUN4:
    print ("===== RUN4 ===== ")
else:
    print ("===== RUN3 ===== ")
if G4HEPEM:
    print ("===== G4HepEm: YES ===== ")
else:
    print ("===== G4HepEm: No   ===== ")



if RUN4:
    ## ------------------------------------------------------------------
    ## RUN4
    ENERGY=14000.0
    from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
    process = cms.Process('gg2ttbar', Phase2C17I13M9)

    process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
    process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
else:
    ## ------------------------------------------------------------------
    ## RUN3
    ENERGY=13600.0
    from Configuration.Eras.Era_Run3_dd4hep_cff import Run3_dd4hep
    process = cms.Process('gg2ttbar', Run3_dd4hep)

    process.load('Configuration.Geometry.GeometryDD4hepExtended2024Reco_cff')
    process.load('IOMC.EventVertexGenerators.VtxSmearedRun3RoundOptics25ns13TeVLowSigmaZ_cfi')
    process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13p6TeVEarly2022Collision_cfi')
    ## ------------------------------------------------------------------


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('SimG4CMS.Calo.CaloSimHitStudy_cfi')

process.Timing = cms.Service('Timing')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000) ### NUMER OF EVENTS TO BE SET
#    input = cms.untracked.int32(50000)
#     input = cms.untracked.int32(100000)
)

process.options.numberOfThreads = 4 ### NUMBER OF THREADS TO BE SET

process.source = cms.Source('EmptySource')

'''
process.generator = cms.EDFilter('Pythia8ConcurrentGeneratorFilter',
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring('gg2ttbar'),
        gg2ttbar = cms.vstring(
            'Tune:pp = -1',
            'Top:gg2ttbar = on',
        ),
    ),
    comEnergy = cms.double(14000.0),
)
process.generator.pythiaHepMCVerbosity = cms.untracked.bool(False)
process.generator.pythiaPylistVerbosity = cms.untracked.int32(0)
'''


##
## ttbar event generation tajken from:
## /Validation/HGCalValidation/scripts/testHGCalSimTTBar_cfg.py
process.generator = cms.EDProducer(
    "FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(11),          # e- 11, pi- -211, mu- 13
        MinEta = cms.double(3.0),
        MaxEta = cms.double(3.0),         # fixed, eta = -ln(tan(theta/2)) (theta 0.165 rad -> tan theta = 1/6)
        MinPhi = cms.double(0.785398163),
        MaxPhi = cms.double(0.785398163), # fixed
        MinE   = cms.double(100.0),
        MaxE   = cms.double(100.0),        # fixed
    ),
    Verbosity = cms.untracked.int32(0),
    AddAntiParticle = cms.bool(False)
)


# process.g4SimHits.FileNameGDML = 'cmsRun4D110.gdml'
##process.g4SimHits.G4Commands = ['/tracking/verbose 1']
process.g4SimHits.Physics.type = 'SimG4Core/Physics/FTFP_BERT_EMH'
##
process.RandomNumberGeneratorService.g4SimHits.initialSeed = 1234
#process.g4SimHits.Physics.G4HepEmActive = G4HEPEM
process.g4SimHits.Physics.Verbosity = cms.untracked.int32(1)
##
# allow log messages from the physics list
process.MessageLogger.PhysicsList = dict()
#
process.g4SimHits.HCalSD.TestNumberingScheme = False
process.g4SimHits.EnableHGCalSubregions = True
process.CaloSimHitStudy.TestNumbering = False
process.MessageLogger.G4cout = dict()
process.MessageLogger.SimG4CoreApplication = dict()


from Configuration.AlCa.GlobalTag import GlobalTag
if RUN4:
    ## ------------------------------------------------------------------
    ## RUN4
    process.TFileService = cms.Service('TFileService',
          fileName = cms.string('phase2_Run4D110.gg2ttbar.14TeV.128k.t32.G4HepEm.root')
    )
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')
else:
    ## ------------------------------------------------------------------
    ## RUN3
    process.TFileService = cms.Service('TFileService',
        fileName = cms.string('run3_dd4hep.2024.gg2ttbar.13p6TeV.128k.t32.G4HepEm.root')
    )
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_realistic', '')
    ## ------------------------------------------------------------------



# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.analysis_step   = cms.EndPath(process.CaloSimHitStudy)


# Schedule definition
process.schedule = cms.Schedule(process.generation_step, process.simulation_step, process.analysis_step)
for path in process.paths:
    getattr(process, path).insert(0, process.generator)


