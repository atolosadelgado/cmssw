import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process("SIM", Phase2C17I13M9)

# Basic services
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.Timing = cms.Service('Timing')

# Geometry, B-field
process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T21", "")

# Input GEN file
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:ttbar_GEN.root")
)

# Load the standard SIM sequence (SimIdeal) which in turn configures g4SimHits
process.load("Configuration.StandardSequences.SimIdeal_cff")

# FIX: ensure HepMC product labels have the correct types for this release
process.g4SimHits.Generator.HepMCProductLabel = cms.InputTag("generator", "unsmeared")
process.g4SimHits.Generator.HepMCProductLabel2 = cms.InputTag("")

# Physics list / verbosity and SD toggles (optional)
process.g4SimHits.Physics.type = cms.string("SimG4Core/Physics/FTFP_BERT_EMH")
process.g4SimHits.Physics.Verbosity = cms.untracked.int32(0)
process.g4SimHits.HCalSD.TestNumberingScheme = cms.bool(False)
process.g4SimHits.EnableHGCalSubregions = True

# Do NOT load obsolete SimFastTiming modules here.

# Execute only the simulation (no output)
process.simulation_step = cms.Path(process.g4SimHits)
process.schedule = cms.Schedule(process.simulation_step)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    numberOfThreads = cms.untracked.uint32(4)
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
