import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")

process.load("Configuration.Generator.TTbar_14TeV_TuneCP5_cfi")

# ---------------------------------------------------
# Número de eventos
# ---------------------------------------------------
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000)
)

process.source = cms.Source("EmptySource")

# ---------------------------------------------------
# Semilla fija para reproducibilidad
# ---------------------------------------------------
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    generator = cms.PSet(
        initialSeed = cms.untracked.uint32(12345),
        engineName  = cms.untracked.string('HepJamesRandom')
    )
)

# ---------------------------------------------------
# Output HepMC
# ---------------------------------------------------
process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string("ttbar_GEN.root"),
    outputCommands = cms.untracked.vstring(
        "keep *_generator_*_*"
    )
)

process.generation_step = cms.Path(process.ProductionFilterSequence)
process.outpath = cms.EndPath(process.output)

# Needed for multithreading
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)
