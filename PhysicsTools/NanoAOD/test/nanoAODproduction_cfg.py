import FWCore.ParameterSet.Config     as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.StandardSequences.Eras import eras

#### Importing parameters
options = VarParsing.VarParsing()
options.register('GlobalTag',
                 "", # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global tag")

options.register('Era',
                 "", # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Era to be considered")

options.register('IsData',
                 False, # default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Is this a data dataset or a MC one")

options.parseArguments()

if options.Era == "" or options.GlobalTag == "":
    raise RuntimeError("ERROR: global tag or/and era could not be retrieved from configuration json files.")

# Creating process
process = cms.Process('NANO', getattr(eras, options.Era.split(",")[0]), getattr(eras, options.Era.split(",")[1]))

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('PhysicsTools.NanoAOD.nano_cff')

if options.IsData:
    process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
else:
    process.load('Configuration.StandardSequences.MagneticField_cff')
    process.load('SimGeneral.MixingModule.mixNoPU_cfi')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# Input source

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:./FC05E0FF-4737-E811-BE05-FA163EE2974F.root'),
    fileNames = cms.untracked.vstring('file:./EC827C0E-B534-E811-8895-008CFAEA2564.root'),
    secondaryFileNames = cms.untracked.vstring()
)
process.options = cms.untracked.PSet(
)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('nanoAOD production'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)


# Output definition
if options.IsData:
    process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(9),
        dataset = cms.untracked.PSet(
            dataTier = cms.untracked.string('NANOAOD'),
            filterName = cms.untracked.string('')
        ),
        fileName = cms.untracked.string('nanoAODoutput.root'),
        outputCommands = process.NANOAODEventContent.outputCommands
    )
else:
    process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(9),
        dataset = cms.untracked.PSet(
            dataTier = cms.untracked.string("NANOAODSIM"),
            filterName = cms.untracked.string('')
        ),
        fileName = cms.untracked.string('nanoAODoutput.root'),
        outputCommands = process.NANOAODSIMEventContent.outputCommands,
    )


# Other statements
from Configuration.AlCa.GlobalTag         import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GlobalTag, '')


# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequence if options.IsData else process.nanoSequenceMC)
process.endjob_step  = cms.EndPath(process.endOfProcess)
if options.IsData: process.NANOAODoutput_step    = cms.EndPath(process.NANOAODoutput)
else:              process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)


# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step, process.endjob_step, process.NANOAODoutput_step if options.IsData else process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers  import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


#### Customisation of the process
from PhysicsTools.NanoAOD.nano_cff        import nanoAOD_customizeMC, nanoAOD_customizeData
if options.IsData: process = nanoAOD_customizeData(process)
else:              process = nanoAOD_customizeMC(process)
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))


# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
