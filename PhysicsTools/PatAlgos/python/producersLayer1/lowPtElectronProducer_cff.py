import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi import *

sourceElectrons = cms.InputTag("lowPtGsfElectrons")

lowPtElectronMatch = electronMatch.clone(
   src = sourceElectrons,
)

patLowPtElectrons = patElectrons.clone(

    # Input collection
    electronSource = sourceElectrons,

    # MC matching
    addGenMatch      = cms.bool(True),
    embedGenMatch    = cms.bool(False),
    genParticleMatch = cms.InputTag("lowPtElectronMatch"),

    # Electron ID
    addElectronID = cms.bool(True),
    electronIDSources = cms.PSet(
        unbiased = cms.InputTag("rekeyLowPtGsfElectronSeedValueMaps:unbiased"),
        ptbiased = cms.InputTag("rekeyLowPtGsfElectronSeedValueMaps:ptbiased"),
        ID       = cms.InputTag("lowPtGsfElectronID"),
    ),

    # User data to add
    userData = cms.PSet(
        userClasses = cms.PSet( src = cms.VInputTag('') ),
        userFloats = cms.PSet( src = cms.VInputTag('') ),
        userInts = cms.PSet( src = cms.VInputTag('') ),
        userCands = cms.PSet( src = cms.VInputTag('') ),
        userFunctions = cms.vstring(),
        userFunctionLabels = cms.vstring()
    ),

    # Embedding of RECO/AOD items
    embedTrack                  = cms.bool(True),
    embedGsfElectronCore        = cms.bool(True),
    embedGsfTrack               = cms.bool(True),
    embedSuperCluster           = cms.bool(True),
    embedSeedCluster            = cms.bool(True),
    embedBasicClusters          = cms.bool(True),
    embedPreshowerClusters      = cms.bool(False),
    embedRecHits                = cms.bool(False),
    embedPflowSuperCluster      = cms.bool(False),
    embedPflowBasicClusters     = cms.bool(False),
    embedPflowPreshowerClusters = cms.bool(False),
    embedPFCandidate            = cms.bool(False),

    # Miscellaneous flags
    useParticleFlow         = cms.bool(False),
    usePfCandidateMultiMap  = cms.bool(False),
    addMVAVariables         = cms.bool(False),
    addEfficiencies         = cms.bool(False),
    addResolutions          = cms.bool(False),
    embedHighLevelSelection = cms.bool(False),
    computeMiniIso          = cms.bool(False),
    addPFClusterIso         = cms.bool(False),
    addPuppiIsolation       = cms.bool(False),
    isoDeposits             = cms.PSet(),
    userIsolation           = cms.PSet(),
    isolationValues         = cms.PSet(),
    isolationValuesNoPFId   = cms.PSet(),
)

makePatLowPtElectronsTask = cms.Task(
    lowPtElectronMatch,
    patLowPtElectrons,
    )

makePatLowPtElectrons = cms.Sequence(makePatLowPtElectronsTask)


# Modifiers
from Configuration.Eras.Modifier_fastSim_cff import fastSim
fastSim.toModify(patLowPtElectrons, embedTrack = False)
from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
from Configuration.Eras.Modifier_run2_miniAOD_94XFall17_cff import run2_miniAOD_94XFall17
(run2_miniAOD_80XLegacy | run2_miniAOD_94XFall17).toModify(patLowPtElectrons,
                                                           electronSource = "gedGsfElectrons",
                                                           genParticleMatch = "electronMatch"
                                                           )

# Schedule rekeying of seed BDT ValueMaps by reco::GsfElectron for run2_miniAOD_UL
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronSeedValueMaps_cff import rekeyLowPtGsfElectronSeedValueMaps
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cfi import lowPtGsfElectronID
makePatLowPtElectronsTask_ = makePatLowPtElectronsTask.copy()
makePatLowPtElectronsTask_.add(rekeyLowPtGsfElectronSeedValueMaps)
makePatLowPtElectronsTask_.add(lowPtGsfElectronID)
run2_miniAOD_UL.toReplaceWith(makePatLowPtElectronsTask,makePatLowPtElectronsTask_)

# Apply the energy regression for UL
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronFinalizer_cfi import lowPtGsfElectrons as _lowPtGsfElectrons
lowPtGsfElectronsPostRegression = _lowPtGsfElectrons.clone(
    previousGsfElectronsTag = sourceElectrons
)

from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
sourceElectronsPostRegression = cms.InputTag("lowPtGsfElectronsPostRegression")
run2_miniAOD_UL.toModify(lowPtGsfLinks, electrons = sourceElectronsPostRegression)
run2_miniAOD_UL.toModify(lowPtElectronMatch, src = sourceElectronsPostRegression)
run2_miniAOD_UL.toModify(patLowPtElectrons, electronSource = sourceElectronsPostRegression)
run2_miniAOD_UL.toModify(rekeyLowPtGsfElectronSeedValueMaps, gsfElectrons = sourceElectronsPostRegression)

makePatLowPtElectronsTaskWithRegression = cms.Task(lowPtGsfElectronsPostRegression, makePatLowPtElectronsTask.copy())
run2_miniAOD_UL.toReplaceWith(makePatLowPtElectronsTask, makePatLowPtElectronsTaskWithRegression)
