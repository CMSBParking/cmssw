import FWCore.ParameterSet.Config as cms

#==============================================================================
# Sequence to make low pT electrons.
# Several steps are cloned and modified to 'open up' the sequence:
# tracker-driven electron seeds, KF track candidates, GSF tracks.
#==============================================================================

# PFTracks (input: generalTracks module)
from RecoParticleFlow.PFTracking.pfTrack_cfi import *
lowPtGsfElePfTracks = pfTrack.clone()
lowPtGsfElePfTracks.TkColList = ['generalTracks']
lowPtGsfElePfTracks.GsfTracksInEvents = False
lowPtGsfElePfTracks.GsfTrackModuleLabel = ''

# Low pT electron seeds
# Below relies on default configuration for generalTracks
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronSeeds_cfi import *

# Electron (KF) track candidates
from TrackingTools.GsfTracking.CkfElectronCandidateMaker_cff import *
lowPtGsfEleTrajectoryFilter = TrajectoryFilterForElectrons.clone()
lowPtGsfEleTrajectoryFilter.minPt = 0.
lowPtGsfEleTrajectoryFilter.minimumNumberOfHits = 3
lowPtGsfEleTrajectoryBuilder = TrajectoryBuilderForElectrons.clone()
lowPtGsfEleTrajectoryBuilder.trajectoryFilter.refToPSet_ = 'lowPtGsfEleTrajectoryFilter'
lowPtGsfEleCkfTrackCandidates = electronCkfTrackCandidates.clone()
lowPtGsfEleCkfTrackCandidates.TrajectoryBuilderPSet.refToPSet_ = 'lowPtGsfEleTrajectoryBuilder'
lowPtGsfEleCkfTrackCandidates.src = 'lowPtGsfElectronSeeds'

# GSF tracks
from TrackingTools.GsfTracking.GsfElectronGsfFit_cff import *
lowPtGsfEleFittingSmoother = GsfElectronFittingSmoother.clone()
lowPtGsfEleFittingSmoother.ComponentName = 'lowPtGsfEleFittingSmoother'
lowPtGsfEleFittingSmoother.MinNumberOfHits = 2
from TrackingTools.GsfTracking.GsfElectronGsfFit_cff import * 
lowPtGsfEleGsfTracks = electronGsfTracks.clone()
lowPtGsfEleGsfTracks.Fitter = 'lowPtGsfEleFittingSmoother'
lowPtGsfEleGsfTracks.src = 'lowPtGsfEleCkfTrackCandidates'

# PFTracks (input: lowPtGsfEleGsfTracks module)
#from RecoParticleFlow.PFTracking.pfTrack_cfi import *
#lowPtGsfElePfTracks = pfTrack.clone()
#lowPtGsfElePfTracks.TkColList = ['generalTracks']
#lowPtGsfElePfTracks.GsfTracksInEvents = True
#lowPtGsfElePfTracks.GsfTrackModuleLabel = 'lowPtGsfEleGsfTracks'

# PFGSFTracks
from RecoParticleFlow.PFTracking.pfTrackElec_cfi import *
lowPtGsfElePfGsfTracks = pfTrackElec.clone()
lowPtGsfElePfGsfTracks.GsfTrackModuleLabel = 'lowPtGsfEleGsfTracks'
lowPtGsfElePfGsfTracks.PFRecTrackLabel = 'lowPtGsfElePfTracks'
lowPtGsfElePfGsfTracks.applyGsfTrackCleaning = False
lowPtGsfElePfGsfTracks.useFifthStepForTrackerDrivenGsf = True

# SuperCluster generator and matching to GSF tracks
# Below relies on the following default configurations:
# RecoParticleFlow/PFClusterProducer/python/particleFlowClusterECALUncorrected_cfi.py
# RecoParticleFlow/PFClusterProducer/python/particleFlowClusterECAL_cff.py
# (particleFlowClusterECAL_cfi is generated automatically)
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronSuperClusters_cfi import *
lowPtGsfElectronSuperClusters.gsfPfRecTracks = 'lowPtGsfElePfGsfTracks'

# Low pT electron cores
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronCores_cfi import *
lowPtGsfElectronCores.gsfPfRecTracks = 'lowPtGsfElePfGsfTracks'
lowPtGsfElectronCores.gsfTracks = 'lowPtGsfEleGsfTracks'
lowPtGsfElectronCores.useGsfPfRecTracks = True

# Low pT electrons
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectrons_cfi import *
lowPtGsfElectrons.gsfElectronCoresTag = 'lowPtGsfElectronCores'
lowPtGsfElectrons.seedsTag = 'lowPtGsfElectronSeeds'
lowPtGsfElectrons.gsfPfRecTracksTag = ''
lowPtGsfElectrons.useGsfPfRecTracks = False

# Full Open sequence 
lowPtGsfElectronSequence = cms.Sequence(lowPtGsfElePfTracks+
                                        lowPtGsfElectronSeeds+
                                        lowPtGsfEleCkfTrackCandidates+
                                        lowPtGsfEleGsfTracks+
                                        #lowPtGsfElePfTracks+
                                        lowPtGsfElePfGsfTracks+
                                        lowPtGsfElectronSuperClusters+
                                        lowPtGsfElectronCores+
                                        lowPtGsfElectrons
                                        )
