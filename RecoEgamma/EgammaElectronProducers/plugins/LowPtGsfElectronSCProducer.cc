#include "RecoEgamma/EgammaElectronProducers/plugins/LowPtGsfElectronSCProducer.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFClusterWidthAlgo.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <iostream>

LowPtGsfElectronSCProducer::LowPtGsfElectronSCProducer( const edm::ParameterSet& cfg )
{
  produces<reco::CaloClusterCollection>();
  produces<reco::SuperClusterCollection>();
  produces< edm::ValueMap<reco::SuperClusterRef> >();
  gsfPfRecTracks_ = consumes<reco::GsfPFRecTrackCollection>( cfg.getParameter<edm::InputTag>("gsfPfRecTracks") );
  ecalClusters_ = consumes<reco::PFClusterCollection>( cfg.getParameter<edm::InputTag>("ecalClusters") );
}

LowPtGsfElectronSCProducer::~LowPtGsfElectronSCProducer()
{}

void LowPtGsfElectronSCProducer::produce( edm::Event& event, const edm::EventSetup& setup )
{


  // Input EcalClusters collection
  edm::Handle<reco::PFClusterCollection> ecalClusters;
  event.getByToken(ecalClusters_,ecalClusters);
  if ( !ecalClusters.isValid() ) { edm::LogError("Problem with ecalClusters handle"); }

  // Input GsfPFRecTracks collection
  edm::Handle<reco::GsfPFRecTrackCollection> gsfPfRecTracks;
  event.getByToken(gsfPfRecTracks_,gsfPfRecTracks);
  if ( !gsfPfRecTracks.isValid() ) { edm::LogError("Problem with gsfPfRecTracks handle"); }

  // Output SuperClusters collection and getRefBeforePut
  auto superClusters = std::make_unique<reco::SuperClusterCollection>( reco::SuperClusterCollection() );
  superClusters->reserve(gsfPfRecTracks->size());
  const reco::SuperClusterRefProd superClustersRefProd = event.getRefBeforePut<reco::SuperClusterCollection>();

  // Output ValueMap container of SuperClusterRef to GsfPFRecTrackRef index
  std::vector<reco::SuperClusterRef> superClustersValueMap;

  // Output CaloClusters collection
  auto caloClusters = std::make_unique<reco::CaloClusterCollection>( reco::CaloClusterCollection() );
  caloClusters->reserve(ecalClusters->size());

  // Temporary map of CaloClusterPtr to CaloClusterRef index
  std::map<reco::CaloClusterPtr,unsigned int> caloClustersMap;

  // Iterate through GsfPfRecTracks and create corresponding SuperClusters
  std::vector<int> matchedClusters;
  for ( size_t igsfpf = 0; igsfpf < gsfPfRecTracks->size(); ++igsfpf ) { 

    // Refs to GSF PF tracks
    reco::GsfPFRecTrackRef gsfpf(gsfPfRecTracks, igsfpf);

    // Temp PFClusters collection to build SC
    std::vector<reco::PFClusterRef> tmpClusters;

    // Find closest "seed cluster" to GSF track extrapolated to ECAL
    const reco::PFTrajectoryPoint& point1 = gsfpf->extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax);
    reco::PFClusterRef best_seed = closestCluster( point1, ecalClusters, matchedClusters );
    if ( best_seed.isNonnull() ) { 
      tmpClusters.push_back(best_seed);
      reco::CaloClusterPtr ptr(edm::refToPtr(best_seed));
      if ( !caloClustersMap.count(ptr) ) {
	caloClusters->push_back(*ptr); // Copy CaloCluster
        caloClustersMap[ptr] = caloClusters->size() - 1;
      }
    }
    
    // Find closest "brem cluster" using brem trajectory extrapolated to ECAL
    const std::vector<reco::PFBrem>& brems = gsfpf->PFRecBrem();
    for ( auto brem : brems ) {
      const reco::PFTrajectoryPoint& point2 = brem.extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax);
      reco::PFClusterRef best_brem = closestCluster( point2, ecalClusters, matchedClusters );
      if ( best_brem.isNonnull() ) { 
	tmpClusters.push_back(best_brem);
	if ( best_seed.isNull() ) { best_seed = best_brem; } // Use brem as seed
	reco::CaloClusterPtr ptr(edm::refToPtr(best_brem));
	caloClusters->push_back(*ptr); // Copy CaloCluster 
        caloClustersMap[ptr] = caloClusters->size() - 1;
      }
    }
  
    // If all else fails, attempt to extrapolate KF track and match to seed PF cluster
    reco::PFTrajectoryPoint point3;
    if ( best_seed.isNull() ) { 
      const reco::PFRecTrackRef& kfpf = gsfpf->kfPFRecTrackRef();
      point3 = kfpf->extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax);
      reco::PFClusterRef best_kf = closestCluster( point3, ecalClusters, matchedClusters );
      if ( best_kf.isNonnull() ) { 
	tmpClusters.push_back(best_kf);
	best_seed = best_kf; // Use KF as seed
	reco::CaloClusterPtr ptr(edm::refToPtr(best_kf));
	caloClusters->push_back(*ptr); // Copy CaloCluster 
        caloClustersMap[ptr] = caloClusters->size() - 1;
      }
    }

    // Now make the SuperCluster
    if ( !best_seed.isNull() ) {

      float posX=0.,posY=0.,posZ=0.;
      float scEnergy=0.;
      for ( const auto clus : tmpClusters ) {
	scEnergy+=clus->correctedEnergy();
	posX+=clus->correctedEnergy()*clus->position().X();
	posY+=clus->correctedEnergy()*clus->position().Y();
	posZ+=clus->correctedEnergy()*clus->position().Z();
      }
      posX/=scEnergy;
      posY/=scEnergy;
      posZ/=scEnergy;
      reco::SuperCluster sc(scEnergy,math::XYZPoint(posX,posY,posZ));
      sc.setCorrectedEnergy(scEnergy);
      sc.setSeed(edm::refToPtr(best_seed));
      std::vector<const reco::PFCluster*> barePtrs;
      for ( const auto clus : tmpClusters ) {
	sc.addCluster(edm::refToPtr(clus));
	barePtrs.push_back(&*clus);
      }
      PFClusterWidthAlgo pfwidth(barePtrs);
      sc.setEtaWidth(pfwidth.pflowEtaWidth());
      sc.setPhiWidth(pfwidth.pflowPhiWidth());
      sc.rawEnergy(); // Cache the value of raw energy

      // Store new SuperCluster 
      superClusters->push_back(sc);
      
      // Populate ValueMap container
      superClustersValueMap.push_back( reco::SuperClusterRef(superClustersRefProd,igsfpf) );
    } else {
      superClustersValueMap.push_back( reco::SuperClusterRef(superClustersRefProd.id()) );
    }
  }

  // Put CaloClusters in event first
  const edm::OrphanHandle<reco::CaloClusterCollection>& caloClustersH = event.put(std::move(caloClusters));

  // Update CaloClusterRefs in SuperClusters
  for ( auto& sc : *superClusters ) {
    sc.setSeed( reco::CaloClusterPtr( caloClustersH, caloClustersMap[sc.seed()] ) );
    reco::CaloClusterPtrVector clusters;
    for ( auto clu : sc.clusters() ) {
      clusters.push_back( reco::CaloClusterPtr( caloClustersH, caloClustersMap[clu] ) );
    }
    sc.setClusters(clusters);
  }

  // Put SuperClusters in event
  event.put(std::move(superClusters));

  // Put ValueMap<SuperClusterRef> in event
  auto ptr = std::make_unique< edm::ValueMap<reco::SuperClusterRef> >( edm::ValueMap<reco::SuperClusterRef>() );
  edm::ValueMap<reco::SuperClusterRef>::Filler filler(*ptr);
  filler.insert(gsfPfRecTracks, superClustersValueMap.begin(), superClustersValueMap.end());
  filler.fill();
  event.put(std::move(ptr));

}

reco::PFClusterRef LowPtGsfElectronSCProducer::closestCluster( const reco::PFTrajectoryPoint& point,
							       const edm::Handle<reco::PFClusterCollection>& clusters,
							       std::vector<int>& matched ) {
  reco::PFClusterRef closest;
  if ( point.isValid() ) {
    float dr2min = 1.e6;
    for ( size_t ii = 0; ii < clusters->size(); ++ii ) {
      if ( std::find( matched.begin(), matched.end(), ii ) == matched.end() ) {
	float dr2 = reco::deltaR2( clusters->at(ii), point.positionREP() );
	if ( dr2 < dr2min ) {
	  closest = reco::PFClusterRef( clusters, ii );
	  dr2min = dr2;
	}
      }
    }
    if ( dr2min < 1.e5 ) { matched.push_back( closest.index() ); }
  }
  return closest;
}

void LowPtGsfElectronSCProducer::fillDescription( edm::ParameterSetDescription& desc ) 
{
  desc.add<edm::InputTag>("gsfPfRecTracks",edm::InputTag("lowPtGsfElePfGsfTracks"));
  desc.add<edm::InputTag>("ecalClusters",edm::InputTag("particleFlowClusterECAL"));
}
