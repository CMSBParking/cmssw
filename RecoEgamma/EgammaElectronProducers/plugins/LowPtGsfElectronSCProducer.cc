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

LowPtGsfElectronSCProducer::LowPtGsfElectronSCProducer( const edm::ParameterSet& cfg ) :
  gsfPfRecTracks_{consumes<reco::GsfPFRecTrackCollection>( cfg.getParameter<edm::InputTag>("gsfPfRecTracks") )},
  ecalClusters_{consumes<reco::PFClusterCollection>( cfg.getParameter<edm::InputTag>("ecalClusters") )}
{
  produces<reco::CaloClusterCollection>();
  produces<reco::SuperClusterCollection>();
  produces< edm::ValueMap<reco::SuperClusterRef> >();
}

LowPtGsfElectronSCProducer::~LowPtGsfElectronSCProducer()
{}

void LowPtGsfElectronSCProducer::produce( edm::Event& event, const edm::EventSetup& setup )
{

  // Input GsfPFRecTracks collection
  edm::Handle<reco::GsfPFRecTrackCollection> gsfPfRecTracks;
  event.getByToken(gsfPfRecTracks_,gsfPfRecTracks);
  if ( !gsfPfRecTracks.isValid() ) { edm::LogError("Problem with gsfPfRecTracks handle"); }

  // Input EcalClusters collection
  edm::Handle<reco::PFClusterCollection> ecalClusters;
  event.getByToken(ecalClusters_,ecalClusters);
  if ( !ecalClusters.isValid() ) { edm::LogError("Problem with ecalClusters handle"); }

  // Output SuperClusters collection and getRefBeforePut
  auto superClusters = std::make_unique<reco::SuperClusterCollection>( reco::SuperClusterCollection() );
  superClusters->reserve(gsfPfRecTracks->size());
  const reco::SuperClusterRefProd superClustersRefProd = event.getRefBeforePut<reco::SuperClusterCollection>();

  // Output ValueMap container of SuperClusterRef to GsfPFRecTrackRef index
  std::vector<reco::SuperClusterRef> superClustersValueMap;

  // Output CaloClusters collection
  auto caloClusters = std::make_unique<reco::CaloClusterCollection>( reco::CaloClusterCollection() );
  caloClusters->reserve(ecalClusters->size());
  //const reco::SuperClusterRefProd caloClustersRefProd = event.getRefBeforePut<reco::SuperClusterCollection>();

  // Index for "best" CaloCluster per GSF track per trajectory point
  std::vector< std::vector<int> > cluster_idx;
  cluster_idx.resize( gsfPfRecTracks->size(), std::vector<int>() );

  // Index for "best" PFCluster per GSF track per trajectory point
  std::vector< std::vector<int> > pfcluster_idx;
  pfcluster_idx.resize( gsfPfRecTracks->size(), std::vector<int>() );

  // dr2min for "best" cluster per GSF track per trajectory point
  std::vector< std::vector<float> > cluster_dr2min;
  cluster_dr2min.resize( gsfPfRecTracks->size(), std::vector<float>() );

  // Construct list of trajectory points from the GSF track and electron brems
  std::vector< std::vector<const reco::PFTrajectoryPoint*> > points;
  points.resize( gsfPfRecTracks->size(), std::vector<const reco::PFTrajectoryPoint*>() );
  for ( size_t itrk = 0; itrk < gsfPfRecTracks->size(); ++itrk ) { 
    reco::GsfPFRecTrackRef trk(gsfPfRecTracks,itrk);
    points[itrk].push_back( &trk->extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax) );
    for ( auto brem : trk->PFRecBrem() ) {
      points[itrk].push_back( &brem.extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax) ); 
    }
    cluster_idx[itrk].resize(points[itrk].size(),-1);
    pfcluster_idx[itrk].resize(points[itrk].size(),-1);
    cluster_dr2min[itrk].resize(points[itrk].size(),1.);
  }

  // Loop over clusters
  for ( size_t iclu = 0; iclu < ecalClusters->size(); ++iclu ) {
    std::pair<int,int> point = std::make_pair(-1,-1);
    float dr2min = 1.;
    // Loop over nested vector of points
    for ( size_t ipoint = 0; ipoint < points.size(); ++ipoint ) {
      for ( size_t jpoint = 0; jpoint < points[ipoint].size(); ++jpoint ) {
	if ( points[ipoint][jpoint]->isValid() ) {
	  float dr2 = reco::deltaR2( ecalClusters->at(iclu), 
				     points[ipoint][jpoint]->positionREP() );
	  if ( dr2 < dr2min ) {
	    // Store nearest point to this cluster
	    dr2min = dr2;
	    point = std::make_pair(ipoint,jpoint);
	  }
	}
      }
    }
    if ( point.first >= 0 && 
	 point.second >= 0 && 
	 dr2min < cluster_dr2min[point.first][point.second] ) {
      // Copy CaloCluster to new collection
      caloClusters->push_back(ecalClusters->at(iclu));
      // Store index for creation of SC later
      cluster_idx[point.first][point.second] = caloClusters->size()-1;
      pfcluster_idx[point.first][point.second] = iclu;
      cluster_dr2min[point.first][point.second] = dr2min;
    }
  }

  // Put CaloClusters in event and get orphan handle
  const edm::OrphanHandle<reco::CaloClusterCollection>& caloClustersH = event.put(std::move(caloClusters));

  // Loop through GSF tracks
  for ( size_t itrk = 0; itrk < gsfPfRecTracks->size(); ++itrk ) { 
    
    // Used to create SC
    float energy = 0.;
    float X = 0., Y = 0., Z = 0.;
    reco::CaloClusterPtr seed;
    reco::CaloClusterPtrVector clusters;
    std::vector<const reco::PFCluster*> barePtrs;

    // Loop through clusters associated with GSF track (via points)
    for ( size_t iclu = 0; iclu < cluster_idx[itrk].size(); ++iclu ) { 
      if ( cluster_idx[itrk][iclu] < 0 ) { continue; }
      reco::CaloClusterPtr clu(caloClustersH,cluster_idx[itrk][iclu]);
      if ( clu.isNull() ) { continue; }
      if ( seed.isNull() ) { seed = clu; }
      clusters.push_back(clu);
      energy += clu->correctedEnergy();
      X += clu->position().X() * clu->correctedEnergy();
      Y += clu->position().Y() * clu->correctedEnergy();
      Z += clu->position().Z() * clu->correctedEnergy();
      reco::PFClusterRef pfclu(ecalClusters,pfcluster_idx[itrk][iclu]);
      if ( pfclu.isNonnull() ) { barePtrs.push_back(&*pfclu); }
    }
    X /= energy;
    Y /= energy;
    Z /= energy;

    // Create SC
    if ( seed.isNonnull() ) {
      reco::SuperCluster sc(energy,math::XYZPoint(X,Y,Z));
      sc.setCorrectedEnergy(energy);
      sc.setSeed(seed);
      std::cout << "ok " << itrk << " " << clusters.size() << " " << energy;// << std::endl;
      sc.setClusters(clusters);
      for ( const auto clu : clusters ) { sc.addCluster(clu); }
      PFClusterWidthAlgo pfwidth(barePtrs);
      sc.setEtaWidth(pfwidth.pflowEtaWidth());
      sc.setPhiWidth(pfwidth.pflowPhiWidth());
      sc.rawEnergy(); // Cache the value of raw energy
      std::cout << " ok " << sc.clusters().size() << " " << sc.correctedEnergy() << std::endl;
      superClusters->push_back(sc);

      // Populate ValueMap container
      superClustersValueMap.push_back( reco::SuperClusterRef(superClustersRefProd,itrk) );
    } else {
      std::cout << "missing " << itrk << std::endl;
      superClustersValueMap.push_back( reco::SuperClusterRef(superClustersRefProd.id()) );
    }

  } // GSF tracks

  // Put SuperClusters in event
  event.put(std::move(superClusters));

  auto ptr = std::make_unique< edm::ValueMap<reco::SuperClusterRef> >( edm::ValueMap<reco::SuperClusterRef>() );
  edm::ValueMap<reco::SuperClusterRef>::Filler filler(*ptr);
  filler.insert(gsfPfRecTracks, superClustersValueMap.begin(), superClustersValueMap.end());
  filler.fill();
  event.put(std::move(ptr));

  //@@

//  // Temporary map of CaloClusterPtr to CaloClusterRef index
//  std::map<reco::CaloClusterPtr,unsigned int> caloClustersMap;
//
//  // Iterate through GsfPfRecTracks and create corresponding SuperClusters
//  std::vector<int> matchedClusters;
//  for ( size_t igsfpf = 0; igsfpf < gsfPfRecTracks->size(); ++igsfpf ) { 
//
//    // Refs to GSF PF tracks
//    reco::GsfPFRecTrackRef gsfpf(gsfPfRecTracks, igsfpf);
//
//    // Temp PFClusters collection to build SC
//    std::vector<reco::PFClusterRef> tmpClusters;
//
//    // Find closest "seed cluster" to GSF track extrapolated to ECAL
//    const reco::PFTrajectoryPoint& point1 = gsfpf->extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax);
//    reco::PFClusterRef best_seed = closestCluster( point1, ecalClusters, matchedClusters );
//    if ( best_seed.isNonnull() ) { 
//      tmpClusters.push_back(best_seed);
//      reco::CaloClusterPtr ptr(edm::refToPtr(best_seed));
//      if ( !caloClustersMap.count(ptr) ) {
//	caloClusters->push_back(*ptr); // Copy CaloCluster
//	caloClustersMap[ptr] = caloClusters->size() - 1;
//      }
//    }
//    
//    // Find closest "brem cluster" using brem trajectory extrapolated to ECAL
//    const std::vector<reco::PFBrem>& brems = gsfpf->PFRecBrem();
//    for ( auto brem : brems ) {
//      const reco::PFTrajectoryPoint& point2 = brem.extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax);
//      reco::PFClusterRef best_brem = closestCluster( point2, ecalClusters, matchedClusters );
//      if ( best_brem.isNonnull() ) { 
//	tmpClusters.push_back(best_brem);
//	if ( best_seed.isNull() ) { best_seed = best_brem; } // Use brem as seed
//	reco::CaloClusterPtr ptr(edm::refToPtr(best_brem));
//	caloClusters->push_back(*ptr); // Copy CaloCluster 
//	caloClustersMap[ptr] = caloClusters->size() - 1;
//      }
//    }
//    
//    // If all else fails, attempt to extrapolate KF track and match to seed PF cluster
//    if ( best_seed.isNull() ) { 
//      const reco::PFRecTrackRef& kfpf = gsfpf->kfPFRecTrackRef();
//      const reco::PFTrajectoryPoint& point3 = kfpf->extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax);
//      reco::PFClusterRef best_kf = closestCluster( point3, ecalClusters, matchedClusters );
//      if ( best_kf.isNonnull() ) { 
//	tmpClusters.push_back(best_kf);
//	best_seed = best_kf; // Use KF as seed
//	reco::CaloClusterPtr ptr(edm::refToPtr(best_kf));
//	caloClusters->push_back(*ptr); // Copy CaloCluster
//	caloClustersMap[ptr] = caloClusters->size() - 1;
//      }
//    }
//    
//    // Now make the SuperCluster
//    if ( !best_seed.isNull() ) {
//
//      float posX=0.,posY=0.,posZ=0.;
//      float scEnergy=0.;
//      for ( const auto clus : tmpClusters ) {
//	scEnergy+=clus->correctedEnergy();
//	posX+=clus->correctedEnergy()*clus->position().X();
//	posY+=clus->correctedEnergy()*clus->position().Y();
//	posZ+=clus->correctedEnergy()*clus->position().Z();
//      }
//      posX/=scEnergy;
//      posY/=scEnergy;
//      posZ/=scEnergy;
//      reco::SuperCluster sc(scEnergy,math::XYZPoint(posX,posY,posZ));
//      sc.setCorrectedEnergy(scEnergy);
//      sc.setSeed(edm::refToPtr(best_seed));
//      std::vector<const reco::PFCluster*> barePtrs;
//      for ( const auto clus : tmpClusters ) {
//	sc.addCluster(edm::refToPtr(clus));
//	barePtrs.push_back(&*clus);
//      }
//      PFClusterWidthAlgo pfwidth(barePtrs);
//      sc.setEtaWidth(pfwidth.pflowEtaWidth());
//      sc.setPhiWidth(pfwidth.pflowPhiWidth());
//      sc.rawEnergy(); // Cache the value of raw energy
//
//      // Store new SuperCluster 
//      superClusters->push_back(sc);
//      
//      // Populate ValueMap container
//      superClustersValueMap.push_back( reco::SuperClusterRef(superClustersRefProd,igsfpf) );
//    } else {
//      superClustersValueMap.push_back( reco::SuperClusterRef(superClustersRefProd.id()) );
//    }
//  }
//
//  // Put CaloClusters in event first
//  const edm::OrphanHandle<reco::CaloClusterCollection>& caloClustersH = event.put(std::move(caloClusters));
//
//  // Update CaloClusterRefs in SuperClusters
//  for ( auto& sc : *superClusters ) {
//    sc.setSeed( reco::CaloClusterPtr( caloClustersH, caloClustersMap[sc.seed()] ) );
//    reco::CaloClusterPtrVector clusters;
//    for ( auto clu : sc.clusters() ) {
//      clusters.push_back( reco::CaloClusterPtr( caloClustersH, caloClustersMap[clu] ) );
//    }
//    sc.setClusters(clusters);
//  }
//
//  // Put SuperClusters in event
//  event.put(std::move(superClusters));
//
//  // Put ValueMap<SuperClusterRef> in event
//  auto ptr = std::make_unique< edm::ValueMap<reco::SuperClusterRef> >( edm::ValueMap<reco::SuperClusterRef>() );
//  edm::ValueMap<reco::SuperClusterRef>::Filler filler(*ptr);
//  filler.insert(gsfPfRecTracks, superClustersValueMap.begin(), superClustersValueMap.end());
//  filler.fill();
//  event.put(std::move(ptr));

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



//@@
//      // Construct list of trajectory points from the GSF track and electron brems
//      std::vector<const reco::PFTrajectoryPoint*> points;
//      points.push_back( &gsfpf->extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax) );
//      const std::vector<reco::PFBrem>& brems = gsfpf->PFRecBrem();
//      for ( auto brem : brems ) {
//	points.push_back( &brem.extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax) ); 
//      }
//
//      // "Best" cluster index and dr2min per trajectory point
//      std::vector<float> cluster_idx(points.size(),-1);
//      std::vector<float> cluster_dr2min(points.size(),1.);
//      std::cout << "points "
//		<< points.size() << " " 
//		<< cluster_idx.size() << " " 
//		<< cluster_dr2min.size() << " " 
//		<< std::endl;
//
//      // Loop over clusters
//      for ( size_t iclu = 0; iclu < ecalClusters->size(); ++iclu ) {
//	int point = -1;
//	float dr2min = 1.;
//	for ( size_t ipoint = 0; ipoint < points.size(); ++ipoint ) {
//	  if( points[ipoint]->isValid() ) {
//	    float dr2 = reco::deltaR2( ecalClusters->at(iclu), points[ipoint]->positionREP() );
//	    if ( dr2 < dr2min ) {
//	      dr2min = dr2;
//	      point = ipoint;
//	    }
//	  }
//	}
//	// Store cluster index if matched and closer than the current best (if any)
//	if ( point >= 0 && dr2min < cluster_dr2min[point] ) {
//	  cluster_idx[point] = iclu;
//	  cluster_dr2min[point] = dr2min;
//	}
//      }
//
//      // Loop through trajectory points
//      for ( size_t ipoint = 0; ipoint < points.size(); ++ipoint ) {
//	if ( cluster_idx[ipoint] >= 0 ) {
//	  // Store ClusterRef to form SC 
//	  tmpClusters.emplace_back(ecalClusters,cluster_idx[ipoint]);
//	  // Identify first CaloCluster found as seed
//	  if ( best_seed.isNull() ) { best_seed = tmpClusters.back(); }
//	  // Copy of CaloCluster for new collection 
//	  reco::CaloClusterPtr ptr(edm::refToPtr(tmpClusters.back()));
//	  if ( !caloClustersMap.count(ptr) ) {
//	    caloClusters->push_back(*ptr); 
//	    caloClustersMap[ptr] = caloClusters->size() - 1;
//	  }
//	}
//	std::cout << "ipoint "
//		  << ipoint << " " 
//		  << cluster_idx[ipoint] << " " 
//		  << cluster_dr2min[ipoint] << " " 
//		  << std::endl;
//      }
//      std::cout << "end "
//		<< tmpClusters.size() << " " 
//		<< caloClusters->size() << " " 
//		<< caloClustersMap.size() << " " 
//		<< best_seed.isNull()
//		<< std::endl;
