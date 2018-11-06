#ifndef LowPtGsfElectronCoreProducer_h
#define LowPtGsfElectronCoreProducer_h

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCoreFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "RecoEgamma/EgammaElectronProducers/plugins/GsfElectronCoreBaseProducer.h"

class LowPtGsfElectronCoreProducer : public GsfElectronCoreBaseProducer {

 public:
  
  explicit LowPtGsfElectronCoreProducer( const edm::ParameterSet& conf );
  
  ~LowPtGsfElectronCoreProducer() override;
  
  void produce( edm::Event&, const edm::EventSetup& ) override;
  
 private:
  
  edm::EDGetTokenT<reco::SuperClusterCollection> superClusters_;
  edm::EDGetTokenT< edm::ValueMap<reco::SuperClusterRef> > superClusterRefs_;

};

#endif // LowPtGsfElectronCoreProducer_h

