#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CondFormats/GeometryObjects/interface/PGeometricDetExtra.h"
#include "Geometry/Records/interface/PGeometricDetExtraRcd.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDetExtra.h"
#include "DetectorDescription/DDCMS/interface/DDCompactView.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <vector>

class PGeometricDetExtraBuilder : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  PGeometricDetExtraBuilder(const edm::ParameterSet&);

  void beginRun(edm::Run const& iEvent, edm::EventSetup const&) override;
  void analyze(edm::Event const& iEvent, edm::EventSetup const&) override {}
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {}

private:
  void putOne(const GeometricDetExtra& gde, PGeometricDetExtra* pgde);
  bool fromDD4hep;
};

PGeometricDetExtraBuilder::PGeometricDetExtraBuilder(const edm::ParameterSet& iConfig) {
  fromDD4hep = iConfig.getUntrackedParameter<bool>("fromDD4hep", false);
}


void PGeometricDetExtraBuilder::beginRun(const edm::Run&, edm::EventSetup const& es) {
  PGeometricDetExtra* pgde = new PGeometricDetExtra;
  edm::Service<cond::service::PoolDBOutputService> mydbservice;
  if (!mydbservice.isAvailable()) {
    edm::LogError("PGeometricDetExtraBuilder") << "PoolDBOutputService unavailable";
    return;
  }

  if (!fromDD4hep) {
    edm::ESTransientHandle<DDCompactView> cpvH;
    es.get<IdealGeometryRecord>().get(cpvH);
  } else {
    edm::ESTransientHandle<cms::DDCompactView> cpvH;
    es.get<IdealGeometryRecord>().get(cpvH);
  }
  edm::ESHandle<std::vector<GeometricDetExtra> > gdeH;
  es.get<IdealGeometryRecord>().get(gdeH);
  const std::vector<GeometricDetExtra>& gdes = (*gdeH);

  std::vector<GeometricDetExtra>::const_iterator git = gdes.begin();
  std::vector<GeometricDetExtra>::const_iterator egit = gdes.end();

  for (; git != egit; ++git) {  // one level below "tracker"
    putOne(*git, pgde);
  }
  if (mydbservice->isNewTagRequest("PGeometricDetExtraRcd")) {
    mydbservice->createNewIOV<PGeometricDetExtra>(
        pgde, mydbservice->beginOfTime(), mydbservice->endOfTime(), "PGeometricDetExtraRcd");
  } else {
    edm::LogError("PGeometricDetExtraBuilder") << "PGeometricDetExtra and PGeometricDetExtraRcd Tag already present";
  }
}

void PGeometricDetExtraBuilder::putOne(const GeometricDetExtra& gde, PGeometricDetExtra* pgde) {
  PGeometricDetExtra::Item item;
  item._geographicalId = gde.geographicalId();
  item._volume = gde.volume();
  item._density = gde.density();
  item._weight = gde.weight();
  item._copy = gde.copyno();
  item._material = gde.material();
  pgde->pgdes_.push_back(item);
}

DEFINE_FWK_MODULE(PGeometricDetExtraBuilder);
