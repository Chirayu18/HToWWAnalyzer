// -*- C++ -*-
//
// Package:    Analyzer/GenLevelInfoAnalyzer
// Class:      GenLevelInfoAnalyzer
//
/**\class GenLevelInfoAnalyzer GenLevelInfoAnalyzer.cc Analyzer/GenLevelInfoAnalyzer/plugins/GenLevelInfoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Chirayu Gupta
//         Created:  Mon, 27 Jan 2025 19:43:30 GMT
//
//

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class GenLevelInfoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit GenLevelInfoAnalyzer(const edm::ParameterSet&);
  ~GenLevelInfoAnalyzer() override;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  int cch=0,cchg=0,cchgg=0,gchc=0,lclhc=0,gghcc=0,other=0,total=0;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenLevelInfoAnalyzer::GenLevelInfoAnalyzer(const edm::ParameterSet& iConfig){
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
#endif
  genParticlesToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  //now do what ever initialization is needed
}

GenLevelInfoAnalyzer::~GenLevelInfoAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
		std::cout<<"cc -> h : "<<float(cch)/total<<std::endl;
		std::cout<<"cc -> hg : "<<float(cchg)/total<<std::endl;
		std::cout<<"cc -> hgg : "<<float(cchgg)/total<<std::endl;
		std::cout<<"gc -> hc : "<<float(gchc)/total<<std::endl;
		std::cout<<"gg -> hcc : "<<float(gghcc)/total<<std::endl;
		std::cout<<"l(g)c -> l(g)hc : "<<float(lclhc)/total<<std::endl;
		std::cout<<"Unidentified Process : "<<float(other)/total<<std::endl;
		std::cout<<"Total events : "<<total<<std::endl;
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void GenLevelInfoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);
  std::vector<int> initialState;
  std::vector<int> finalState;

    for (const auto& particle : *genParticles) {
	if(particle.isHardProcess())
	{
		if(particle.pdgId() == 24 ) // Higgs decay encountered. No need to continue beyond this
			break;

		if(particle.status()==21)
			initialState.push_back(particle.pdgId());
		if(particle.status()>21)
			finalState.push_back(particle.pdgId());
		//std::cout << "PDG ID: " << particle.pdgId()
		//<< ", Status: " << particle.status()
		//<< ", Pt: " << particle.pt()
		//<< ", Eta: " << particle.eta()
		//<< ", Phi: " << particle.phi()
		//<< ", Mass: " << particle.mass()
		//<< ", Number of Daughters: " << particle.numberOfDaughters()
		//<< ", Number of Mothers: " << particle.numberOfMothers()
		//<< std::endl;
	}
  }


	//std::cout << "Initial State: ";
	//for (const int &num : initialState) {
	    //std::cout << num << " ";
	//}
	//std::cout << std::endl;
	//std::cout << "Final State: ";
	//for (const int &num : finalState) {
	    //std::cout << num << " ";
	//}
	//std::cout << std::endl;
	total++;
	if( initialState == std::vector<int>{4,-4} || initialState == std::vector<int>{-4,4})
	{
		if(finalState.size() == 1 )
		{
			//std::cout<<"cc -> h"<<std::endl;
			cch++;
		}
		else if(finalState.size() == 2)
		{
			//std::cout<<"cc -> hg"<<std::endl;
			cchg++;
		}
		else if(finalState.size() == 3)
		{
			//std::cout<<"cc -> hgg"<<std::endl;
			cchgg++;
		}
		else
		{
			//std::cout<<"Unidentified Process"<<std::endl;
			other++;
		}
	}
	else if( initialState == std::vector<int>{21,-4} || initialState == std::vector<int>{21,4}|| initialState == std::vector<int>{4,21}|| initialState == std::vector<int>{-4,21})
	{
		if( finalState == std::vector<int>{25,-4} || finalState == std::vector<int>{25,4}|| finalState == std::vector<int>{4,25}|| finalState == std::vector<int>{-4,25})
		{
			//std::cout<<"gc -> hc"<<std::endl;
			gchc++;
		}
		else
		{
			//std::cout<<"l(g)c -> l(g)hc"<<std::endl;
			lclhc++;
		}


	}
	else if( initialState == std::vector<int>{21,21})
	{
		//std::cout<<"gg -> hcc"<<std::endl;
		gghcc++;
	}
	else
	{
		//std::cout<<"l(g)c -> l(g)hc"<<std::endl;
		lclhc++;
	}


	//std::cout << " ======================================== "<<std::endl;
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void GenLevelInfoAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void GenLevelInfoAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenLevelInfoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenLevelInfoAnalyzer);
