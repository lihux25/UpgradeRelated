#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Common/interface/Ptr.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"

//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

typedef unsigned int size;
static const int sgnfnDof = 2;

using namespace edm;
using namespace std;

//For sorting by pt
struct GreaterByPtCandPtr
{
  bool operator()( const edm::Ptr<reco::Candidate> & t1, const edm::Ptr<reco::Candidate> & t2 ) const
  {return t1->pt() > t2->pt();}
};

class anaPFlow : public edm::EDFilter{

  public:

    explicit anaPFlow(const edm::ParameterSet & iConfig);
    ~anaPFlow();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob();
    virtual bool beginRun(edm::Run&, const edm::EventSetup&);
    virtual bool endRun(edm::Run&, const edm::EventSetup&);

    size run, event, ls; bool isdata;
    edm::InputTag vtxSrc_;
    edm::Handle<edm::View<reco::Vertex> > vertices;
    size nVtxPUcut_, vtxSize;
    void loadEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);

    edm::InputTag jetSrc_;
    edm::Handle<edm::View<reco::Jet> > jets;
    virtual void loadRecoJets(const edm::Event& iEvent);

    edm::InputTag muonSrc_;
    edm::Handle<edm::View<reco::Muon > > muons;
    edm::InputTag eleSrc_;
    edm::Handle<edm::View<reco::GsfElectron> > electrons;
    size nMuons, nElectrons;
    int muon1Charge, muon2Charge;
    int ele1Charge, ele2Charge;
    virtual void loadLeptons(const edm::Event& iEvent);

    edm::InputTag metSrc_;
    edm::Handle<edm::View<reco::MET> > met;
    virtual void loadMETMHT(const edm::Event& iEvent);

    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    int npv; double avg_npv; int nm1, n0, np1; double tru_npv;
    virtual void loadGenInfo(const edm::Event& iEvent);

    bool debug_;

    bool doFillTree_;

    bool isData;

    bool printOnce_;

    TTree *outTree;
    std::vector<TLorentzVector> *jetsLVec_TR;

    std::vector<TLorentzVector> *muonsLVec_TR;
    std::vector<double> *muonsAux_TR;
    std::vector<TLorentzVector> *elesLVec_TR;
    std::vector<double> *elesAux_TR;

    void setTreeDefaultVars();
};

anaPFlow::anaPFlow(const edm::ParameterSet & iConfig) {

   isData = true;
 
   jetSrc_ = iConfig.getParameter<edm::InputTag>("jetSrc");
   muonSrc_ = iConfig.getParameter<edm::InputTag>("muonSrc");
   eleSrc_ = iConfig.getParameter<edm::InputTag>("eleSrc");

   vtxSrc_ = iConfig.getParameter<edm::InputTag>("vtxSrc");

   metSrc_ = iConfig.getParameter<edm::InputTag>("metSrc");
 
   debug_ = iConfig.getParameter<bool>("debug");

   doFillTree_ = iConfig.getParameter<bool>("doFillTree");

   npv = -1; avg_npv = -1; nm1 = -1; n0 = -1; np1 = -1; tru_npv = -1;

   jetsLVec_TR = new std::vector<TLorentzVector>;
   muonsLVec_TR = new std::vector<TLorentzVector>; muonsAux_TR = new std::vector<double>;
   elesLVec_TR = new std::vector<TLorentzVector>; elesAux_TR = new std::vector<double>;

   setTreeDefaultVars();

   edm::Service<TFileService> fs;

   if( doFillTree_ ){

      outTree = fs->make<TTree>("AUX", "aux info");
      outTree->Branch("run", &run, "run/I");
      outTree->Branch("event", &event, "event/I");
      outTree->Branch("lumi", &ls, "lumi/I");
      outTree->Branch("npv", &npv, "npv/I");
      outTree->Branch("avg_npv", &avg_npv, "avg_npv/D");
      outTree->Branch("nm1", &nm1, "nm1/I");
      outTree->Branch("n0", &n0, "n0/I");
      outTree->Branch("np1", &np1, "np1/I");
      outTree->Branch("tru_npv", &tru_npv, "tru_npv/D");
      outTree->Branch("vtxSize", &vtxSize, "vtxSize/I");
      outTree->Branch("nJets", &nJets, "nJets/I");
      outTree->Branch("nMuons", &nMuons, "nMuons/I");
      outTree->Branch("nElectrons", &nElectrons, "nElectrons/I");
      outTree->Branch("met", &met_TR, "met/D");
      outTree->Branch("metphi", &metphi_TR, "metphi/D");
      outTree->Branch("jetsLVec", "std::vector<TLorentzVector>", &jetsLVec_TR, 32000, 0);
      outTree->Branch("muonsLVec", "std::vector<TLorentzVector>", &muonsLVec_TR, 32000, 0);
      outTree->Branch("muonsAux", "std::vector<double>", &muonsAux_TR, 32000, 0);
      outTree->Branch("elesLVec", "std::vector<TLorentzVector>", &elesLVec_TR, 32000, 0);
      outTree->Branch("elesAux", "std::vector<double>", &elesAux_TR, 32000, 0);
   }
}

anaPFlow::~anaPFlow() {

}

void anaPFlow::setTreeDefaultVars(){

   mt_TR= -99; metphi_TR= -99;
   jetsLVec_TR->clear();

   muonsLVec_TR->clear(); elesLVec_TR->clear();
   muonsAux_TR->clear(); elesAux_TR->clear();

   npv = -1; avg_npv = -1; nm1 = -1; n0 = -1; np1 = -1; tru_npv = -1;
}

// ------------ method called on each new Event  ------------
bool anaPFlow::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

   setTreeDefaultVars();

//   iSetup.getData( pdt_ );

// Event setup
   loadEventInfo(iEvent, iSetup);

   loadGenInfo(iEvent);
   loadRecoJets(iEvent);
   loadLeptons(iEvent);
   loadMETMHT(iEvent);

   met_TR = (*met)[0].pt(); metphi_TR = (*met)[0].phi();

   for(size ij=0; ij<nJets; ij++){
      TLorentzVector perJetLVec;
      perJetLVec.SetPtEtaPhiE( (*jets)[ij].pt(), (*jets)[ij].eta(), (*jets)[ij].phi(), (*jets)[ij].energy() );
      jetsLVec_TR->push_back(perJetLVec);
   }

   for(size im=0; im<nMuons; im++){
      TLorentzVector perMuonLVec;
      perMuonLVec.SetPtEtaPhiE( (*muons)[im].pt(), (*muons)[im].eta(), (*muons)[im].phi(), (*muons)[im].energy() );
      muonsLVec_TR->push_back(perMuonLVec);
      muonsAux_TR->push_back( (*muons)[im].charge() );
   }

   for(size ie=0; ie<nElectrons; ie++){
      TLorentzVector perEleLVec;
      perEleLVec.SetPtEtaPhiE( (*electrons)[ie].pt(), (*electrons)[ie].eta(), (*electrons)[ie].phi(), (*electrons)[ie].energy() );
      elesLVec_TR->push_back(perEleLVec);
      elesAux_TR->push_back( (*electrons)[ie].charge() );
   }

   if( doFillTree_ ){
      outTree->Fill(); 
   }

   return true;

}

// ------------ method called once each job just before starting event loop  ------------
void anaPFlow::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void anaPFlow::endJob() {
}

// ------------ method called once each run just before starting event loop  ------------
bool anaPFlow::beginRun(edm::Run &run, const edm::EventSetup& iSetup) {
  return true;
}

// ------------ method called once each run just after starting event loop  ------------
bool anaPFlow::endRun(edm::Run &run, const edm::EventSetup& iSetup) {
  return true;
}


void anaPFlow::loadEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup){

// Determine if it's data
   if( !iEvent.isRealData() ) isData = false;

// Get run, event, lumi info
   run = iEvent.id().run();
   event = iEvent.id().event();
   ls = iEvent.luminosityBlock();

// Get vertices
   iEvent.getByLabel(vtxSrc_, vertices); vtxSize = vertices->size();
   
// Get event weight
   iEvent.getByLabel(evtWeightInput_, evtWeight_);

}

void anaPFlow::loadGenInfo(const edm::Event& iEvent){

// MC generate level related info
   scalePDF = -1; pthat = -1;
   if (!isData) {
      iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator PVI;

      npv = -1; avg_npv = 0; nm1 = -1; n0 = -1; np1 = -1; tru_npv = -1;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

         int BX = PVI->getBunchCrossing();

         avg_npv += double(PVI->getPU_NumInteractions());

         if(BX == -1) { 
            nm1 = PVI->getPU_NumInteractions();
         }
         if(BX == 0) { 
            n0 = PVI->getPU_NumInteractions();
         }
         if(BX == 1) { 
            np1 = PVI->getPU_NumInteractions();
         }

         if(BX == 0) {
            npv = PVI->getPU_NumInteractions();
            tru_npv = PVI->getTrueNumInteractions();
            continue; // No touching of this "continue", since I'm not sure why it's here originally
         }
      }
      avg_npv /= 3.0;
   }
}

void anaPFlow::loadRecoJets(const edm::Event& iEvent){
   iEvent.getByLabel(jetSrc_, jets); nJets = jets->size();
}

void anaPFlow::loadLeptons(const edm::Event& iEvent){
   iEvent.getByLabel(muonSrc_, muons); nMuons = muons->size();
   iEvent.getByLabel(eleSrc_, electrons); nElectrons = electrons->size();
}

void anaPFlow::loadPhotons(const edm::Event& iEvent){
}

void anaPFlow::loadMETMHT(const edm::Event& iEvent){
   iEvent.getByLabel(metSrc_, met);
}

void anaPFlow::loadHT(const edm::Event& iEvent){
}

void anaPFlow::loadAUX(const edm::Event& iEvent){
}

void anaPFlow::fillGenDecayInfo(const edm::Event& iEvent){
   if( !isData ){
   }
}

int anaPFlow::getConsMatchedJetIdx(const std::vector<pat::Jet> & patjets, const TLorentzVector tomatchLVec, const double tomatchCharge, const double minDRcut){

   int consMatchedJetsIdx = -1;
   double tomatchEta = tomatchLVec.Eta(), tomatchPhi = tomatchLVec.Phi();
   for(int ij=0; ij<(int)patjets.size(); ij++){

      const std::vector<reco::PFCandidatePtr> & constituents = patjets[ij].getPFConstituents();
      const unsigned int numConstituents = constituents.size();

      double minDRCon = 999.; int selConIdx = -1;
      for (unsigned int iCon = 0; iCon < numConstituents; ++iCon){
         const reco::PFCandidatePtr& constituent = constituents[iCon];
         if( constituent->charge() != tomatchCharge ) continue;
         const double tomatchDR = reco::deltaR(constituent->eta(), constituent->phi(), tomatchEta, tomatchPhi);
         if( minDRCon > tomatchDR ){
            minDRCon = tomatchDR;
            selConIdx = (int)iCon;
         }
      }
      if( selConIdx ){/*empty to avoid a compiling error*/}

      if( minDRCon < minDRcut ){
         consMatchedJetsIdx = ij;
         break;
      }
   }

   return consMatchedJetsIdx;
}

int anaPFlow::find_idx(const reco::Candidate & target){
   int pickedIdx = -1;
   for(size_t ig=0; ig<genParticles->size(); ig++){
      const reco::GenParticle& gen = genParticles->at(ig);
      if( target.p4() == gen.p4() && target.vertex() == gen.vertex() && target.charge() == gen.charge() ){
         pickedIdx = (int)ig;
         break;
      }
   }
   return pickedIdx;
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(anaPFlow);
