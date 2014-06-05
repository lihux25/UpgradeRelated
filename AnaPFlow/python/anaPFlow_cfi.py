import FWCore.ParameterSet.Config as cms

histAndTree = cms.EDFilter('histAndTree',
   muonSrc             = cms.InputTag("myPatMuonsPFIDIso"),
   eleSrc              = cms.InputTag("patElectronsIDIso"),
   jetSrc              = cms.InputTag("myPatJetsPFchsPt30noOverlapLeptons"),
   
   forVetoMuonSrc      = cms.InputTag("myPatMuonsPFIDIso"),
   forVetoElectronSrc  = cms.InputTag("patElectronsIDIso"),
   forVetoIsoTrkSrc    = cms.InputTag("trackIsolation"),

   genJetSrc           = cms.InputTag("ak5GenJets"),
   genMETSrc           = cms.InputTag("genMetCalo"),
   genParticleSrc      = cms.InputTag("genParticles"),

   recoGenJetsDR       = cms.double(0.3),

   vtxSrc              = cms.InputTag("goodVertices"),
   evtWeightInput      = cms.InputTag("weightProducer:weight"),
   nVtxPUcut           = cms.uint32(4),

   metSrc              = cms.InputTag("patMETsPF"),
   type1metSrc         = cms.InputTag("patMETsPF"),
   mhtSrc              = cms.InputTag("mymhtPFchs"),
   mht_forSgnfSrc      = cms.InputTag("mymhtPFforSgnf"),
   htSrc               = cms.InputTag("myhtPFchs"),

   dPhis_CUT_vec_Src   = cms.InputTag("jetMHTDPhiForSkimsRA2:jetMHTDPhiVec"),
   nJets_CUT_Src       = cms.InputTag("nJetsForSkimsRA2:nJets"),

   nLeptonsSels        = cms.vint32(-1, -1),

   debug               = cms.bool(False),

   doFillTree          = cms.bool(True),

   externalBitToTree_Src = cms.InputTag(""),

   doSgnf              = cms.bool(True),

   doFillGenTopInfo    = cms.bool(False),
   gentopIdxVecSrc     = cms.InputTag("topDecaysFilter:topIdxVec"),
   gentopDausIdxVecSrc = cms.InputTag("topDecaysFilter:topDausIdxVec"),
   genWDausIdxVecSrc   = cms.InputTag("topDecaysFilter:WDausIdxVec"),
   gentauDausIdxVecSrc = cms.InputTag("topDecaysFilter:tauDausIdxVec"),
   decayTypeVecSrc     = cms.InputTag("topDecaysFilter:decayTypeVec"),

   hepTaggerJetSrc     = cms.InputTag("HEPTopTagJets"),

   fillGenInfo         = cms.bool(False),

   genDecayStrVecSrc   = cms.InputTag("printDecay:decayStr"),
   genDecayChainPartIdxVecSrc = cms.InputTag("printDecay:decayChainPartIdxVec"),

   smsModelFileNameStrSrc = cms.InputTag("smsModelFilter:fileNameStr"),
   smsModelStrSrc = cms.InputTag("smsModelFilter:smsModelStr"),
   smsModelMotherMassSrc = cms.InputTag("smsModelFilter:smsMotherMass"),
   smsModelDaughterMassSrc = cms.InputTag("smsModelFilter:smsDaughterMass"),

   genJetsInputTags = cms.VInputTag(),
)                               
