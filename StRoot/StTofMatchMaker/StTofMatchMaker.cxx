#include "StTofMatchMaker.h"

#include <TMath.h>

#include <algorithm>
#include <fstream>
#include <vector>

#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StThreeVectorF.hh"
#include "Stiostream.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "phys_constants.h"

#include "StRoot/CentCorrTool/CentCorrTool.h"
#include "StRoot/MeanDcaTool/MeanDcaTool.h"
#include "StRoot/TpcShiftTool/TpcShiftTool.h"
#include "StRoot/TriggerTool/TriggerTool.h"
#include "StRoot/StCFMult/StCFMult.h"

StTofMatchMaker::StTofMatchMaker(
	const char* name, 
	StPicoDstMaker* picoMaker,
    const char* outName
) : StMaker(name) {
	mOutputName = outName;
	mPicoDstMaker = picoMaker;
	mPicoDst = 0;
}

StTofMatchMaker::~StTofMatchMaker() {}

Int_t StTofMatchMaker::Init() {
  	mFileOut = new TFile(mOutputName, "recreate");
	const char* pname[2] = {"Pro", "Pbar"};
	for (int ptype=0; ptype<2; ptype++) {
		for (int iCent=0; iCent<nCent; iCent++) {
			for (int iVz=0; iVz<nVz; iVz++) {
				teff_y_pt[iCent][iVz][ptype] = new TEfficiency(
					Form("TofEff_cent%d_vz%d_%s", iCent, iVz, pname[ptype]),
					";y;p_{T} [GeV/c];Efficiency", 30, -1.5, 1.5, 21, 0.0, 2.1
				);
			}
		}
	}
	// initialize costume modules

	// mean dca tool
	mtDca = new MeanDcaTool();
	mtDca->ReadParams();

	// centrality tool
	mtCent = new CentCorrTool();
	mtCent->EnableIndianMethod(true);
	mtCent->ReadParams();

	// multiplicity and shift tool
	mtShift = new TpcShiftTool();
	mtShift->Init();
	mtMult = new StCFMult();
	mtMult->ImportShiftTool(mtShift);

	// trigger tool
	mtTrg = new TriggerTool();

	// set cuts
	mCut_dca = 1.0;
	mCut_nHitsFit = 20;
    mCut_nSig = 2.0;

	return kStOK;
}

//---------------------------------------------------------
Int_t StTofMatchMaker::vz_split(double vz) {
	if (-30 < vz && vz < -10) {
		return 0;
	} else if (-10 < vz && vz < 10) {
		return 1;
	} else if (10 < vz && vz < 30) {
		return 2;
	} else if (-50 < vz && vz < -30) {
		return 3;
	} else if (30 < vz && vz < 50) {
		return 4;
	} else {
		return -1;
	}
}

//---------------------------------------------------------
Int_t StTofMatchMaker::Finish() {
	mFileOut->cd();
	for (int ptype=0; ptype<2; ptype++) {
		for (int iCent=0; iCent<nCent; iCent++) {
			for (int iVz=0; iVz<nVz; iVz++) {
				teff_y_pt[iCent][iVz][ptype]->Write();
			}
		}
	}
	mFileOut->Close();
	return kStOK;
}

void StTofMatchMaker::Clear(Option_t* opt) {}

//---------------------------------------------------------------
Int_t StTofMatchMaker::Make() {
	if (!mPicoDstMaker) {
		LOG_WARN << " No PicoDstMaker! Skip! " << endm;
		return kStWarn;
	}

	mPicoDst = mPicoDstMaker->picoDst();
	if (!mPicoDst) {
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStWarn;
	}

	if (!mPicoDst) {
		return kStOK;
	}

	// Load event
	event = (StPicoEvent*)mPicoDst->event();
	if (!event) {
		cerr << "Error opening picoDst Event, skip!" << endl;
		return kStOK;
	}

	TVector3 pVtx = event->primaryVertex();
	Double_t vx = pVtx.X();
	Double_t vy = pVtx.Y();
	Double_t vz = pVtx.Z();

	if (fabs(vx) < 1.e-5 && 
		fabs(vy) < 1.e-5 &&
		fabs(vz) < 1.e-5) {
		return kStOK;
	}

	// using Ashish's shifted vr cut
	// -> see: https://drupal.star.bnl.gov/STAR/system/files/Vr_xy_N_Vzcut.txt
	vx = vx - 0.0417;
	vy = vy + 0.2715;
	Double_t vr = sqrt(vx * vx + vy * vy);

	if (vr >= 1.0 || fabs(vz) > 50.0) {
		return kStOK;
	}

	Int_t runId = event->runId();
	Int_t trgid = mtTrg->GetTriggerID(event);
	if (trgid < 0) { return kStOK; }

	mtMult->make(mPicoDst);
	Int_t refMult = mtMult->mRefMult;
	Int_t tofMult = mtMult->mTofMult;
	Int_t nTofMatch = mtMult->mNTofMatch;
	Int_t nTofBeta = mtMult->mNTofBeta;

	Int_t refMult3 = mtMult->mRefMult3;
	refMult3 = mtCent->GetIndianRefMult3Corr(
		refMult, refMult3, tofMult, nTofMatch, nTofBeta,
		vz, false
	);
	if (refMult3 < 0) { return kStOK; }
	Int_t cent = mtCent->GetCentrality9(refMult3);
	if (cent < 0 || cent >= 9) { return kStOK; }

	// check DCA
	if (!mtDca->Make(mPicoDst)) { return kStOK; }
	if (mtDca->IsBadMeanDcaZEvent(mPicoDst) || mtDca->IsBadMeanDcaXYEvent(mPicoDst)) {
		return kStOK;
	}

	// track loop
  	Int_t nTracks = mPicoDst->numberOfTracks();
	const Float_t mField = event->bField();

	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
		picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
		if (!picoTrack) { continue; }

		if (!picoTrack->isPrimary()) { continue; }

		// Float_t dca = fabs(picoTrack->gDCA(vx, vy, vz));
		StPicoPhysicalHelix helix = picoTrack->helix(mField);
        Float_t dca = fabs(helix.geometricSignedDistance(pVtx));

		TVector3 pmomentum = picoTrack->pMom();
		Float_t p = pmomentum.Mag();
		Float_t pt = pmomentum.Perp();
		Float_t pz = pmomentum.Z();
		Float_t eta = pmomentum.PseudoRapidity();

		Float_t EP = sqrt(p * p + 0.938272 * 0.938272);
		Float_t YP = 0.5 * log((EP + pz) / (EP - pz));
		if (isnan(YP)) { continue; }

    	Int_t nHitsFit = picoTrack->nHitsFit();
    	Int_t nHitsdEdx = picoTrack->nHitsDedx();
    	Int_t nHitsPoss = picoTrack->nHitsMax();
    	Float_t nHitsRatio = nHitsFit*1.0 / nHitsPoss;
    	Float_t nSigProton = picoTrack->nSigmaProton();
    	Int_t charge = (Int_t)picoTrack->charge();

		Int_t btofMatchFlag = 0;
		Int_t tofId = picoTrack->bTofPidTraitsIndex();
        Double_t beta = -1.0;
        Double_t btofYLocal = -999.0;
        if (tofId >= 0) {
            StPicoBTofPidTraits* tofPid = mPicoDst->btofPidTraits(tofId);
            btofMatchFlag = tofPid->btofMatchFlag();
            if (tofPid) {
                beta = tofPid->btofBeta();
                btofYLocal = tofPid->btofYLocal();
                if (beta < 1e-4) { // recalculate time of flight
                    Double_t tof = tofPid->btof();
                    TVector3 btofHitPos = tofPid->btofHitPos();
                    const StThreeVectorF* btofHitsPosSt = new StThreeVectorF(
                        btofHitPos.X(),btofHitPos.Y(),btofHitPos.Z()
                    );
                    const StThreeVectorF* vtxPosSt = new StThreeVectorF(
                        vx, vy, vz
                    );
                    Double_t L = tofPathLength(vtxPosSt, btofHitsPosSt, helix.curvature());
                    beta = tof > 0 ? L / (tof * (C_C_LIGHT/1.e9)) : std::numeric_limits<Float_t>::quiet_NaN(); // note: quiet nan will never pass > N or < N
                }
            }
        }

		if (nHitsdEdx < 5 || nHitsFit < mCut_nHitsFit || nHitsRatio < 0.52) { // SYS ERR FLAG
			continue;
		}
		if (dca > mCut_dca) { // SYS ERR FLAG
			continue;
		}

		// ignore TPC PID here
        nSigProton -= mtShift->GetShift(runId, pt, eta);
		if (fabs(nSigProton) > mCut_nSig) { // SYS ERR FLAG
			continue;
		}

    	bool acc = btofMatchFlag > 0 && beta > 0 && fabs(btofYLocal) < 1.8;
    	Int_t pidx = (Int_t)(charge < 0);  // 0 for proton, 1 for antiproton
		// fill efficiency as function of eta and y
		Int_t vzBin = vz_split(vz);
		if (vzBin < 0) { continue; }
		teff_y_pt[cent][vzBin][pidx]->Fill(acc, YP, pt);

	}  // picotracks loop end
	return kStOK;
}
