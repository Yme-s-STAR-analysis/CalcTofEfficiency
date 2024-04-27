#ifndef _StTofMatchMaker_head
#define _StTofMatchMaker_head
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TString.h"
#include "TVector3.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoDstMaker;
class TH1F;
class TH2F;
class TProfile;
class TTree;
class TH2D;
class TEfficiency;

class StCFMult;
class TpcShiftTool;
class TriggerTool;
class MeanDcaTool;
class CentCorrTool;
class VtxShiftTool;


class StTofMatchMaker : public StMaker {
	public:
		StTofMatchMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName="tofMatchTree.root");
		virtual ~StTofMatchMaker();

		virtual Int_t Init();
		virtual Int_t Make();
		virtual void  Clear(Option_t *opt="");
		virtual Int_t Finish();

		void set_cut_dca(double value) {
			mCut_dca = value;
		}

		void set_cut_nHitsFit(int value) {
			mCut_nHitsFit = value;
		}
	  
		void set_cut_nSig(double value) {
			mCut_nSig = value;
		}

		Int_t vz_split(double vz);

	private:
		StPicoDstMaker *mPicoDstMaker;
		StPicoDst      *mPicoDst;
		StPicoEvent	   *event;
		StPicoTrack    *picoTrack;

		StCFMult* mtMult;
		TpcShiftTool* mtShift;
		CentCorrTool* mtCent;
		MeanDcaTool* mtDca;
		TriggerTool* mtTrg;
		VtxShiftTool* mtVtx;

		static const int nCent = 9;
		static const int nVz = 5;

		// quality cuts
		double mCut_dca;
		int mCut_nHitsFit;
        double mCut_nSig;
	
		TEfficiency* teff_y_pt[nCent][nVz][2];

		TString mOutputName;
		TFile* mFileOut;

		ClassDef(StTofMatchMaker, 1)
};

ClassImp(StTofMatchMaker)

#endif
