#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;
class StTofMatchMaker;

StChain *chain;
void readPicoDst(const Char_t *inputFile = "file.list", TString JobIdName = "tof")
{
    Int_t nEvents = 15000000;

    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();

    gSystem->Load("StUtilities");
    gSystem->Load("StPicoEvent");
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StTofMatchMaker");
    gSystem->Load("StCFMult");
    gSystem->Load("TpcShiftTool");
    gSystem->Load("CentCorrTool");
    gSystem->Load("MeanDcaTool");
    gSystem->Load("TriggerTool");

    chain = new StChain();

    TString Name = JobIdName;
    Name.Append(".root");
    StPicoDstMaker *picoMaker = new StPicoDstMaker(2, inputFile, "picoDst");
    StTofMatchMaker *anaMaker = new StTofMatchMaker("ana", picoMaker, Name);

    anaMaker->set_cut_dca(1.0);
    anaMaker->set_cut_nHitsFit(20);
    anaMaker->set_cut_nSig(2.0);

    chain->Init();
    cout << "chain->Init();" << endl;
    int total = picoMaker->chain()->GetEntries();
    cout << " Total entries = " << total << endl;
    if (nEvents > total)
        nEvents = total;
    for (Int_t i = 0; i < nEvents; i++)
    {
        if (i % 10000 == 0)
            cout << "Working on eventNumber " << i << endl;

        chain->Clear();
        int iret = chain->Make(i);

        if (iret)
        {
            cout << "Bad return code!" << iret << endl;
            break;
        }

        total++;
    }

    cout << "****************************************** " << endl;
    cout << "Work done... now its time to close up shop!" << endl;
    cout << "****************************************** " << endl;
    chain->Finish();
    cout << "****************************************** " << endl;
    cout << "total number of events  " << nEvents << endl;
    cout << "****************************************** " << endl;

    delete chain;
}
