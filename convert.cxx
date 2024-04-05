void convert() {
    TFile* tfin = new TFile("TofEffRaw.root");
    TFile* tfout = new TFile("TofEff.root", "recreate");

    const int nCent = 9;
    const int nVz = 3;
    const int nType = 2;
    const char* pNames[nType] = {"Pro", "Pbar"};

    TEfficiency* teff[nCent][nVz][nType];
    TH2D* th2[nCent][nVz][nType];

    tfout->cd();
    for (int i=0; i<nCent; i++) {
        for (int j=0; j<nType; j++) {
            for (int k=0; k<nVz; k++) {
                teff[i][k][j] = (TEfficiency*)tfin->Get(Form("TofEff_cent%d_vz%d_%s", i, k, pNames[j]));
                th2[i][k][j] = (TH2D*)teff[i][k][j]->CreateHistogram();
                th2[i][k][j]->SetName(Form("TofEff_cent%d_vz%d_%s", i, k, pNames[j]));
                th2[i][k][j]->Write();
            }
        }
    }
    tfin->Close();
    tfout->Close();
}
