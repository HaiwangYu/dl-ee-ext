#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"

#include <vector>
#include <iostream>

using namespace std;


TH1 *resolution_hist(TH2 *in)
{
    TString basename(in->GetName());
    auto current_dir = in->GetDirectory();
    in->FitSlicesY();
    auto mean = (TH1 *) current_dir->Get(basename + "_1");
    auto sigma = (TH1 *) current_dir->Get(basename + "_2");
    auto out = (TH1 *) mean->Clone(basename + "res");
    for (int ix = 1; ix < out->GetNbinsX(); ++ix) {
        out->SetBinError(ix, sigma->GetBinContent(ix));
    }
    return out;
}

float rms(TH1 *h)
{
    float sum = 0;
    float nentries_valid = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        sum += h->GetBinCenter(i) * h->GetBinCenter(i) * h->GetBinContent(i);
        nentries_valid += h->GetBinContent(i);
    }
    return sqrt(sum / nentries_valid);
}

float rms(const std::vector<float> &v)
{
    float sum = 0;
    for (auto x : v) {
        // if (fabs(x) > 1) continue;
        sum += x * x;
    }
    return sqrt(sum / v.size());
}

void add_corr(
    const char *input = "sel.root",
    const char *output = "corr.root") {

    auto *fin = TFile::Open(input, "read");
    auto *T_sel = (TTree *) fin->Get("T_sel");
    float E_l_rec, E_l_rec_corr, E_nu_rec, E_nu_rec_corr, E_l_tru, E_nu_tru;
    int ke_info;
    float vtx_diff;
    bool truth_isFC;
    bool match_isFC;
    T_sel->SetBranchAddress("E_l_rec", &E_l_rec);
    T_sel->SetBranchAddress("E_nu_rec", &E_nu_rec);
    T_sel->SetBranchAddress("E_l_tru", &E_l_tru);
    T_sel->SetBranchAddress("E_nu_tru", &E_nu_tru);
    T_sel->SetBranchAddress("ke_info", &ke_info);
    T_sel->SetBranchAddress("vtx_diff", &vtx_diff);
    T_sel->SetBranchAddress("truth_isFC", &truth_isFC);
    T_sel->SetBranchAddress("match_isFC", &match_isFC);
    
    TFile *fout = new TFile(output, "recreate");
    TTree *T_corr = new TTree("T_corr", "T_corr");
    T_corr->Branch("E_l_rec", &E_l_rec, "E_l_rec/F");
    T_corr->Branch("E_l_rec_cor", &E_l_rec_corr, "E_l_rec_corr/F");
    T_corr->Branch("E_nu_rec", &E_nu_rec, "E_nu_rec/F");
    T_corr->Branch("E_nu_rec_cor", &E_nu_rec_corr, "E_nu_rec_corr/F");
    T_corr->Branch("E_l_tru", &E_l_tru, "E_l_tru/F");
    T_corr->Branch("E_nu_tru", &E_nu_tru, "E_nu_tru/F");
    T_corr->Branch("ke_info", &ke_info, "ke_info/I");
    T_corr->Branch("vtx_diff", &vtx_diff, "vtx_diff/F");
    T_corr->Branch("truth_isFC", &truth_isFC, "truth_isFC/O");
    T_corr->Branch("match_isFC", &match_isFC, "match_isFC/O");

    TH2F *h_nu_diff_vs_nu_rec = new TH2F("h_nu_diff_vs_nu_rec", "h_nu_diff_vs_nu_rec", 50, 0, 5, 100, -5, 5);
    int nentries = T_sel->GetEntries();
    // nentries = 10000;
    for (int ientry = 0; ientry < nentries; ++ientry) {
        T_sel->GetEntry(ientry);
        if (match_isFC == false) continue;
        h_nu_diff_vs_nu_rec->Fill(E_nu_rec, E_nu_tru - E_nu_rec);
    }

    const int LINE_WIDTH = 2;

    // TCanvas *c2 = new TCanvas("c2", "c2");
    // c2->SetGrid();
    // c2->SetLogz();
    h_nu_diff_vs_nu_rec->SetTitle(";E_{\nu}^{rec} [GeV];true-reco");
    h_nu_diff_vs_nu_rec->SetMinimum(1);
    // h_nu_diff_vs_nu_rec->Draw("colz");

    // resolution_hist
    auto h_nu_diff_vs_nu_rec_res = resolution_hist(h_nu_diff_vs_nu_rec);
    // h_nu_diff_vs_nu_rec_res->Draw("e,same");

    std::vector<float> rms_rec;
    std::vector<float> rms_corr;
    for (int ientry = 0; ientry < nentries; ++ientry) {
        T_sel->GetEntry(ientry);
        if (match_isFC == false) continue;
        const int ibin = h_nu_diff_vs_nu_rec_res->FindBin(E_nu_rec);
        const float corr = h_nu_diff_vs_nu_rec_res->GetBinContent(ibin);
        E_nu_rec_corr = E_nu_rec + corr;
        T_corr->Fill();
        rms_rec.push_back((E_nu_tru - E_nu_rec) / E_nu_tru);
        rms_corr.push_back((E_nu_tru - E_nu_rec_corr) / E_nu_tru);
    }
    std::cout << "rms_rec: " << rms(rms_rec) << std::endl;
    std::cout << "rms_corr: " << rms(rms_corr) << std::endl;
    fout->cd();
    T_corr->Write();
    h_nu_diff_vs_nu_rec->Write();
    h_nu_diff_vs_nu_rec_res->Write();
    fout->Close();
}