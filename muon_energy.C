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

enum LIMITS {
    MAX_TRACKS = 30000,
};

bool is_in_fv(const TLorentzVector &pos, const float margin = 20)
{
    if (pos.X() > 3 + margin && pos.X() < 253 - margin && pos.Y() > -113 + margin && pos.Y() < 113 - margin &&
        pos.Z() > 3 + margin && pos.Z() < 1034 - margin) {
        return true;
    }

    return false;
}

bool is_in_fv(const float x, const float y, const float z, const float margin = 20)
{
    TLorentzVector pos(x, y, z, 0);
    return is_in_fv(pos, margin);
}

bool is_good_proton(const TLorentzVector &start, const TLorentzVector &end, const float proton_daugher_ke,
                    const float max_proton_daugher_ke = 0)
{
    if (is_in_fv(start, 0) && is_in_fv(end, 0) && proton_daugher_ke <= max_proton_daugher_ke) {
        return true;
    }
    return false;
}

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

bool particle_match(const TLorentzVector &pa, const TLorentzVector &pb, const TLorentzVector &ma, const TLorentzVector &mb)
{
    // pos + direction
    // if (fabs(pa.X() - pb.X() + 0.1) > 0.5) return false;
    // if (fabs(pa.Y() - pb.Y()) > 1.) return false;
    // if (fabs(pa.Z() - pb.Z()) > 1.) return false;
    // if (fabs(ma.Theta() - mb.Theta()) > 0.5) return false;
    // if (fabs(ma.Phi() - mb.Phi()) > 0.5) return false;

    // dist
    auto dist = (pa - pb).Vect();
    if (dist.Mag() > 1.) return false;

    // if (fabs(ma.Theta() - mb.Theta()) > 1.0) return false;
    // if (fabs(ma.Phi() - mb.Phi()) > 1.0) return false;

    return true;
}

void muon_energy(
    // const char *input = "data/nue_overlay_run2.root"
    const char *input = "checkout_prodgenie_bnb_nu_overlay_run1_PF.root")
{
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

    auto *tf = TFile::Open(input, "read");
    auto *dir = (TDirectoryFile *) tf->Get("wcpselection");

    // old
    // auto *T_PFDump = (TTree *) dir->Get("T_PFDump");
    // auto *T_PFeval = (TTree *) dir->Get("T_PFeval");
    // T_PFDump->AddFriend(T_PFeval);

    // new
    auto *T_PFDump = (TTree *) dir->Get("T_PFeval");

    auto *T_BDTvars = (TTree *) dir->Get("T_BDTvars");
    T_PFDump->AddFriend(T_BDTvars);
    auto *T_eval = (TTree *) dir->Get("T_eval");
    T_PFDump->AddFriend(T_eval);
    auto *T_KINEvars = (TTree *) dir->Get("T_KINEvars");
    T_PFDump->AddFriend(T_KINEvars);

    const int target_pdg = 13;  // 11, 13, 2212

    int truth_Ntrack;
    int truth_id[MAX_TRACKS];
    int truth_pdg[MAX_TRACKS];
    int truth_process[MAX_TRACKS];
    int truth_mother[MAX_TRACKS];
    float truth_startXYZT[MAX_TRACKS][4];
    float truth_endXYZT[MAX_TRACKS][4];
    float truth_startMomentum[MAX_TRACKS][4];
    float truth_endMomentum[MAX_TRACKS][4];
    std::vector<std::vector<int> > *truth_daughters = 0;

    int reco_Ntrack;
    int reco_id[MAX_TRACKS];
    int reco_pdg[MAX_TRACKS];
    int reco_process[MAX_TRACKS];
    int reco_mother[MAX_TRACKS];
    float reco_startXYZT[MAX_TRACKS][4];
    float reco_endXYZT[MAX_TRACKS][4];
    float reco_startMomentum[MAX_TRACKS][4];
    float reco_endMomentum[MAX_TRACKS][4];
    std::vector<std::vector<int> > *reco_daughters = 0;

    bool truth_isCC;
    float truth_nu_momentum[4];
    float truth_corr_nuvtxX;
    float truth_corr_nuvtxY;
    float truth_corr_nuvtxZ;
    float reco_nuvtxX;
    float reco_nuvtxY;
    float reco_nuvtxZ;

    float numu_cc_flag;

    float stm_clusterlength;
    bool truth_isFC;

    float kine_reco_Enu;
    std::vector<float> *kine_energy_particle = 0;
    std::vector<int> *kine_energy_info = 0;
    std::vector<int> *kine_particle_type = 0;

    T_PFDump->SetBranchAddress("truth_Ntrack", &truth_Ntrack);
    T_PFDump->SetBranchAddress("truth_id", &truth_id);
    T_PFDump->SetBranchAddress("truth_pdg", &truth_pdg);
    T_PFDump->SetBranchAddress("truth_mother", &truth_mother);
    T_PFDump->SetBranchAddress("truth_startXYZT", &truth_startXYZT);
    T_PFDump->SetBranchAddress("truth_endXYZT", &truth_endXYZT);
    T_PFDump->SetBranchAddress("truth_startMomentum", &truth_startMomentum);
    T_PFDump->SetBranchAddress("truth_daughters", &truth_daughters);

    T_PFDump->SetBranchAddress("reco_Ntrack", &reco_Ntrack);
    T_PFDump->SetBranchAddress("reco_pdg", &reco_pdg);
    T_PFDump->SetBranchAddress("reco_mother", &reco_mother);
    T_PFDump->SetBranchAddress("reco_startXYZT", &reco_startXYZT);
    T_PFDump->SetBranchAddress("reco_endXYZT", &reco_endXYZT);
    T_PFDump->SetBranchAddress("reco_startMomentum", &reco_startMomentum);
    T_PFDump->SetBranchAddress("reco_daughters", &reco_daughters);

    T_PFDump->SetBranchAddress("truth_isCC", &truth_isCC);
    T_PFDump->SetBranchAddress("truth_nu_momentum", &truth_nu_momentum);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxX", &truth_corr_nuvtxX);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxY", &truth_corr_nuvtxY);
    T_PFDump->SetBranchAddress("truth_corr_nuvtxZ", &truth_corr_nuvtxZ);
    T_PFDump->SetBranchAddress("reco_nuvtxX", &reco_nuvtxX);
    T_PFDump->SetBranchAddress("reco_nuvtxY", &reco_nuvtxY);
    T_PFDump->SetBranchAddress("reco_nuvtxZ", &reco_nuvtxZ);

    T_PFDump->SetBranchAddress("numu_cc_flag", &numu_cc_flag);  // T_BDTvars

    T_PFDump->SetBranchAddress("stm_clusterlength", &stm_clusterlength);  // T_eval
    T_PFDump->SetBranchAddress("truth_isFC", &truth_isFC);                // T_eval

    T_PFDump->SetBranchAddress("kine_energy_particle", &kine_energy_particle);
    T_PFDump->SetBranchAddress("kine_energy_info", &kine_energy_info);
    T_PFDump->SetBranchAddress("kine_particle_type", &kine_particle_type);
    T_PFDump->SetBranchAddress("kine_reco_Enu", &kine_reco_Enu);

    TH1F *h_ndaughter_truth = new TH1F("h_ndaughter_truth", "h_ndaughter_truth", 100, -0.5, 99.5);
    TH1F *h_ndaughter_reco = new TH1F("h_ndaughter_reco", "h_ndaughter_reco", 100, -0.5, 99.5);

    TH1F *h_reco_m_truth = new TH1F("h_reco_m_truth", "h_reco_m_truth", 100, -1, 1);

    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 100, 0, 2, 100, 0, 2); // energy
    TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 20, 0, 2, 40, 0, 2); // energy bias
    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 20, 0, 3, 100, -TMath::Pi(), TMath::Pi()); //
    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 10, 0, 3, 100, -3, 3);  // angle
    // TH2F *h_reco_v_truth = new TH2F("h_reco_v_truth", "h_reco_v_truth", 20, 0, 3, 100, -1, 1);  // tmp

    TH1F *h_truth_e_all = new TH1F("h_truth_e_all", "h_truth_e_all", 60, 0, 3);
    TH1F *h_truth_e_match = new TH1F("h_truth_e_match", "h_truth_e_match", 60, 0, 3);

    TH2F *h_dtheta_theta =
        new TH2F("h_dtheta_theta", "h_dtheta_theta", 200, -TMath::Pi(), TMath::Pi(), 200, -TMath::Pi(), TMath::Pi());

    int counter_all = 0;
    int counter_pass = 0;
    for (int ientry = 0; ientry < T_PFDump->GetEntries(); ++ientry) {
    // for (int ientry = 0; ientry < 10000; ++ientry) {
        T_PFDump->GetEntry(ientry);
        if (ientry % 1000 == 0) cout << "processing: " << ientry / 10000. * 100 << "%" << endl;
        if (numu_cc_flag < 0 || stm_clusterlength < 15) continue;  // generic nu selection
        if (!truth_isCC) continue;
        if (truth_isFC!=true) continue; // FV cut
        // if (!is_in_fv(truth_corr_nuvtxX, truth_corr_nuvtxY, truth_corr_nuvtxZ)) continue; // alternative FV

        TVector3 truth_nuvtx(truth_corr_nuvtxX, truth_corr_nuvtxY, truth_corr_nuvtxZ);
        TVector3 reco_nuvtx(reco_nuvtxX, reco_nuvtxY, reco_nuvtxZ);
        if ((reco_nuvtx - truth_nuvtx).Mag() > 1.0) continue;

        TLorentzVector target_pos_start_truth;
        TLorentzVector target_pos_end_truth;
        TLorentzVector target_mom_start_truth;
        int target_index_truth = -1;
        float current_max_energy_truth = 0;
        // find the leading primary target particle
        std::map<int, int> map_id_itruth;
        for (int itruth = 0; itruth < truth_Ntrack; ++itruth) {
            map_id_itruth[truth_id[itruth]] = itruth;
            h_ndaughter_truth->Fill(truth_daughters->at(itruth).size());
            if (truth_mother[itruth] != 0) continue;
            if (truth_pdg[itruth] != target_pdg) continue;
            if (truth_startMomentum[itruth][3] < current_max_energy_truth) continue;
            target_index_truth = itruth;
            current_max_energy_truth = truth_startMomentum[itruth][3];
            target_pos_start_truth.SetXYZT(truth_startXYZT[itruth][0], truth_startXYZT[itruth][1],
                                           truth_startXYZT[itruth][2], truth_startXYZT[itruth][3]);
            target_pos_end_truth.SetXYZT(truth_endXYZT[itruth][0], truth_endXYZT[itruth][1], truth_endXYZT[itruth][2],
                                         truth_endXYZT[itruth][3]);
            target_mom_start_truth.SetXYZT(truth_startMomentum[itruth][0], truth_startMomentum[itruth][1],
                                           truth_startMomentum[itruth][2], truth_startMomentum[itruth][3]);
            // cout << ientry << ", truth: " << target_mom_start_truth.E() << endl;
        }
        // if (!is_in_fv(target_pos_start_truth)) continue;

        // leading reco target
        TLorentzVector pos_reco;
        TLorentzVector mom_reco;
        float current_max_energy_reco = 0;
        for (int ireco = 0; ireco < reco_Ntrack; ++ireco) {
            h_ndaughter_reco->Fill(reco_daughters->at(ireco).size());
            if (reco_mother[ireco] != 0) continue;
            if (reco_pdg[ireco] != target_pdg) continue;
            if (reco_startMomentum[ireco][3] < current_max_energy_reco) continue;
            current_max_energy_reco = reco_startMomentum[ireco][3];
            pos_reco.SetXYZT(reco_startXYZT[ireco][0], reco_startXYZT[ireco][1], reco_startXYZT[ireco][2],
                             reco_startXYZT[ireco][3]);
            mom_reco.SetXYZT(reco_startMomentum[ireco][0], reco_startMomentum[ireco][1], reco_startMomentum[ireco][2],
                             reco_startMomentum[ireco][3]);
            // cout << ientry << ", reco: " << mom_reco.E() << ", " << endl;
        }

        if (current_max_energy_truth == 0) continue;
        ++counter_all;

        // h_truth_e_all->Fill(truth_nu_momentum[3]);
        h_truth_e_all->Fill(target_mom_start_truth.E());

        if (current_max_energy_reco == 0) continue;

        // h_reco_m_truth->Fill((reco_nuvtx - truth_nuvtx).Mag());
        if (!particle_match(pos_reco, target_pos_start_truth, mom_reco, target_mom_start_truth)) continue;
        // if (fabs((reco_nuvtx - truth_nuvtx).Mag()) > 1.0) continue;
        ++counter_pass;
        // h_truth_e_match->Fill(truth_nu_momentum[3]);
        h_truth_e_match->Fill(target_mom_start_truth.E());

        h_reco_m_truth->Fill(mom_reco.E() - target_mom_start_truth.E());
        // h_reco_m_truth->Fill(pos_reco.X() - target_pos_start_truth.X());
        // h_reco_v_truth->Fill(target_mom_start_truth.E(), mom_reco.E());

        // find ke_energy_info
        float reco_ke = mom_reco.E()-mom_reco.M();
        int selected_ke_energy_info = -1;
        for (size_t ikine = 0; ikine < kine_energy_particle->size(); ++ikine) {
            if(fabs(reco_ke*1000.-(*kine_energy_particle)[ikine])<1) {
                selected_ke_energy_info = (*kine_energy_info)[ikine];
                break;
            }
        }
        if(selected_ke_energy_info==1) {
            h_reco_v_truth->Fill(target_mom_start_truth.E()-target_mom_start_truth.M(),reco_ke/(target_mom_start_truth.E()-target_mom_start_truth.M()));
            // h_reco_v_truth->Fill(target_mom_start_truth.E()-target_mom_start_truth.M(),reco_ke);
            // h_reco_v_truth->Fill(truth_nu_momentum[3],(kine_reco_Enu/1000.)/truth_nu_momentum[3]);
        } else if (selected_ke_energy_info==0) {
            h_reco_v_truth->Fill(target_mom_start_truth.E()-target_mom_start_truth.M(),reco_ke/(target_mom_start_truth.E()-target_mom_start_truth.M()));
            // h_reco_v_truth->Fill(target_mom_start_truth.E()-target_mom_start_truth.M(),reco_ke);
            // h_reco_v_truth->Fill(truth_nu_momentum[3],(kine_reco_Enu/1000.)/truth_nu_momentum[3]);
        }
        // angle
        // h_reco_v_truth->Fill(target_mom_start_truth.E(),mom_reco.Theta()-target_mom_start_truth.Theta());
        // h_reco_v_truth->Fill(target_mom_start_truth.E(), mom_reco.Phi() - target_mom_start_truth.Phi());
        // pos
        // h_reco_v_truth->Fill(target_pos_start_truth.Y(), pos_reco.Y());

        // angle
        // auto dtheta = TMath::ACos(target_mom_start_truth.Vect().Dot(mom_reco.Vect()) /
        //                           target_mom_start_truth.Vect().Mag() / mom_reco.Vect().Mag());
        // TLorentzVector target_mom_start_truth_yzx;
        // dtheta = mom_reco.Phi() - target_mom_start_truth.Phi();
        // target_mom_start_truth_yzx.SetXYZT(target_mom_start_truth.Z(), target_mom_start_truth.Y(),
        //                                    target_mom_start_truth.X(), target_mom_start_truth.T());
        // h_dtheta_theta->Fill(TMath::Pi() / 2 - target_mom_start_truth_yzx.Theta(),
        //                      TMath::Pi() / 2 - target_mom_start_truth.Theta());
    }
    cout << "target: truth " << counter_all << "; reco " << counter_pass
         << "; ratio: " << 100. * counter_pass / counter_all << "%" << endl;

    const int LINE_WIDTH = 2;

    TCanvas *c0 = new TCanvas("c0", "c0");
    h_truth_e_all->SetLineColor(kBlack);
    h_truth_e_all->SetLineWidth(LINE_WIDTH);
    h_truth_e_all->SetTitle(";E^{truth} [GeV]");
    h_truth_e_all->Draw();
    h_truth_e_match->SetLineColor(kRed);
    h_truth_e_match->SetLineWidth(LINE_WIDTH);
    h_truth_e_match->Draw("same");

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->SetGrid();
    auto pEff = new TEfficiency(*h_truth_e_match, *h_truth_e_all);
    // pEff->GetXaxis()->SetRangeUser(0,1);
    pEff->SetMarkerStyle(20);
    pEff->SetTitle(";E^{truth} [GeV];Efficiency");
    pEff->Draw("ap");
    // auto h_truth_e_ratio = (TH1F *) h_truth_e_match->Clone("h_truth_e_ratio");
    // h_truth_e_ratio->Divide(h_truth_e_all);
    // h_truth_e_ratio->SetLineColor(kRed);
    // h_truth_e_ratio->SetLineWidth(LINE_WIDTH);
    // h_truth_e_ratio->SetTitle(";E^{truth} [GeV];Efficiency");
    // h_truth_e_ratio->SetStats(0);
    // h_truth_e_ratio->Draw("");

    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->SetGrid();
    c2->SetLogz();
    // h_reco_v_truth->SetTitle(";Z^{truth} [cm];Z^{reco}-Z^{truth} [cm]");
    // h_reco_v_truth->SetTitle(";#phi^{truth} [rad];#phi^{reco}-#phi^{truth} [rad]");
    h_reco_v_truth->SetTitle(";E^{truth} [GeV];reco/true");
    // h_reco_v_truth->SetTitle(";E^{truth} [GeV];E^{reco} [GeV]");
    // h_reco_v_truth->SetTitle(";E^{truth} [GeV];#Delta #theta");
    // h_reco_v_truth->SetTitle(";E^{truth} [GeV];#Delta #phi");
    // h_reco_v_truth->SetStats(0);
    h_reco_v_truth->SetMinimum(1);
    h_reco_v_truth->Draw("colz");

    // resolution_hist
    auto h_reco_v_truth_1 = resolution_hist(h_reco_v_truth);
    h_reco_v_truth_1->Draw("e,same");

    TCanvas *c3 = new TCanvas("c3", "c3");
    c3->SetGrid();
    // h_reco_m_truth->SetTitle(";|Vtx^{reco}-Vtx^{truth}| [cm]");
    // h_reco_m_truth->SetTitle(";X^{reco}-X^{truth} [cm]");
    // h_reco_m_truth->SetTitle(";#phi^{reco}-#phi^{truth} [rad]");
    h_reco_m_truth->SetTitle(";E^{reco}-E^{truth} [GeV]");
    // h_reco_m_truth->SetStats(0);
    h_reco_m_truth->Draw();

    // TCanvas *c4 = new TCanvas("c4", "c4");
    // h_ndaughter_truth->Draw();
    // TCanvas *c5 = new TCanvas("c5", "c5");
    // h_ndaughter_reco->Draw();

    TCanvas *c6 = new TCanvas("c6", "c6");
    c6->SetGrid();
    h_dtheta_theta->SetTitle("; #theta_{YZ}; #theta_{XY}");
    // h_reco_m_truth->SetStats(0);
    h_dtheta_theta->Draw("colz");
}
