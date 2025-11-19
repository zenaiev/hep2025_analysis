// matching.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TDecompLU.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPad.h"
#include "Math/Vector4D.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldIds.h"

// g++ matching.cxx -o matching `root-config --libs --cflags` -I /home/nik2305/Documents/HEP_univ/ -I /home/nik2305/RooUnfold/src -lRooUnfold -L /home/nik2305/RooUnfold/build
// export LD_LIBRARY_PATH=/home/nik2305/RooUnfold/build:$LD_LIBRARY_PATH
using namespace std;
using ROOT::Math::PxPyPzMVector;

// --- helpers ---
TH1D* clone_and_set_style(TH1D* h) {
    TH1D* c = (TH1D*)h->Clone();
    return c;
}
TH1* tunfold(const TH1D* h_rec, const RooUnfoldResponse* h_resp_tunf) {
    auto unfold_reg = new RooUnfoldTUnfold(h_resp_tunf, h_rec);
    unfold_reg->OptimiseTau();
    TH1* h_unfold_reg = unfold_reg->Hunfold();
    h_unfold_reg->SetLineColor(3);
    return(h_unfold_reg);
}
// Unfold: uses response matrix (TH2D) to invert and apply to rec histogram
TH1D* unfold(const TH1D* h_rec, const RooUnfoldResponse* h_resp_tunf) {
    const TH2* h_resp = h_resp_tunf->Hresponse();
    int nbinsX = h_resp->GetNbinsX();
    int nbinsY = h_resp->GetNbinsY();
    h_resp->GetXaxis();
    if (nbinsX != nbinsY) {
        cerr << "Warning: response matrix not square (X=" << nbinsX << " Y=" << nbinsY << "). Attempting invert anyway may fail." << endl;
    }
    int N = nbinsX; // assume square

    // Build response matrix as in Python: mat_resp[j,i] = h_resp(i+1,j+1)
    TMatrixD mat_resp(N, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // careful: Python used mat_resp[j, i] = h_resp.GetBinContent(i+1, j+1)
            mat_resp(j, i) = h_resp->GetBinContent(i+1, j+1);
        }
    }

    // Invert using LU
    TDecompLU lu(mat_resp);
    TMatrixD mat_inv(N, N);
    if (!lu.Invert(mat_inv)) {
        cerr << "Matrix inversion failed!" << endl;
        return nullptr;
    }

    // Build rec vector (N x 1)
    TMatrixD vec_rec(N, 1);
    for (int i = 0; i < N; ++i) {
        vec_rec(i, 0) = h_rec->GetBinContent(i+1);
    }

    // Multiply: xsec = mat_inv * vec_rec
    TMatrixD mat_xsec = mat_inv * vec_rec;

    // Create histogram clone and fill content
    TH1D* h_xsec = (TH1D*)h_rec->Clone(Form("%s_xsec", h_rec->GetName()));
    for (int i = 0; i < N; ++i) {
        h_xsec->SetBinContent(i+1, mat_xsec(i, 0));
    }
    h_xsec->SetLineColor(8);
    h_xsec->SetLineStyle(2);
    return h_xsec;
}

// calc_recgen: produce gen and rec 1D histograms from 2D response
pair<TH1D*, TH1D*> calc_recgen(const RooUnfoldResponse* h2_resp_tunf, const char* name_prefix) {
    const TH2* h2_resp = h2_resp_tunf->Hresponse();
    int nbinsX = h2_resp->GetNbinsX();
    int nbinsY = h2_resp->GetNbinsY();

    // Create gen histogram: nbins = nbinsX, x axis from h2 x-axis
    double xlow = h2_resp->GetXaxis()->GetBinLowEdge(1);
    double xup = h2_resp->GetXaxis()->GetBinLowEdge(nbinsX + 1);
    TH1D* h_gen = new TH1D(Form("%s_gen", name_prefix), "", nbinsX, xlow, xup);

    // fill including under/overflow as in Python: loop bin index 0..nbinsX+1
    for (int binx = 0; binx <= nbinsX + 1; ++binx) {
        double sum = 0;
        for (int ibiny = 0; ibiny <= nbinsY + 1; ++ibiny) {
            sum += h2_resp->GetBinContent(binx, ibiny);
        }
        // set content to bin index 'binx' in histogram -> in ROOT hist bins are 1..nbins, with 0 = under, nbins+1 = overflow
        // For binx in 0..nbinsX+1 we map to histogram bin index = binx (same)
        if (binx >= 0 && binx <= nbinsX+1) {
            h_gen->SetBinContent(binx, sum);
        }
    }

    // rec histogram: nbins = nbinsY, axis from y-axis
    double ylow = h2_resp->GetYaxis()->GetBinLowEdge(1);
    double yup = h2_resp->GetYaxis()->GetBinLowEdge(nbinsY + 1);
    TH1D* h_rec = new TH1D(Form("%s_rec", name_prefix), "", nbinsY, ylow, yup);

    for (int biny = 0; biny <= nbinsY + 1; ++biny) {
        double sum = 0;
        for (int ibinx = 0; ibinx <= nbinsX + 1; ++ibinx) {
            sum += h2_resp->GetBinContent(ibinx, biny);
        }
        if (biny >= 0 && biny <= nbinsY + 1) {
            h_rec->SetBinContent(biny, sum);
        }
    }

    h_gen->SetLineColor(2);
    h_rec->SetLineColor(4);
    return {h_gen, h_rec};
}

// Normalize response rows (sum over Y for each X) similar to norm_response in Python
/*void norm_response(RooUnfoldResponse* h2_tunf) {
    const TH2* h2 = h2_tunf->Hresponse();
    h2->Sumw2();
    int nbinsX = h2->GetNbinsX();
    int nbinsY = h2->GetNbinsY();
    // loop ibinx 0..nbinsX+1 (include under/overflow)
    for (int ibinx = 0; ibinx <= nbinsX + 1; ++ibinx) {
        double norm = 0;
        for (int ibiny = 0; ibiny <= nbinsY + 1; ++ibiny) {
            norm += h2->GetBinContent(ibinx, ibiny);
        }
        if (norm != 0) {
            for (int ibiny = 0; ibiny <= nbinsY + 1; ++ibiny) {
                double c = h2->GetBinContent(ibinx, ibiny) / norm;
                h2->SetBinContent(ibinx, ibiny, c);
            }
        }
    }
}
*/
// purity & stability calculation


// draw ratios similar to Python function draw_ratios; draws on current pad
void draw_ratios(const TH1* h_gen, const TH1* h_rec, const TH1* h_xsec) {
    // h_gen_r = h_gen.Clone(); then divide by h_gen ??? original code was odd:
    // In python they did h_gen_r.Sumw2(); h_gen_r.Divide(h_gen) which leaves ones.
    // We'll reproduce shapes: draw reference histogram (h_gen) and then ratios rec/gen and xsec/gen.
    TH1* base = (TH1*)h_gen->Clone(Form("%s_ratio_base", h_gen->GetName()));
    base->Sumw2();
    base->Divide(base);
    // to produce ratio histograms:
    TH1* rec_r = (TH1*)h_rec->Clone(Form("%s_rec_ratio", h_rec->GetName()));
    rec_r->Sumw2();
    rec_r->Divide(h_gen);

    TH1* xsec_r = (TH1*)h_xsec->Clone(Form("%s_xsec_ratio", h_xsec->GetName()));
    xsec_r->Sumw2();
    xsec_r->Divide(h_gen);

    base->SetMinimum(0.5);
    base->SetMaximum(1.5);
    base->Draw("hist");
    rec_r->Draw("hist e1 same");
    xsec_r->Draw("hist same");

    // cleanup temporaries
    //delete base;
    //delete rec_r;
    //delete xsec_r;
}

// --- main ---
int main() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".2f");

    const int njets_max = 25;
    const double jet_delta_match = 0.5;

    // Open file and tree
    TFile* f = TFile::Open("ttbarSel_merged.root");
    if (!f || f->IsZombie()) {
        cerr << "Cannot open file ttbarSel_merged.root" << endl;
        return 1;
    }
    TTree* t = (TTree*)f->Get("tree");
    if (!t) {
        cerr << "Cannot find tree 'tree' in file" << endl;
        return 1;
    }
    t->SetBranchStatus("*", 0);

    // Branch arrays
    float mcB[4] = {0,0,0,0};
    float mcBbar[4] = {0,0,0,0};
    int Njet = 0;
    float jetPt[njets_max]; for (int i=0;i<njets_max;++i) jetPt[i]=0;
    float jetEta[njets_max]; for (int i=0;i<njets_max;++i) jetEta[i]=0;
    float jetPhi[njets_max]; for (int i=0;i<njets_max;++i) jetPhi[i]=0;

    // Activate and set branch addresses (names must match your tree)
    t->SetBranchStatus("mcB", 1);
    t->SetBranchAddress("mcB", mcB);
    t->SetBranchStatus("mcBbar", 1);
    t->SetBranchAddress("mcBbar", mcBbar);

    t->SetBranchStatus("Njet", 1);
    t->SetBranchAddress("Njet", &Njet);
    t->SetBranchStatus("jetPt", 1);
    t->SetBranchAddress("jetPt", jetPt);
    t->SetBranchStatus("jetEta", 1);
    t->SetBranchAddress("jetEta", jetEta);
    t->SetBranchStatus("jetPhi", 1);
    t->SetBranchAddress("jetPhi", jetPhi);

    // Histograms (response matrices)
    RooUnfoldResponse* h2_jet_resp_eta = new RooUnfoldResponse(8, -2.4, 2.4, 4, -2.4, 2.4);
    RooUnfoldResponse* h2_jet_resp_phi = new RooUnfoldResponse(8, 0, TMath::Pi(), 4, 0, TMath::Pi());
    RooUnfoldResponse* h2_jet_resp_pt  = new RooUnfoldResponse(8, 0, 400., 4, 0, 400.);

    Long64_t n = t->GetEntries();
    cout << "Number of events: " << n << endl;

    auto calculate_delta_r = [&](int gencharge, const float genp4[4], int nrec, const float rec_eta[], const float rec_phi[], const float rec_pt[]) {
        PxPyPzMVector genVector(genp4[0], genp4[1], genp4[2], genp4[3]);
        double best_delta_r = jet_delta_match;
        int best_rec = -1;
        for (int irec = 0; irec < nrec && irec < njets_max; ++irec) {
            if (rec_pt[irec] * gencharge < 0) continue;
            double delta_phi = fabs(genVector.Phi() - rec_phi[irec]);
            if (delta_phi > TMath::Pi()) delta_phi = 2 * TMath::Pi() - delta_phi;
            double delta_r = sqrt( pow(genVector.Eta() - rec_eta[irec], 2) + delta_phi*delta_phi );
            if (delta_r < best_delta_r) {
                best_delta_r = delta_r;
                best_rec = irec;
            }
        }
        if (best_rec >= 0) {
            h2_jet_resp_eta->Fill(genVector.Eta(), rec_eta[best_rec]);
            h2_jet_resp_phi->Fill(genVector.Phi(), rec_phi[best_rec]);
            h2_jet_resp_pt->Fill(genVector.Pt(), rec_pt[best_rec]);
        }
    };

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        calculate_delta_r(+1, mcB, Njet, jetEta, jetPhi, jetPt);
        calculate_delta_r(+1, mcBbar, Njet, jetEta, jetPhi, jetPt);
    }

    // Compute gen/rec projections and normalize response matrices
    auto p_eta = calc_recgen(h2_jet_resp_eta, "jet_eta");
    auto p_phi = calc_recgen(h2_jet_resp_phi, "jet_phi");
    auto p_pt  = calc_recgen(h2_jet_resp_pt,  "jet_pt");

    TH1D* h_gen_jet_eta = p_eta.first;
    TH1D* h_rec_jet_eta = p_eta.second;
    TH1D* h_gen_jet_phi = p_phi.first;
    TH1D* h_rec_jet_phi = p_phi.second;
    TH1D* h_gen_jet_pt  = p_pt.first;
    TH1D* h_rec_jet_pt  = p_pt.second;

    //norm_response(h2_jet_resp_eta);
    //norm_response(h2_jet_resp_phi);
   // norm_response(h2_jet_resp_pt);

    // Unfold


    // Draw response matrices to PDF
    TCanvas* c_resp = new TCanvas("c_resp_jet", "Jet response matrix", 900, 300);
    c_resp->Divide(3,1);
    c_resp->cd(1);
    auto h2_jet_resp_eta_draw = h2_jet_resp_eta->Hresponse();
    h2_jet_resp_eta_draw->Draw("colz");
    c_resp->cd(2);
    auto h2_jet_resp_phi_draw = h2_jet_resp_phi->Hresponse();
    h2_jet_resp_phi_draw->Draw("colz");
    c_resp->cd(3);
    auto h2_jet_resp_pt_draw = h2_jet_resp_pt->Hresponse();
    h2_jet_resp_pt_draw->Draw("colz");
    c_resp->SaveAs("response_jet.pdf");

    TH1* h_xsec_jet_eta = tunfold(h_rec_jet_eta, h2_jet_resp_eta);
    TH1* h_xsec_jet_phi = tunfold(h_rec_jet_phi, h2_jet_resp_phi);
    TH1* h_xsec_jet_pt  = tunfold(h_rec_jet_pt,  h2_jet_resp_pt);
    
    // Distributions and ratios canvas
    TCanvas* c_distr_jet = new TCanvas("c_distr_jet", "Distributions jet", 900, 600);
    c_distr_jet->Divide(3,2);

    c_distr_jet->cd(1);
    h_gen_jet_eta->SetMinimum(0);
    h_gen_jet_eta->Draw("hist");
    h_rec_jet_eta->Draw("hist e1 same");
    if (h_xsec_jet_eta) h_xsec_jet_eta->Draw("hist same");

    c_distr_jet->cd(2);
    h_gen_jet_phi->SetMinimum(0);
    h_gen_jet_phi->Draw("hist");
    h_rec_jet_phi->Draw("hist e1 same");
    if (h_xsec_jet_phi) h_xsec_jet_phi->Draw("hist same");

    c_distr_jet->cd(3);
    gPad->SetLogy();
    // set maximum as in Python: 1.05 * max of three hist maxima
    double maxval = 0;
    maxval = max(maxval, h_gen_jet_pt->GetMaximum());
    maxval = max(maxval, h_rec_jet_pt->GetMaximum());
    if (h_xsec_jet_pt) maxval = max(maxval, h_xsec_jet_pt->GetMaximum());
    h_gen_jet_pt->SetMaximum(1.05 * maxval);
    h_gen_jet_pt->SetMinimum(1);
    h_gen_jet_pt->Draw("hist");
    h_rec_jet_pt->Draw("hist e1 same");
    if (h_xsec_jet_pt) h_xsec_jet_pt->Draw("hist same");

    c_distr_jet->cd(4);
    draw_ratios(h_gen_jet_eta, h_rec_jet_eta, h_xsec_jet_eta);
    c_distr_jet->cd(5);
    draw_ratios(h_gen_jet_phi, h_rec_jet_phi, h_xsec_jet_phi);
    c_distr_jet->cd(6);
    draw_ratios(h_gen_jet_pt, h_rec_jet_pt, h_xsec_jet_pt);
    c_distr_jet->cd(4);
    gPad->SetLogy(0);
    draw_ratios(h_gen_jet_eta, h_rec_jet_eta, h_xsec_jet_eta);

    // PHI ratio
    c_distr_jet->cd(5);
    gPad->SetLogy(0);
    draw_ratios(h_gen_jet_phi, h_rec_jet_phi, h_xsec_jet_phi);

    // PT ratio (now guaranteed to exist)
    c_distr_jet->cd(6);
    gPad->SetLogy(0);
    draw_ratios(h_gen_jet_pt, h_rec_jet_pt, h_xsec_jet_pt);

    c_distr_jet->SaveAs("distr_jet.pdf");
    
    // Optionally: save histograms to file
    TFile out("unfolding_output.root", "RECREATE");
    h2_jet_resp_eta->Write();
    h2_jet_resp_phi->Write();
    h2_jet_resp_pt->Write();
    h_gen_jet_eta->Write();
    h_rec_jet_eta->Write();
    if (h_xsec_jet_eta) h_xsec_jet_eta->Write();
    h_gen_jet_phi->Write();
    h_rec_jet_phi->Write();
    if (h_xsec_jet_phi) h_xsec_jet_phi->Write();
    h_gen_jet_pt->Write();
    h_rec_jet_pt->Write();
    if (h_xsec_jet_pt) h_xsec_jet_pt->Write();
    out.Close();

    // cleanup (some pointers are owned by ROOT and will be cleaned on file close)
    delete c_resp;
    delete c_distr_jet;
    delete h2_jet_resp_eta;
    delete h2_jet_resp_phi;
    delete h2_jet_resp_pt;
    delete h_gen_jet_eta;
    delete h_rec_jet_eta;
    delete h_gen_jet_phi;
    delete h_rec_jet_phi;
    delete h_gen_jet_pt;
    delete h_rec_jet_pt;
    if (h_xsec_jet_eta) delete h_xsec_jet_eta;
    if (h_xsec_jet_phi) delete h_xsec_jet_phi;
    if (h_xsec_jet_pt) delete h_xsec_jet_pt;

    f->Close();
    delete f;

    cout << "Done." << endl;
    return 0;
}
