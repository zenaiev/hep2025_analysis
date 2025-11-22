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

using namespace std;
using ROOT::Math::PxPyPzMVector;
// fish: g++ matching.cxx -o matching (root-config --cflags --libs | string split " ") -I /home/nik2305/Documents/HEP_univ/ -I /home/nik2305/RooUnfold/src -lRooUnfold -L /home/nik2305/RooUnfold/build
//       set -x LD_LIBRARY_PATH /home/nik2305/RooUnfold/build $LD_LIBRARY_PATH
// bash: g++ matching.cxx -o matching `root-config --libs --cflags` -I /home/nik2305/Documents/HEP_univ/ -I /home/nik2305/RooUnfold/src -lRooUnfold -L /home/nik2305/RooUnfold/build
//       export LD_LIBRARY_PATH=/home/nik2305/RooUnfold/build:$LD_LIBRARY_PATH
// ./matching

TH1* tunfold(const TH1D* h_rec, const RooUnfoldResponse* h_resp_tunf, double tau, int color, int style) {
    auto unfold_reg = new RooUnfoldTUnfold(h_resp_tunf, h_rec);
    
    // reg parameter
    // 0.0 = unreg
    // >0.0 = reg
    unfold_reg->SetRegParm(tau); 
    
    TH1* h_res = unfold_reg->Hunfold();
    
    h_res->SetName(Form("%s_tau_%.2f", h_rec->GetName(), tau));
    h_res->SetLineColor(color);
    h_res->SetLineStyle(style);
    h_res->SetLineWidth(1);
    return h_res;
}

// draw ratios: reg and unreg
void draw_ratios(const TH1* h_gen, const TH1* h_rec, const TH1* h_unreg, const TH1* h_reg) {
    TH1* base = (TH1*)h_gen->Clone(Form("%s_ratio_base", h_gen->GetName()));
    base->Sumw2();
    base->Divide(base); 
    
    // Ratio reg
    TH1* r_unreg = nullptr;
    if (h_unreg) {
        r_unreg = (TH1*)h_unreg->Clone(Form("%s_unreg_r", h_unreg->GetName()));
        r_unreg->Divide(h_gen);
    }

    // Ratio unreg
    TH1* r_reg = nullptr;
    if (h_reg) {
        r_reg = (TH1*)h_reg->Clone(Form("%s_reg_r", h_reg->GetName()));
        r_reg->Divide(h_gen);
    }

    base->SetMinimum(0.5);
    base->SetMaximum(1.5);
    base->GetXaxis()->SetTitle(h_gen->GetXaxis()->GetTitle());
    base->GetYaxis()->SetTitle("Ratio to Gen");
    base->GetYaxis()->SetNdivisions(505);
    base->GetYaxis()->SetLabelSize(0.05);
    base->GetXaxis()->SetLabelSize(0.05);
    base->GetYaxis()->SetTitleSize(0.05);
    base->GetXaxis()->SetTitleSize(0.05);

    gPad->SetBottomMargin(0.15); 
    gPad->SetLeftMargin(0.18);    
    gPad->SetRightMargin(0.05);   


    base->SetMinimum(0.5);
    base->SetMaximum(1.5);
    
    base->GetXaxis()->SetTitle(h_gen->GetXaxis()->GetTitle());
    base->GetYaxis()->SetTitle("Ratio to Gen");

    base->GetYaxis()->SetNdivisions(505);

    base->GetYaxis()->SetLabelSize(0.05);
    base->GetXaxis()->SetLabelSize(0.05);
    base->GetYaxis()->SetTitleSize(0.06);
    base->GetXaxis()->SetTitleSize(0.06);

    base->GetXaxis()->SetTitleOffset(1.1); 
    base->GetYaxis()->SetTitleOffset(1.3); 
    
    base->Draw("hist"); // Лінія на 1
    if (r_unreg) r_unreg->Draw("hist same");
    if (r_reg)   r_reg->Draw("hist same");
}

// --- main ---
int main() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".2f");
    gStyle->SetOptTitle(0);

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

    // Activate and set branch addresses
    t->SetBranchStatus("mcB", 1); t->SetBranchAddress("mcB", mcB);
    t->SetBranchStatus("mcBbar", 1); t->SetBranchAddress("mcBbar", mcBbar);
    t->SetBranchStatus("Njet", 1); t->SetBranchAddress("Njet", &Njet);
    t->SetBranchStatus("jetPt", 1); t->SetBranchAddress("jetPt", jetPt);
    t->SetBranchStatus("jetEta", 1); t->SetBranchAddress("jetEta", jetEta);
    t->SetBranchStatus("jetPhi", 1); t->SetBranchAddress("jetPhi", jetPhi);

    // Histograms (response matrices)
    RooUnfoldResponse* h2_jet_resp_eta = new RooUnfoldResponse(10, -2.4, 2.4, 10, -2.4, 2.4);
    RooUnfoldResponse* h2_jet_resp_phi = new RooUnfoldResponse(10, 0, TMath::Pi(), 10, 0, TMath::Pi());
    RooUnfoldResponse* h2_jet_resp_pt  = new RooUnfoldResponse(10, 0, 400., 10, 0, 400.);

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
            h2_jet_resp_eta->Fill(rec_eta[best_rec], genVector.Eta());
            h2_jet_resp_phi->Fill(rec_phi[best_rec], genVector.Phi());
            h2_jet_resp_pt->Fill(rec_pt[best_rec], genVector.Pt());
        }
    };

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        calculate_delta_r(+1, mcB, Njet, jetEta, jetPhi, jetPt);
        calculate_delta_r(+1, mcBbar, Njet, jetEta, jetPhi, jetPt);
    }

    TH1D* h_rec_jet_eta = (TH1D*)h2_jet_resp_eta->Hmeasured()->Clone("jet_eta_rec");
    TH1D* h_rec_jet_phi = (TH1D*)h2_jet_resp_phi->Hmeasured()->Clone("jet_phi_rec");
    TH1D* h_rec_jet_pt  = (TH1D*)h2_jet_resp_pt->Hmeasured()->Clone("jet_pt_rec");

    TH1D* h_gen_jet_eta = (TH1D*)h2_jet_resp_eta->Htruth()->Clone("jet_eta_gen");
    TH1D* h_gen_jet_phi = (TH1D*)h2_jet_resp_phi->Htruth()->Clone("jet_phi_gen");
    TH1D* h_gen_jet_pt  = (TH1D*)h2_jet_resp_pt->Htruth()->Clone("jet_pt_gen");

    h_gen_jet_eta->GetXaxis()->SetTitle("#eta");
    h_gen_jet_phi->GetXaxis()->SetTitle("#phi");
    h_gen_jet_pt->GetXaxis()->SetTitle("p_{T} [GeV]");

    h_gen_jet_eta->SetLineColor(kRed);    
    h_rec_jet_eta->SetLineColor(kBlue);
    h_rec_jet_eta->SetLineWidth(2); 
    h_rec_jet_eta->SetLineStyle(2);
    h_gen_jet_eta->SetLineWidth(2);
    
    h_gen_jet_phi->SetLineColor(kRed);     
    h_gen_jet_phi->SetLineWidth(2);
    h_rec_jet_phi->SetLineColor(kBlue);    
    h_rec_jet_phi->SetLineWidth(2); 
    h_rec_jet_phi->SetLineStyle(2);
    
    h_gen_jet_pt->SetLineColor(kRed);      
    h_gen_jet_pt->SetLineWidth(2);
    h_rec_jet_pt->SetLineColor(kBlue);     
    h_rec_jet_pt->SetLineWidth(2);  
    h_rec_jet_pt->SetLineStyle(2);

    // Unfold
    // 1. Unregularized (Tau = 0.0)
    TH1* h_unreg_eta = tunfold(h_rec_jet_eta, h2_jet_resp_eta, 0.0, 6, 1);
    TH1* h_unreg_phi = tunfold(h_rec_jet_phi, h2_jet_resp_phi, 0.0, 6, 1);
    TH1* h_unreg_pt  = tunfold(h_rec_jet_pt,  h2_jet_resp_pt,  0.0, 6, 1);

    // 2. Regularized (Tau != 0.0)
    TH1* h_reg_eta = tunfold(h_rec_jet_eta, h2_jet_resp_eta, 0.001, 8, 1);
    TH1* h_reg_phi = tunfold(h_rec_jet_phi, h2_jet_resp_phi, 0.001, 8, 1);
    TH1* h_reg_pt  = tunfold(h_rec_jet_pt,  h2_jet_resp_pt,  0.0005, 8, 1);
    
    // Distributions and ratios canvas
    TCanvas* c_distr_jet = new TCanvas("c_distr_jet", "Distributions jet", 1000, 800);
    c_distr_jet->Divide(3,2);

    TLegend* leg = new TLegend(0.6, 0.7, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->AddEntry(h_gen_jet_eta, "Gen", "l");
    leg->AddEntry(h_rec_jet_eta, "Rec", "l");
    leg->AddEntry(h_unreg_eta, "Unreg", "l");
    leg->AddEntry(h_reg_eta,   "Reg", "l");

    //PLOT ETA ---
    c_distr_jet->cd(1);
    h_gen_jet_eta->SetMinimum(0);
    h_gen_jet_eta->Draw("hist");
    h_rec_jet_eta->Draw("hist same");
    if (h_unreg_eta) h_unreg_eta->Draw("hist same");
    if (h_reg_eta)   h_reg_eta->Draw("hist same");

    //PLOT PHI ---
    c_distr_jet->cd(2);
    h_gen_jet_phi->SetMinimum(0);
    h_gen_jet_phi->Draw("hist");
    h_rec_jet_phi->Draw("hist same");
    if (h_unreg_phi) h_unreg_phi->Draw("hist same");
    if (h_reg_phi)   h_reg_phi->Draw("hist same");

    //PLOT PT ---
    c_distr_jet->cd(3);
    gPad->SetLogy();
    double maxval = 0;
    maxval = max(maxval, h_gen_jet_pt->GetMaximum());
    maxval = max(maxval, h_rec_jet_pt->GetMaximum());
    if (h_unreg_pt) maxval = max(maxval, h_unreg_pt->GetMaximum());
    h_gen_jet_pt->SetMaximum(2.0 * maxval); 
    h_gen_jet_pt->SetMinimum(1);
    h_gen_jet_pt->Draw("hist");
    h_rec_jet_pt->Draw("hist same");
    if (h_unreg_pt) h_unreg_pt->Draw("hist same");
    if (h_reg_pt)   h_reg_pt->Draw("hist same");
    leg->Draw();

    //RATIOS
    c_distr_jet->cd(4);
    gPad->SetLogy(0);
    draw_ratios(h_gen_jet_eta, h_rec_jet_eta, h_unreg_eta, h_reg_eta);

    c_distr_jet->cd(5);
    gPad->SetLogy(0);
    draw_ratios(h_gen_jet_phi, h_rec_jet_phi, h_unreg_phi, h_reg_phi);

    c_distr_jet->cd(6);
    gPad->SetLogy(0);
    draw_ratios(h_gen_jet_pt, h_rec_jet_pt, h_unreg_pt, h_reg_pt);

    c_distr_jet->SaveAs("distr_jet.pdf");
    
    // cleanup
    delete c_distr_jet;
    delete h2_jet_resp_eta; delete h2_jet_resp_phi; delete h2_jet_resp_pt;
    delete h_gen_jet_eta; delete h_rec_jet_eta;
    delete h_gen_jet_phi; delete h_rec_jet_phi;
    delete h_gen_jet_pt;  delete h_rec_jet_pt;
    
    if (h_unreg_eta) delete h_unreg_eta;
    if (h_unreg_phi) delete h_unreg_phi;
    if (h_unreg_pt)  delete h_unreg_pt;
    
    if (h_reg_eta) delete h_reg_eta;
    if (h_reg_phi) delete h_reg_phi;
    if (h_reg_pt)  delete h_reg_pt;

    f->Close();
    delete f;

    cout << "Done." << endl;
    return 0;
}