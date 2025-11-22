// matchingv2.cpp
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

// --- COMMANDS TO RUN ---
// Fish:
// g++ matchingv2.cxx -o matchingv2 (root-config --cflags --libs | string split " ") -I /home/nik2305/Documents/HEP_univ/ -I /home/nik2305/RooUnfold/src -lRooUnfold -L /home/nik2305/RooUnfold/build
// set -x LD_LIBRARY_PATH /home/nik2305/RooUnfold/build $LD_LIBRARY_PATH
// ./matchingv2

// Bash:
// g++ matchingv2.cxx -o matchingv2 `root-config --libs --cflags` -I /home/nik2305/Documents/HEP_univ/ -I /home/nik2305/RooUnfold/src -lRooUnfold -L /home/nik2305/RooUnfold/build
// export LD_LIBRARY_PATH=/home/nik2305/RooUnfold/build:$LD_LIBRARY_PATH
// ./matchingv2

TH1* tunfold(const TH1D* h_rec, const RooUnfoldResponse* h_resp_tunf, double tau, int color, int style) {
    auto unfold_reg = new RooUnfoldTUnfold(h_resp_tunf, h_rec);
    // reg parameter
    // 0.0 = unreg
    // >0.0 = reg
    unfold_reg->SetRegParm(tau); 
    
    TH1* h_res = unfold_reg->Hunfold();
    
    h_res->SetName(Form("%s_tau_%.5f", h_rec->GetName(), tau));
    h_res->SetLineColor(color);
    h_res->SetLineStyle(style);
    h_res->SetLineWidth(2);
    return h_res;
}

void draw_ratios(const TH1* h_gen, const TH1* h_rec, const TH1* h_unreg, const TH1* h_reg) {
    
    TH1* base = (TH1*)h_gen->Clone(Form("%s_ratio_base", h_gen->GetName()));
    base->Sumw2();
    base->Divide(base); 
    base->SetLineColor(kBlack);
    base->SetLineStyle(2); 
    base->SetLineWidth(1);
    

    TH1* r_unreg = nullptr;
    if (h_unreg) {
        r_unreg = (TH1*)h_unreg->Clone(Form("%s_unreg_r", h_unreg->GetName()));
        r_unreg->Divide(h_gen);
        
    }

   
    TH1* r_reg = nullptr;
    if (h_reg) {
        r_reg = (TH1*)h_reg->Clone(Form("%s_reg_r", h_reg->GetName()));
        r_reg->Divide(h_gen);
    }


    gPad->SetBottomMargin(0.15); 
    gPad->SetLeftMargin(0.15);    
    gPad->SetRightMargin(0.05);   

    base->SetMinimum(-1.0); 
    base->SetMaximum(3.0);
    
    base->GetXaxis()->SetTitle(h_gen->GetXaxis()->GetTitle());
    base->GetYaxis()->SetTitle("Ratio to Gen");
    base->GetYaxis()->SetNdivisions(505);


    base->GetYaxis()->SetLabelSize(0.05);
    base->GetXaxis()->SetLabelSize(0.05);
    base->GetYaxis()->SetTitleSize(0.06);
    base->GetXaxis()->SetTitleSize(0.06);

  
    base->GetXaxis()->SetTitleOffset(1.1); 
    base->GetYaxis()->SetTitleOffset(1.1);
    

    base->Draw("hist"); 
    

    if (r_unreg) r_unreg->Draw("same E1"); 
    if (r_reg)   r_reg->Draw("same E1");
}

int main() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat(".2f");
    gStyle->SetOptTitle(0);

    const int njets_max = 25;
    const double jet_delta_match = 0.5;


    TFile* f = TFile::Open("ttbarSel_merged.root");
    if (!f || f->IsZombie()) { cerr << "Error opening file" << endl; return 1; }
    TTree* t = (TTree*)f->Get("tree");
    if (!t) { cerr << "Error finding tree" << endl; return 1; }
    t->SetBranchStatus("*", 0);

    // Variables
    float mcB[4] = {0}, mcBbar[4] = {0};
    int Njet = 0;
    float jetPt[njets_max] = {0}, jetEta[njets_max] = {0}, jetPhi[njets_max] = {0};

    t->SetBranchStatus("mcB", 1); t->SetBranchAddress("mcB", mcB);
    t->SetBranchStatus("mcBbar", 1); t->SetBranchAddress("mcBbar", mcBbar);
    t->SetBranchStatus("Njet", 1); t->SetBranchAddress("Njet", &Njet);
    t->SetBranchStatus("jetPt", 1); t->SetBranchAddress("jetPt", jetPt);
    t->SetBranchStatus("jetEta", 1); t->SetBranchAddress("jetEta", jetEta);
    t->SetBranchStatus("jetPhi", 1); t->SetBranchAddress("jetPhi", jetPhi);

    // Binning
    
    // ETA
    int nRec_eta = 8; 
    double min_eta = -2.4; 
    double max_eta = 2.4;
    int nGen_eta = 8; 
    
    // PHI
    int nRec_phi = 8; 
    double min_phi = -TMath::Pi(); 
    double max_phi = TMath::Pi();
    int nGen_phi = 8;

    // PT
    int nRec_pt = 12; 
    double min_pt = 30.0; 
    double max_pt = 330.0;
    int nGen_pt = 6;

    // Response 
    RooUnfoldResponse* h2_jet_resp_eta = new RooUnfoldResponse(nRec_eta, min_eta, max_eta, nGen_eta, min_eta, max_eta);
    RooUnfoldResponse* h2_jet_resp_phi = new RooUnfoldResponse(nRec_phi, min_phi, max_phi, nGen_phi, min_phi, max_phi);
    RooUnfoldResponse* h2_jet_resp_pt  = new RooUnfoldResponse(nRec_pt, min_pt, max_pt, nGen_pt, min_pt, max_pt);

    //Independent Histograms (Test Half)
    TH1D* h_rec_jet_eta = new TH1D("jet_eta_rec", "", nRec_eta, min_eta, max_eta); h_rec_jet_eta->Sumw2();
    TH1D* h_gen_jet_eta = new TH1D("jet_eta_gen", "", nGen_eta, min_eta, max_eta); h_gen_jet_eta->Sumw2();

    TH1D* h_rec_jet_phi = new TH1D("jet_phi_rec", "", nRec_phi, min_phi, max_phi); h_rec_jet_phi->Sumw2();
    TH1D* h_gen_jet_phi = new TH1D("jet_phi_gen", "", nGen_phi, min_phi, max_phi); h_gen_jet_phi->Sumw2();

    TH1D* h_rec_jet_pt  = new TH1D("jet_pt_rec",  "", nRec_pt, min_pt, max_pt); h_rec_jet_pt->Sumw2();
    TH1D* h_gen_jet_pt  = new TH1D("jet_pt_gen",  "", nGen_pt, min_pt, max_pt); h_gen_jet_pt->Sumw2();

    Long64_t n = t->GetEntries();
    Long64_t n_half = n / 2; 
    cout << "Total events: " << n << endl;
    cout << "Training (Response): " << n_half << endl;
    cout << "Test (Pseudo-Data): " << n - n_half << endl;


    auto calculate_delta_r = [&](int gencharge, const float genp4[4], int nrec, 
                                 const float rec_eta[], const float rec_phi[], const float rec_pt[], 
                                 bool is_training) {
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
            bool in_pt_range  = (genVector.Pt() >= min_pt && genVector.Pt() < max_pt && rec_pt[best_rec] >= min_pt && rec_pt[best_rec] < max_pt);
            bool in_eta_range = (genVector.Eta() >= min_eta && genVector.Eta() < max_eta && rec_eta[best_rec] >= min_eta && rec_eta[best_rec] < max_eta);
            
            if (in_pt_range && in_eta_range) {
                if (is_training) {
                    h2_jet_resp_eta->Fill(rec_eta[best_rec], genVector.Eta());
                    h2_jet_resp_phi->Fill(rec_phi[best_rec], genVector.Phi());
                    h2_jet_resp_pt->Fill(rec_pt[best_rec], genVector.Pt());
                } else {
                    h_rec_jet_eta->Fill(rec_eta[best_rec]);
                    h_gen_jet_eta->Fill(genVector.Eta());
                    h_rec_jet_phi->Fill(rec_phi[best_rec]);
                    h_gen_jet_phi->Fill(genVector.Phi());
                    h_rec_jet_pt->Fill(rec_pt[best_rec]);
                    h_gen_jet_pt->Fill(genVector.Pt());
                }
            }
        }
    };

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        bool is_training = (i < n_half);
        calculate_delta_r(+1, mcB, Njet, jetEta, jetPhi, jetPt, is_training);
        calculate_delta_r(+1, mcBbar, Njet, jetEta, jetPhi, jetPt, is_training);
    }

    h_gen_jet_eta->GetXaxis()->SetTitle("#eta");
    h_gen_jet_phi->GetXaxis()->SetTitle("#phi");
    h_gen_jet_pt->GetXaxis()->SetTitle("p_{T} [GeV]");

    h_gen_jet_eta->SetLineColor(kRed); h_gen_jet_eta->SetLineWidth(2);
    h_gen_jet_phi->SetLineColor(kRed); h_gen_jet_phi->SetLineWidth(2);
    h_gen_jet_pt->SetLineColor(kRed);  h_gen_jet_pt->SetLineWidth(2);

    h_rec_jet_eta->SetLineColor(kBlue); h_rec_jet_eta->SetLineWidth(2); h_rec_jet_eta->SetLineStyle(2);
    h_rec_jet_phi->SetLineColor(kBlue); h_rec_jet_phi->SetLineWidth(2); h_rec_jet_phi->SetLineStyle(2);
    h_rec_jet_pt->SetLineColor(kBlue);  h_rec_jet_pt->SetLineWidth(2);  h_rec_jet_pt->SetLineStyle(2);

    //Unfolding 
   
    TH1* h_unreg_eta = tunfold(h_rec_jet_eta, h2_jet_resp_eta, 0.0, kGreen+2, 1);
    h_unreg_eta->SetMarkerStyle(20); h_unreg_eta->SetMarkerSize(0.8);

    TH1* h_unreg_phi = tunfold(h_rec_jet_phi, h2_jet_resp_phi, 0.0, kGreen+2, 1);
    h_unreg_phi->SetMarkerStyle(20); h_unreg_phi->SetMarkerSize(0.8);

    TH1* h_unreg_pt  = tunfold(h_rec_jet_pt,  h2_jet_resp_pt,  0.0, kGreen+2, 1);
    h_unreg_pt->SetMarkerStyle(20);  h_unreg_pt->SetMarkerSize(0.8);

    TH1* h_reg_eta = tunfold(h_rec_jet_eta, h2_jet_resp_eta, 0.001, kMagenta, 1);
    h_reg_eta->SetMarkerStyle(20); h_reg_eta->SetMarkerSize(0.8);

    TH1* h_reg_phi = tunfold(h_rec_jet_phi, h2_jet_resp_phi, 0.001, kMagenta, 1);
    h_reg_phi->SetMarkerStyle(20); h_reg_phi->SetMarkerSize(0.8);

    TH1* h_reg_pt  = tunfold(h_rec_jet_pt,  h2_jet_resp_pt,  0.01, kMagenta, 1);
    h_reg_pt->SetMarkerStyle(20);  h_reg_pt->SetMarkerSize(0.8);
    

    TCanvas* c_distr_jet = new TCanvas("c_distr_jet", "Distributions jet", 1000, 800);
    c_distr_jet->Divide(3,2);

    // Plot ETA
    c_distr_jet->cd(1); gPad->SetLeftMargin(0.15);
    h_gen_jet_eta->SetMinimum(0);
    h_gen_jet_eta->GetYaxis()->SetTitleOffset(1.6);
    h_gen_jet_eta->Draw("hist");            
    h_rec_jet_eta->Draw("hist same E1");    
    if (h_unreg_eta) h_unreg_eta->Draw("same E1");
    if (h_reg_eta)   h_reg_eta->Draw("same E1");

  
    TLegend* leg = new TLegend(0.6, 0.7, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->AddEntry(h_gen_jet_eta, "Gen (Test)", "l");
    leg->AddEntry(h_rec_jet_eta, "Rec (Test)", "lep");
    leg->AddEntry(h_unreg_eta, "Unreg", "lep");
    leg->AddEntry(h_reg_eta,   "Reg", "lep");
    leg->Draw();

    //Plot PHI 
    c_distr_jet->cd(2); gPad->SetLeftMargin(0.15);
    h_gen_jet_phi->SetMinimum(0);
    h_gen_jet_phi->GetYaxis()->SetTitleOffset(1.6);
    h_gen_jet_phi->Draw("hist");
    h_rec_jet_phi->Draw("hist same E1");
    if (h_unreg_phi) h_unreg_phi->Draw("same E1");
    if (h_reg_phi)   h_reg_phi->Draw("same E1");
    
    // Local Legend for Phi (Bottom)
    TLegend* leg_phi = new TLegend(0.55, 0.15, 0.89, 0.35);
    leg_phi->SetBorderSize(0);
    leg_phi->AddEntry(h_gen_jet_eta, "Gen", "l");
    leg_phi->AddEntry(h_rec_jet_eta, "Rec", "lep");
    leg_phi->AddEntry(h_unreg_eta, "Unreg", "lep");
    leg_phi->AddEntry(h_reg_eta,   "Reg", "lep");
    leg_phi->Draw();

    //Plot PT
    c_distr_jet->cd(3); gPad->SetLeftMargin(0.15); gPad->SetLogy();
    double maxval = max(h_gen_jet_pt->GetMaximum(), h_rec_jet_pt->GetMaximum());
    if (h_unreg_pt) maxval = max(maxval, h_unreg_pt->GetMaximum());
    h_gen_jet_pt->SetMaximum(10.0 * maxval); 
    h_gen_jet_pt->SetMinimum(1);
    h_gen_jet_pt->GetYaxis()->SetTitleOffset(1.6);
    h_gen_jet_pt->Draw("hist");
    h_rec_jet_pt->Draw("hist same E1");
    if (h_unreg_pt) h_unreg_pt->Draw("same E1");
    if (h_reg_pt)   h_reg_pt->Draw("same E1");
    leg->Draw(); 

    //Ratios
    c_distr_jet->cd(4); gPad->SetLogy(0);
    draw_ratios(h_gen_jet_eta, h_rec_jet_eta, h_unreg_eta, h_reg_eta);

    c_distr_jet->cd(5); gPad->SetLogy(0);
    draw_ratios(h_gen_jet_phi, h_rec_jet_phi, h_unreg_phi, h_reg_phi);

    c_distr_jet->cd(6); gPad->SetLogy(0);
    draw_ratios(h_gen_jet_pt, h_rec_jet_pt, h_unreg_pt, h_reg_pt);

    c_distr_jet->SaveAs("distr_jetv2.pdf");
    
    // Cleanup
    delete c_distr_jet;
    delete h2_jet_resp_eta; delete h2_jet_resp_phi; delete h2_jet_resp_pt;
    delete h_rec_jet_eta; delete h_gen_jet_eta;
    delete h_rec_jet_phi; delete h_gen_jet_phi;
    delete h_rec_jet_pt;  delete h_gen_jet_pt;
    
    if (h_unreg_eta) delete h_unreg_eta; if (h_unreg_phi) delete h_unreg_phi; if (h_unreg_pt) delete h_unreg_pt;
    if (h_reg_eta) delete h_reg_eta; if (h_reg_phi) delete h_reg_phi; if (h_reg_pt) delete h_reg_pt;

    f->Close();
    delete f;
    cout << "Done." << endl;
    return 0;
}