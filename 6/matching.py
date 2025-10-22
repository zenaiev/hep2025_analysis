import ROOT
import array

def add_branch(t, branch_name, branch_type, branch_size):
  br = array.array(branch_type, [0]*branch_size)
  t.SetBranchStatus(branch_name, 1)
  t.SetBranchAddress(branch_name, br)
  return br

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat('.2f')
nrepeat = 1
f = ROOT.TFile('ttbarSel_merged.root')
t = f.Get('tree')
t.SetBranchStatus("*", 0)
nleptons_max = 5
njets_max = 25
lep_delta_match = 0.1
jet_delta_match = 0.5
mcLp = add_branch(t, 'mcLp', 'f', 4)
mcLm = add_branch(t, 'mcLm', 'f', 4)
mcB = add_branch(t, 'mcB', 'f', 4)
mcBbar = add_branch(t, 'mcBbar', 'f', 4)
Nmu = add_branch(t, 'Nmu', 'i', 1)
muPt = add_branch(t, 'muPt', 'f', nleptons_max)
muEta = add_branch(t, 'muEta', 'f', nleptons_max)
muPhi = add_branch(t, 'muPhi', 'f', nleptons_max)
Nel = add_branch(t, 'Nel', 'i', 1)
elPt = add_branch(t, 'elPt', 'f', nleptons_max)
elEta = add_branch(t, 'elEta', 'f', nleptons_max)
elPhi = add_branch(t, 'elPhi', 'f', nleptons_max)
Njet = add_branch(t, 'Njet', 'i', 1)
jetPt = add_branch(t, 'jetPt', 'f', njets_max)
jetEta = add_branch(t, 'jetEta', 'f', njets_max)
jetPhi = add_branch(t, 'jetPhi', 'f', njets_max)

h_lep_delta_r = ROOT.TH1D('h_lep_delta_r', 'lepton delta R', 500, 0., 5.)
h_jet_delta_r = ROOT.TH1D('h_jet_delta_r', 'jet delta R', 500, 0., 5.)
g_lep_eta = ROOT.TGraph()
g_lep_phi = ROOT.TGraph()
g_lep_pt = ROOT.TGraph()
g_jet_eta = ROOT.TGraph()
g_jet_phi = ROOT.TGraph()
g_jet_pt = ROOT.TGraph()
h2_lep_resp_eta = ROOT.TH2D('h2_lep_resp_eta', 'lepton response matrix eta', 15, -2.4, 2.4, 15, -2.4, 2.4)
h2_lep_resp_phi = ROOT.TH2D('h2_lep_resp_phi', 'lepton response matrix phi', 15, 0, ROOT.TMath.Pi(), 15, 0, ROOT.TMath.Pi())
h2_lep_resp_pt = ROOT.TH2D('h2_lep_resp_pt', 'lepton response matrix pt', 15, 0, 400., 15, 0, 400.)
h2_jet_resp_eta = ROOT.TH2D('h2_jet_resp_eta', 'jet response matrix eta', 15, -2.4, 2.4, 15, -2.4, 2.4)
h2_jet_resp_phi = ROOT.TH2D('h2_jet_resp_phi', 'jet response matrix phi', 15, 0, ROOT.TMath.Pi(), 15, 0, ROOT.TMath.Pi())
h2_jet_resp_pt = ROOT.TH2D('h2_jet_resp_pt', 'jet response matrix pt', 15, 0, 400., 15, 0, 400.)

n = t.GetEntries()
print(f'Number of events: {n}')
for irepeat in range(nrepeat):
  for i in range(n):
    t.GetEntry(i)
    def calculate_delta_r(gencharge, genp4, nrec, rec_eta, rec_phi, rec_pt, is_lepton=True):
      genVector = ROOT.Math.PxPyPzMVector(genp4[0], genp4[1], genp4[2], genp4[3])
      best_delta_r = lep_delta_match
      best_rec = -1
      for irec in range(nrec):
        if rec_pt[irec] * gencharge < 0:
          continue
        delta_phi = abs(genVector.Phi() - rec_phi[irec])
        if delta_phi > ROOT.TMath.Pi():
          delta_phi = 2 * ROOT.TMath.Pi() - delta_phi
        delta_r = ((genVector.Eta() - rec_eta[irec])**2 + delta_phi**2)**0.5
        if delta_r < best_delta_r:
          best_delta_r = delta_r
          best_rec = irec
        if is_lepton:
          h_lep_delta_r.Fill(delta_r)
        else:
          h_jet_delta_r.Fill(delta_r)
      if best_rec >= 0:
        if is_lepton:
          g_lep_eta.AddPoint(genVector.Eta(), rec_eta[best_rec])
          g_lep_phi.AddPoint(genVector.Phi(), rec_phi[best_rec])
          g_lep_pt.AddPoint(genVector.Pt(), abs(rec_pt[best_rec]))
          h2_lep_resp_eta.Fill(genVector.Eta(), rec_eta[best_rec])
          h2_lep_resp_phi.Fill(genVector.Phi(), rec_phi[best_rec])
          h2_lep_resp_pt.Fill(genVector.Pt(), rec_pt[best_rec])
        else:
          g_jet_eta.AddPoint(genVector.Eta(), rec_eta[best_rec])
          g_jet_phi.AddPoint(genVector.Phi(), rec_phi[best_rec])
          g_jet_pt.AddPoint(genVector.Pt(), abs(rec_pt[best_rec]))
          h2_jet_resp_eta.Fill(genVector.Eta(), rec_eta[best_rec])
          h2_jet_resp_phi.Fill(genVector.Phi(), rec_phi[best_rec])
          h2_jet_resp_pt.Fill(genVector.Pt(), rec_pt[best_rec])
    calculate_delta_r(1, mcLp, Nel[0], elEta, elPhi, elPt)
    calculate_delta_r(-1, mcLm, Nel[0], elEta, elPhi, elPt)
    calculate_delta_r(1, mcLp, Nmu[0], muEta, muPhi, muPt)
    calculate_delta_r(-1, mcLm, Nmu[0], muEta, muPhi, muPt)
    calculate_delta_r(1, mcB, Njet[0], jetEta, jetPhi, jetPt, is_lepton=False)
    calculate_delta_r(1, mcBbar, Njet[0], jetEta, jetPhi, jetPt, is_lepton=False)

h_lep_delta_r.SetLineColor(2)
c = ROOT.TCanvas('c', 'lepton delta R', 800, 400)
c.Divide(2,1)
c.cd(1)
h_lep_delta_r.Draw()
c.cd(2)
ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()
h_lep_delta_r.Draw()
c.SaveAs('deltaR.pdf')

h_jet_delta_r.SetLineColor(2)
c = ROOT.TCanvas('c_jet', 'jet delta R', 800, 400)
c.Divide(2,1)
c.cd(1)
h_jet_delta_r.Draw()
c.cd(2)
ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()
h_jet_delta_r.Draw()
c.SaveAs('deltaR_jet.pdf')

cc = ROOT.TCanvas('cc', 'lepton delta eta, phi, pT', 900, 300)
cc.Divide(3,1)
cc.cd(1)
g_lep_eta.Draw('ap')
cc.cd(2)
g_lep_phi.Draw('ap')
cc.cd(3)
g_lep_pt.Draw('ap')
cc.SaveAs('delta_etaphipt.pdf')

cc = ROOT.TCanvas('cc_jet', 'jet delta eta, phi, pT', 900, 300)
cc.Divide(3,1)
cc.cd(1)
g_jet_eta.Draw('ap')
cc.cd(2)
g_jet_phi.Draw('ap')
cc.cd(3)
g_jet_pt.Draw('ap')
cc.SaveAs('delta_etaphipt_jet.pdf')
    
def norm_response(h2):
  h2.Sumw2()
  for ibinx in range(h2.GetNbinsX()):
    norm = sum(h2.GetBinContent(ibinx+1, ibiny+1) for ibiny in range(h2.GetNbinsY()))
    if norm != 0:
      for ibiny in range(h2.GetNbinsY()):
        h2.SetBinContent(ibinx+1, ibiny+1, h2.GetBinContent(ibinx+1, ibiny+1) / norm)
norm_response(h2_lep_resp_eta)
norm_response(h2_lep_resp_phi)
norm_response(h2_lep_resp_pt)
c_resp = ROOT.TCanvas('c_resp', 'Lepton response matrix', 900, 300)
c_resp.Divide(3,1)
c_resp.cd(1)
h2_lep_resp_eta.Draw('colztext')
c_resp.cd(2)
h2_lep_resp_phi.Draw('colztext')
c_resp.cd(3)
h2_lep_resp_pt.Draw('colztext')
c_resp.SaveAs('response_lep.pdf')

norm_response(h2_jet_resp_eta)
norm_response(h2_jet_resp_phi)
norm_response(h2_jet_resp_pt)
c_resp = ROOT.TCanvas('c_resp_jet', 'Jet response matrix', 900, 300)
c_resp.Divide(3,1)
c_resp.cd(1)
h2_jet_resp_eta.Draw('colztext')
c_resp.cd(2)
h2_jet_resp_phi.Draw('colztext')
c_resp.cd(3)
h2_jet_resp_pt.Draw('colztext')
c_resp.SaveAs('response_jet.pdf')

def calc_pursta(h_resp):
  h_pur = ROOT.TH1D('', '', h_resp.GetNbinsY(), h_resp.GetYaxis().GetBinLowEdge(1), h_resp.GetYaxis().GetBinLowEdge(h_resp.GetNbinsY() + 1))
  for j in range(h_resp.GetNbinsY()):
    allrec = sum(h_resp.GetBinContent(irec+1, j+1) for irec in range(h_resp.GetNbinsX()))
    if allrec == 0:
      allrec = 1
    h_pur.SetBinContent(j+1, h_resp.GetBinContent(j+1, j+1) / allrec)
  h_sta = ROOT.TH1D('', '', h_resp.GetNbinsX(), h_resp.GetXaxis().GetBinLowEdge(1), h_resp.GetXaxis().GetBinLowEdge(h_resp.GetNbinsX() + 1))
  for i in range(h_resp.GetNbinsX()):
    allgen = sum(h_resp.GetBinContent(i+1, jrec+1) for jrec in range(h_resp.GetNbinsY()))
    if allgen == 0:
      allgen = 1
    h_sta.SetBinContent(i+1, h_resp.GetBinContent(i+1, i+1) / allgen)
  h_pur.SetLineColor(2)
  h_sta.SetLineColor(4)
  for h in [h_pur,h_sta]:
    h.SetMinimum(0)
    h.SetMaximum(1.1)
  return h_pur, h_sta

c_pursta_lep = ROOT.TCanvas('c_pursta_lep', 'Distributions lep', 900, 300)
c_pursta_lep.Divide(3,1)
c_pursta_lep.cd(1)
h_pur_eta_lep, h_sta_eta_lep = calc_pursta(h2_lep_resp_eta)
h_pur_eta_lep.Draw('hist')
h_sta_eta_lep.Draw('hist same')
c_pursta_lep.cd(2)
h_pur_phi_lep, h_sta_phi_lep = calc_pursta(h2_lep_resp_phi)
h_pur_phi_lep.Draw('hist')
h_sta_phi_lep.Draw('hist same')
c_pursta_lep.cd(3)
h_pur_pt_lep, h_sta_pt_lep = calc_pursta(h2_lep_resp_pt)
h_pur_pt_lep.Draw('hist')
h_sta_pt_lep.Draw('hist same')
c_pursta_lep.SaveAs('pursta_lep.pdf')

c_pursta_jet = ROOT.TCanvas('c_pursta_jet', 'Distributions jet', 900, 300)
c_pursta_jet.Divide(3,1)
c_pursta_jet.cd(1)
h_pur_eta_jet, h_sta_eta_jet = calc_pursta(h2_jet_resp_eta)
h_pur_eta_jet.Draw('hist')
h_sta_eta_jet.Draw('hist same')
c_pursta_jet.cd(2)
h_pur_phi_jet, h_sta_phi_jet = calc_pursta(h2_jet_resp_phi)
h_pur_phi_jet.Draw('hist')
h_sta_phi_jet.Draw('hist same')
c_pursta_jet.cd(3)
h_pur_pt_jet, h_sta_pt_jet = calc_pursta(h2_jet_resp_pt)
h_pur_pt_jet.Draw('hist')
h_sta_pt_jet.Draw('hist same')
c_pursta_jet.SaveAs('pursta_jet.pdf')
