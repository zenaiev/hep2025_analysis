import ROOT
import array

def add_branch(t, branch_name, branch_type, branch_size):
  br = array.array(branch_type, [0]*branch_size)
  t.SetBranchStatus(branch_name, 1)
  t.SetBranchAddress(branch_name, br)
  return br

#ROOT.gStyle.SetOptStat(0)
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
        else:
          g_jet_eta.AddPoint(genVector.Eta(), rec_eta[best_rec])
          g_jet_phi.AddPoint(genVector.Phi(), rec_phi[best_rec])
          g_jet_pt.AddPoint(genVector.Pt(), abs(rec_pt[best_rec]))
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
    
