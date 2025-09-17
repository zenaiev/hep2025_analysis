# need to download file MERGED/ntuples-mc/TTJets_TuneZ2_7TeV-madgraph-tauola/010003/ttbarSel_merged.root
# via this link: https://cernbox.cern.ch/s/UmbXF1XxVrT4whQ

import ROOT
import time

f = ROOT.TFile.Open('ttbarSel_merged.root')
tree = f.Get('tree')

# top quark px
if 0:
  tree.Draw("mcT[0]>>(100,-500.,500.)")

# top quark px, py, pz, m
if 0:
  c = ROOT.TCanvas('c', 'top quark kinematics', 900, 900)
  c.Divide(2,2)
  c.cd(1)
  tree.Draw("mcT[0]>>h0(100,-500.,500.)")
  c.cd(2)
  tree.Draw("mcT[1]>>h1(100,-500.,500.)")
  c.cd(3)
  tree.Draw("mcT[2]>>h2(100,-1000.,1000.)")
  c.cd(4)
  tree.Draw("mcT[3]>>h3(100,172.,173.)")
  c.SaveAs('c.pdf')

# efficiency calculation
if 0:
  #c = ROOT.TCanvas('c', 'top quark kinematics', 900, 900)
  #c.cd(1)
  cut = "sqrt(pow(mcLp[0],2.)+pow(mcLp[1],2.))>20&&sqrt(pow(mcLm[0],2.)+pow(mcLm[1],2.))>20"
  p = "sqrt(pow(mcLp[0],2.)+pow(mcLp[1],2.)+pow(mcLp[2],2.))"
  cut += f"&&0.5*fabs(log(({p}+mcLp[2])/({p}-mcLp[2])))<2.4"
  p = "sqrt(pow(mcLm[0],2.)+pow(mcLm[1],2.)+pow(mcLm[2],2.))"
  cut += f"&&0.5*fabs(log(({p}+mcLm[2])/({p}-mcLm[2])))<2.4"
  tree.Draw("mcT[0]>>h0(100,-500.,500.)", cut, 'goff')
  #c.SaveAs('c.pdf')
  h0 = ROOT.gDirectory.Get('h0')
  print(f'all events = {tree.GetEntries()}, selected = {h0.GetEntries()}, eff = {h0.GetEntries()/tree.GetEntries():.2f}')

# efficiency calculation using RDataFrame
if 1:
  ROOT.EnableImplicitMT()
  rdf = ROOT.RDataFrame('tree', 'ttbarSel_merged.root')
  rdf = rdf.Define('mcLp_pt', 'sqrt(pow(mcLp[0],2.)+pow(mcLp[1],2.))')
  rdf = rdf.Define('mcLp_p', 'sqrt(pow(mcLp_pt,2.)+pow(mcLp[2],2.))')
  rdf = rdf.Define('mcLp_eta', '0.5*log((mcLp_p+mcLp[2])/(mcLp_p-mcLp[2]))')
  rdf = rdf.Define('mcLm_pt', 'sqrt(pow(mcLm[0],2.)+pow(mcLm[1],2.))')
  rdf = rdf.Define('mcLm_p', 'sqrt(pow(mcLm_pt,2.)+pow(mcLm[2],2.))')
  rdf = rdf.Define('mcLm_eta', '0.5*log((mcLm_p+mcLm[2])/(mcLm_p-mcLm[2]))')
  rdf_selected = rdf.Filter('mcLp_pt>20&&fabs(mcLp_eta)<2.4&&mcLm_pt>20&&fabs(mcLm_eta)<2.4')
  nevents = rdf.Count()
  nevents_selected = rdf_selected.Count()
  print(f'all events = {nevents.GetValue()}, selected = {nevents_selected.GetValue()}, eff = {nevents_selected.GetValue()/nevents.GetValue():.2f}')
  print(f'number of runs: rdf = {rdf.GetNRuns()}, rdf = {rdf_selected.GetNRuns()}')


#time.sleep(10000)
