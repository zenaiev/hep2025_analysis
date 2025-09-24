import ROOT
import time

f = ROOT.TFile.Open('ttbarSel_merged.root')
tree = f.Get('tree')

# list to keep objects
_keep = []

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

# create data frame
#ROOT.EnableImplicitMT()
rdf = ROOT.RDataFrame('tree', 'ttbarSel_merged.root')

# efficiency calculation using RDataFrame
if 0:
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
  #print(f'all events = {rdf.Count().GetValue()}, selected = {rdf_selected.Count().GetValue()}, eff = {rdf_selected.Count().GetValue()/rdf.Count().GetValue():.2f}')
  print(f'number of runs: rdf = {rdf.GetNRuns()}, rdf_selected = {rdf_selected.GetNRuns()}')
  if 0:
    rdf_selected.Snapshot('tree', 'outttbar.root', ['mcLp_pt', 'mcLp_eta', 'mcLm_pt', 'mcLm_eta'])
  if 1:
    ROOT.RDF.SaveGraph(rdf_selected, "rdf_selected.dot")

# plot top quark momentum components using RDataFrame::Histo1D()
# using little functions to void repeating the same code
if 0:
  def draw_one_component(icanvas, icomponent, title, xmin=-500., xmax=500.):
    c.cd(icanvas)
    h = rdf.Define('var', f"mcT[{icomponent}]").Histo1D(ROOT.RDF.TH1DModel("h", title, 100, xmin, xmax), "var")
    h.Draw()
    _keep.append(h)
  def draw_one_component_kwargs(icanvas, icomponent, title, **kwargs):
    xmin = kwargs.get('xmin', -500.)
    xmax = kwargs.get('xmax', 500.)
    ylog = kwargs.get('ylog')
    c.cd(icanvas)
    h = rdf.Define('var', f"mcT[{icomponent}]").Histo1D(ROOT.RDF.TH1DModel("h", title, 100, xmin, xmax), "var")
    if ylog is not None:
      ROOT.gPad.SetLogy()
    h.Draw()
    _keep.append(h)
  c = ROOT.TCanvas('c', 'top quark kinematics', 900, 900)
  c.Divide(2,2)
  draw_one_component(1, 0, 'top px')
  c.cd(2)
  draw_one_component(2, 1, 'top py')
  c.cd(3)
  #draw_one_component(3, 2, 'top pz', xmin=-1000., xmax=1000.)
  draw_one_component_kwargs(3, 2, 'top pz', xmin=-1000., xmax=1000.)
  c.cd(4)
  draw_one_component_kwargs(4, 3, 'top m', ylog=True)
  c.SaveAs('c.pdf')

# plot correlations of top and antitop quark momenta using RDataFrame::Graph()
if 1:
  c = ROOT.TCanvas('c', 'top quark kinematics', 900, 900)
  c.Divide(2,2)
  c.cd(1)
  g = rdf.Define("x", "mcT[0]").Define("y", "mcTbar[0]").Graph("x", "y")
  g.Draw("AP")
  _keep.append(g)
  c.cd(2)
  g = rdf.Define("x", "mcT[1]").Define("y", "mcTbar[1]").Graph("x", "y")
  g.Draw("AP")
  _keep.append(g)
  c.cd(3)
  g = rdf.Define("x", "mcT[2]").Define("y", "mcTbar[2]").Graph("x", "y")
  g.Draw("AP")
  _keep.append(g)
  c.SaveAs('c.pdf')

#time.sleep(10000)
