import ROOT
ROOT.EnableImplicitMT(3)

#rdf = ROOT.RDataFrame('tree', 'ttbarSel_merged.root')
rdf = ROOT.RDataFrame('tree', ['ttbarSel_merged.root']*1000)
rdf = rdf.Define('x', 'mcT[0]')
rdf = rdf.Define('y', 'mcTbar[0]')
rdf = rdf.Define('xx', 'x*x')
rdf = rdf.Define('yy', 'y*y')
rdf = rdf.Define('xy', 'x*y')
mean_x = rdf.Mean('x')
mean_y = rdf.Mean('y')
mean_xx = rdf.Mean('xx')
mean_yy = rdf.Mean('yy')
mean_xy = rdf.Mean('xy')
sigma_x = (mean_xx.GetValue() - mean_x.GetValue()**2)**0.5
sigma_y = (mean_yy.GetValue() - mean_y.GetValue()**2)**0.5
#cor = (mean_xy - mean_x) / (sigma_x * sigma_y)
cor = (mean_xy.GetValue() - mean_x.GetValue()) / (sigma_x * sigma_y)
print(f'cor = {cor}')
print(f'number of runs: rdf = {rdf.GetNRuns()}')
