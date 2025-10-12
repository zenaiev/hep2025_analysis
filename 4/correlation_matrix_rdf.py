import ROOT
#ROOT.EnableImplicitMT(3)

rdf = ROOT.RDataFrame('tree', 'ttbarSel_merged.root')
#rdf = ROOT.RDataFrame('tree', ['ttbarSel_merged.root']*1000)

#variables = ['mcT[0]', 'mcT[1]', 'mcT[2]']
#variables += ['mcTbar[0]', 'mcTbar[1]', 'mcTbar[2]']
variables = ['(mcNu[0]+mcNubar[0])', '(mcNu[1]+mcNubar[1])', '(mcNu[2]+mcNubar[2])']
#variables += ['metPx', 'metPy']
variables += ['(mcLp[0]+mcLm[0])', '(mcLp[1]+mcLm[1])', '(mcLp[2]+mcLm[2])']
variables += ['(mcB[0]+mcBbar[0])', '(mcB[1]+mcBbar[1])', '(mcB[2]+mcBbar[2])']
means = [0] * len(variables)
covs = [[0 for _ in range(len(variables))] for _ in range(len(variables))]
# book mean values
for ivariable in range(len(variables)):
  means[ivariable] = rdf.Define('TMP', variables[ivariable]).Mean('TMP')
  for ivariable2 in range(ivariable, len(variables)):
    covs[ivariable][ivariable2] = rdf.Define('TMP', f'{variables[ivariable]}*{variables[ivariable2]}').Mean('TMP')
    covs[ivariable2][ivariable] = covs[ivariable][ivariable2]
# calculate mean value
for ivariable in range(len(variables)):
  means[ivariable] = means[ivariable].GetValue()
  for ivariable2 in range(len(variables)):
    covs[ivariable][ivariable2] = covs[ivariable][ivariable2].GetValue()
# calculate correlation matrix
cors = [[0 for _ in range(len(variables))] for _ in range(len(variables))]
for ivariable in range(len(variables)):
  print((covs[ivariable][ivariable] - means[ivariable]**2)**0.5)
  for ivariable2 in range(len(variables)):
    numer = covs[ivariable][ivariable2] - means[ivariable] * means[ivariable2]
    denom = ((covs[ivariable][ivariable] - means[ivariable]**2) * (covs[ivariable2][ivariable2] - means[ivariable2]**2))**0.5
    cors[ivariable][ivariable2] = numer / denom
print(means)
print(covs)
print(cors)
offset = 14
def make_title(name):
  return name.replace(')','').replace('(','').replace('mc','')
print(f'{"":{offset}s}', end='')
for i in range(len(variables)):
  print(f'{make_title(variables[i]):>{offset}s}', end='')
print()
for i in range(len(variables)):
  print(f'{make_title(variables[i]):{offset}s}', end='')
  for ii in range(len(variables)):
    if abs(cors[i][ii]) > 0.05:
      print(f'{cors[i][ii]:{offset}.2f}', end='')
    else:
      print(f'{" ":{offset}s}', end='')
  print()

print(f'number of runs: rdf = {rdf.GetNRuns()}')
