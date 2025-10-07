import ROOT
import array

nrepeat = 1
#nrepeat = 1000
f = ROOT.TFile('ttbarSel_merged.root')
t = f.Get('tree')
t.SetBranchStatus("*", 0)
t.SetBranchStatus("mcT", 1)
t.SetBranchStatus("mcTbar", 1)
mcT = array.array('f', [0.]*4)
t.SetBranchAddress('mcT', mcT)
mcTbar = array.array('f', [0.]*4)
t.SetBranchAddress('mcTbar', mcTbar)
n = t.GetEntries()
mean_x = 0
mean_y = 0
mean_xx = 0
mean_yy = 0
mean_xy = 0
for irepeat in range(nrepeat):
  for i in range(n):
    t.GetEntry(i)
    #print(mcT, mcTbar)
    mean_x += mcT[0]
    mean_y += mcTbar[0]
    mean_xx += mcT[0]**2
    mean_yy += mcTbar[0]**2
    mean_xy += mcT[0]*mcTbar[0]
mean_x /= n
mean_y /= n
mean_xx /= n
mean_yy /= n
mean_xy /= n
sigma_x = (mean_xx - mean_x**2)**0.5
sigma_y = (mean_yy - mean_y**2)**0.5
cor = (mean_xy - mean_x * mean_y) / (sigma_x * sigma_y)
print(f'cor = {cor}')
