import ROOT, os, math, sys, glob
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import mplhep as hep
plt.style.use(hep.style.CMS)

from Functions_W3pi import *

#################################################################
# this version plots the W mass from gen level and from the HLT
# PF candidates (the 3 pions matched to gen level pions)
#################################################################

files = ['/data/evernazz/2025_05_15/CMSSW_14_0_9/src/PrivateProduction_1000ev/HLT_inNANOAODSIM.root']

dataframe_files = ROOT.vector(str)()
for f in files:
    file = ROOT.TFile.Open(f)
    if file.Get("Events"):
        dataframe_files.push_back(f)
df = ROOT.RDataFrame("Events", dataframe_files)#.Range(0, 10) ############# [FIXME] 

output_dir = "Plots_W3pi_2_1000ev"
os.system(f"mkdir -p {output_dir}")

#####################################################
# Plot all pions at gen level

df = df.Define("gen_pion_pt", "GenPart_pt[GenPart_pdgId == 211 || GenPart_pdgId == -211]")
df = df.Define("gen_pion_eta", "GenPart_eta[GenPart_pdgId == 211 || GenPart_pdgId == -211]")

df = df.Define("hlt_pion_pt", "hltPFCandidate_pt[hltPFCandidate_pdgId == 211 || hltPFCandidate_pdgId == -211]")
df = df.Define("hlt_pion_eta", "hltPFCandidate_eta[hltPFCandidate_pdgId == 211 || hltPFCandidate_pdgId == -211]")

h_gen_pt = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 100, 0, 200), "gen_pion_pt")
c_gen_pt = np.array([h_gen_pt.GetBinContent(i)/h_gen_pt.Integral() for i in range(1, h_gen_pt.GetNbinsX()+1)])
e_gen_pt = np.array([h_gen_pt.GetBinLowEdge(i) for i in range(1, h_gen_pt.GetNbinsX()+2)])
h_hlt_pt = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 100, 0, 200), "hlt_pion_pt")
c_hlt_pt = np.array([h_hlt_pt.GetBinContent(i)/h_hlt_pt.Integral() for i in range(1, h_hlt_pt.GetNbinsX()+1)])
e_hlt_pt = np.array([h_hlt_pt.GetBinLowEdge(i) for i in range(1, h_hlt_pt.GetNbinsX()+2)])
fig, ax = plt.subplots(figsize=(10, 10))
hep.histplot(c_gen_pt, bins=e_gen_pt, label=f"All gen-level pions ({int(h_gen_pt.Integral())})", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_hlt_pt, bins=e_hlt_pt, label=f"All hlt-level pions ({int(h_hlt_pt.Integral())})", ax=ax, histtype='step', linewidth=2)
ax.set_xlabel("$p_T$ [GeV]")
ax.set_ylabel("Events")
ax.legend(loc='upper right', fontsize=16)
plt.grid()
print(f" ### INFO: Saving plot {output_dir}/All_pions_pt.png")
plt.savefig(f"{output_dir}/All_pions_pt.png")
plt.savefig(f"{output_dir}/All_pions_pt.pdf")

h_gen_eta = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 100, -6, 6), "gen_pion_eta")
c_gen_eta = np.array([h_gen_eta.GetBinContent(i)/h_gen_eta.Integral() for i in range(1, h_gen_eta.GetNbinsX()+1)])
e_gen_eta = np.array([h_gen_eta.GetBinLowEdge(i) for i in range(1, h_gen_eta.GetNbinsX()+2)])
h_hlt_eta = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 100, -6, 6), "hlt_pion_eta")
c_hlt_eta = np.array([h_hlt_eta.GetBinContent(i)/h_hlt_eta.Integral()  for i in range(1, h_hlt_eta.GetNbinsX()+1)])
e_hlt_eta = np.array([h_hlt_eta.GetBinLowEdge(i) for i in range(1, h_hlt_eta.GetNbinsX()+2)])
fig, ax = plt.subplots(figsize=(10, 10))
hep.histplot(c_gen_eta, bins=e_gen_eta, label=f"All gen-level pions ({int(h_gen_eta.Integral())})", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_hlt_eta, bins=e_hlt_eta, label=f"All hlt-level pions ({int(h_hlt_eta.Integral())})", ax=ax, histtype='step', linewidth=2)
ax.set_xlabel("$\eta$")
ax.set_ylabel("Events")
ax.legend(loc='upper left', fontsize=16)
plt.grid()
print(f" ### INFO: Saving plot {output_dir}/All_pions_eta.png")
plt.savefig(f"{output_dir}/All_pions_eta.png")
plt.savefig(f"{output_dir}/All_pions_eta.pdf")

#####################################################

# ### Debug
# df = df.Define("event_summary", 
#         "event_summary(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother, GenPart_mass, GenPart_pt, GenPart_eta)")
# print(df.AsNumpy(["event_summary"]))

#####################################################

### Get 3 pions from W decay
df = df.Define("gen_pion_indices", "get_gen_pion_idices(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") 
df = df.Define("gen_pion_idx_0", "gen_pion_indices[0]").Define("gen_pion_idx_1", "gen_pion_indices[1]").Define("gen_pion_idx_2", "gen_pion_indices[2]")
df = df.Define("gen_W_mass", "get_W_mass(gen_pion_indices, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
h_gen = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 100, 1, 150), "gen_W_mass")
c_gen = np.array([h_gen.GetBinContent(i) for i in range(1, h_gen.GetNbinsX()+1)])
e_gen = np.array([h_gen.GetBinLowEdge(i) for i in range(1, h_gen.GetNbinsX()+2)])

### Filter events where one pion has eta > 4
df = df.Filter("std::abs(GenPart_eta[gen_pion_idx_0]) < 4 && std::abs(GenPart_eta[gen_pion_idx_1]) < 4 && std::abs(GenPart_eta[gen_pion_idx_2]) < 4")
h_gen_filter = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 100, 1, 150), "gen_W_mass")
c_gen_filter = np.array([h_gen_filter.GetBinContent(i) for i in range(1, h_gen_filter.GetNbinsX()+1)])
e_gen_filter = np.array([h_gen_filter.GetBinLowEdge(i) for i in range(1, h_gen_filter.GetNbinsX()+2)])

### For each gen pion, among the hlt PF candidates with dR < 0.4 and dPt < 10 GeV, find the closest one in dR
df = df.Define("hlt_pion_indices",
               "match_gen_to_hlt(gen_pion_indices, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, " \
               "hltPFCandidate_pdgId, hltPFCandidate_pt, hltPFCandidate_eta, hltPFCandidate_phi, hltPFCandidate_mass)")
df = df.Define("hlt_pion_idx_0", "hlt_pion_indices[0]").Define("hlt_pion_idx_1", "hlt_pion_indices[1]").Define("hlt_pion_idx_2", "hlt_pion_indices[2]")
df = df.Define("hlt_W_mass", "get_W_mass(hlt_pion_indices, hltPFCandidate_pt, hltPFCandidate_eta, hltPFCandidate_phi, hltPFCandidate_mass)")
h_hlt = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 100, 1, 150), "hlt_W_mass")
c_hlt = np.array([h_hlt.GetBinContent(i) for i in range(1, h_hlt.GetNbinsX()+1)])
e_hlt = np.array([h_hlt.GetBinLowEdge(i) for i in range(1, h_hlt.GetNbinsX()+2)])

fig, ax = plt.subplots(figsize=(10, 10))
hep.histplot(c_gen, bins=e_gen, label="Gen-level", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_gen_filter, bins=e_gen_filter, label="Gen-level within acceptance", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_hlt, bins=e_hlt, label="HLT-level", ax=ax, histtype='step', linewidth=2)
ax.set_xlabel("W mass [GeV]")
ax.set_ylabel("Number of events")
ax.legend(loc='upper left', fontsize=16)
plt.grid()
print(f" ### INFO: Saving plot {output_dir}/W_mass.png")
plt.savefig(f"{output_dir}/W_mass.png")
plt.savefig(f"{output_dir}/W_mass.pdf")

#####################################################

### Plot gen pions pt (not ordered in pt, so not really meaningful)

df = df.Define("gen_pion_pt0", "GenPart_pt[gen_pion_idx_0]").Define("gen_pion_pt1", "GenPart_pt[gen_pion_idx_1]").Define("gen_pion_pt2", "GenPart_pt[gen_pion_idx_2]") 
h_gen_pt0 = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 100, 0, 150), "gen_pion_pt0")
h_gen_pt1 = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 100, 0, 150), "gen_pion_pt1")
h_gen_pt2 = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 100, 0, 150), "gen_pion_pt2")
c_gen_pt0 = np.array([h_gen_pt0.GetBinContent(i) for i in range(1, h_gen_pt0.GetNbinsX()+1)])
c_gen_pt1 = np.array([h_gen_pt1.GetBinContent(i) for i in range(1, h_gen_pt1.GetNbinsX()+1)])
c_gen_pt2 = np.array([h_gen_pt2.GetBinContent(i) for i in range(1, h_gen_pt2.GetNbinsX()+1)])
e_gen_pt0 = np.array([h_gen_pt0.GetBinLowEdge(i) for i in range(1, h_gen_pt0.GetNbinsX()+2)])
e_gen_pt1 = np.array([h_gen_pt1.GetBinLowEdge(i) for i in range(1, h_gen_pt1.GetNbinsX()+2)])
e_gen_pt2 = np.array([h_gen_pt2.GetBinLowEdge(i) for i in range(1, h_gen_pt2.GetNbinsX()+2)])

fig, ax = plt.subplots(figsize=(10, 10))
hep.histplot(c_gen_pt0, bins=e_gen_pt0, label="Gen pion 0", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_gen_pt1, bins=e_gen_pt1, label="Gen pion 1", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_gen_pt2, bins=e_gen_pt2, label="Gen pion 2", ax=ax, histtype='step', linewidth=2)
ax.set_xlabel("Gen pion pt [GeV]")
ax.set_ylabel("Number of events")
ax.legend()
plt.grid()
print(f" ### INFO: Saving plot {output_dir}/gen_pion_pt.png")
plt.savefig(f"{output_dir}/gen_pion_pt.png")
plt.savefig(f"{output_dir}/gen_pion_pt.pdf")
