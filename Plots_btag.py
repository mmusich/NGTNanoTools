import ROOT, os, math, sys, glob
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pdb
from sklearn.metrics import roc_curve, auc
import mplhep as hep
plt.style.use(hep.style.CMS)

from Functions_btag import *

# files = ['/eos/user/e/evernazz/www/NGT/NanoAOD/TTbar/Phase2_HLT_inNANOAODSIM.root']
files = ['/data/evernazz/2025_05_15/CMSSW_15_1_0_pre3/src/TTbar/Phase2_HLT_inNANOAODSIM.root']

dataframe_files = ROOT.vector(str)()
for f in files:
    file = ROOT.TFile.Open(f)
    if file.Get("Events"):
        dataframe_files.push_back(f)
df = ROOT.RDataFrame("Events", dataframe_files)

output_dir = "Plots_btag_2"
os.system(f"mkdir -p {output_dir}")

#####################################################
### Compute b tagging scores vs all for all HLT jets

df = df.Define("hltAK4PuppiJet_DeepFlavour_prob_b_vs_all", 
    "(hltAK4PuppiJet_DeepFlavour_prob_b + hltAK4PuppiJet_DeepFalvour_prob_bb + hltAK4PuppiJet_DeepFalvour_prob_lepb) / " \
    "(hltAK4PuppiJet_DeepFlavour_prob_b + hltAK4PuppiJet_DeepFalvour_prob_bb + hltAK4PuppiJet_DeepFalvour_prob_lepb + " \
    "hltAK4PuppiJet_DeepFalvour_prob_c + hltAK4PuppiJet_DeepFalvour_prob_g + hltAK4PuppiJet_DeepFalvour_prob_uds)")

#####################################################
### For each hlt jet, find index of the matched gen jet 

df = df.Define("matched_gen_jet_indices",
               "match_jets_hlt_to_gen(GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, " \
               "hltAK4PuppiJet_pt, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi, hltAK4PuppiJet_mass)")

df = df.Define("hltAK4PuppiJet_mask_matched", "(matched_gen_jet_indices != -1)")
df = df.Define("hltAK4PuppiJet_mask_pt_eta", "(hltAK4PuppiJet_pt > 20) && (hltAK4PuppiJet_eta > -2.5) && (hltAK4PuppiJet_eta < 2.5)")

#####################################################
### Plot b tagging scores for HLT jets matched to gen jets
for score in ['DeepFlavour_prob_b', 'DeepFalvour_prob_bb', 'DeepFalvour_prob_lepb', \
    'DeepFalvour_prob_c', 'DeepFalvour_prob_g', 'DeepFalvour_prob_uds', 'DeepFlavour_prob_b_vs_all']:

    h_all = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 50, 0.0, 1.0), f"hltAK4PuppiJet_{score}")
    c_all = np.array([h_all.GetBinContent(i) for i in range(1, h_all.GetNbinsX()+1)])
    e_all = np.array([h_all.GetBinLowEdge(i) for i in range(1, h_all.GetNbinsX()+2)])

    df = df.Define(f"hltAK4PuppiJet_{score}_matched", f"hltAK4PuppiJet_{score}[hltAK4PuppiJet_mask_matched]")
    h_matched = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 50, 0.0, 1.0), f"hltAK4PuppiJet_{score}_matched")
    c_matched = np.array([h_matched.GetBinContent(i) for i in range(1, h_matched.GetNbinsX()+1)])
    e_matched = np.array([h_matched.GetBinLowEdge(i) for i in range(1, h_matched.GetNbinsX()+2)])

    df = df.Define(f"hltAK4PuppiJet_{score}_matched_pt_eta", f"hltAK4PuppiJet_{score}[hltAK4PuppiJet_mask_matched && hltAK4PuppiJet_mask_pt_eta]")
    h_matched_pt_eta = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 50, 0.0, 1.0), f"hltAK4PuppiJet_{score}_matched_pt_eta")
    c_matched_pt_eta = np.array([h_matched_pt_eta.GetBinContent(i) for i in range(1, h_matched_pt_eta.GetNbinsX()+1)])
    e_matched_pt_eta = np.array([h_matched_pt_eta.GetBinLowEdge(i) for i in range(1, h_matched_pt_eta.GetNbinsX()+2)])

    fig, ax = plt.subplots(figsize=(10, 10))
    hep.histplot(c_all, bins=e_all, label="All HLT jets", ax=ax, histtype='step', linewidth=2)
    hep.histplot(c_matched, bins=e_matched, label="Matched HLT jets", ax=ax, histtype='step', linewidth=2)
    hep.histplot(c_matched_pt_eta, bins=e_matched_pt_eta, label="Matched HLT jets $p_T > 20$ GeV, $|\eta|<2.5$", ax=ax, histtype='step', linewidth=2)
    ax.set_xlabel(f"b-tag {score}")
    ax.set_ylabel("Number of Jets")
    ax.legend(fontsize=16)
    plt.grid()
    print(f" ### INFO: Saving plot in {output_dir}/Matched_hlt_jets_{score}.png")
    plt.savefig(f"{output_dir}/Matched_hlt_jets_{score}.png")
    plt.savefig(f"{output_dir}/Matched_hlt_jets_{score}.pdf")

#####################################################
### For each hlt jet, define whether it is matched to a b or non-b gen jet

df = df.Define("hltAK4PuppiJet_mask_matched_to_b",
    "(matched_gen_jet_indices != -1) " \
    "&& (hltAK4PuppiJet_pt > 20) && (hltAK4PuppiJet_eta > -2.5) && (hltAK4PuppiJet_eta < 2.5) " \
    "&& (Take(GenJet_hadronFlavour, matched_gen_jet_indices) == 5)")

# df = df.Define("hltAK4PuppiJet_mask_matched_to_non_b",
#     "(matched_gen_jet_indices != -1) " \
#     "&& (hltAK4PuppiJet_pt > 20) && (hltAK4PuppiJet_eta > -2.5) && (hltAK4PuppiJet_eta < 2.5) " \
#     "&& (Take(GenJet_hadronFlavour, matched_gen_jet_indices) != 5)")

#####################################################
### Plot b tagging scores for HLT jets matched to b or non-b gen jets
for score in ['DeepFlavour_prob_b', 'DeepFalvour_prob_bb', 'DeepFalvour_prob_lepb', \
    'DeepFalvour_prob_c', 'DeepFalvour_prob_g', 'DeepFalvour_prob_uds', 'DeepFlavour_prob_b_vs_all']:
    df = df.Define(f"hltAK4PuppiJet_{score}_matched_to_b", f"hltAK4PuppiJet_{score}[hltAK4PuppiJet_mask_matched_to_b]")
    df = df.Define(f"hltAK4PuppiJet_{score}_matched_to_non_b", f"hltAK4PuppiJet_{score}[!hltAK4PuppiJet_mask_matched_to_b]")
    h_matched = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 50, 0.0, 1.0), f"hltAK4PuppiJet_{score}_matched_to_b")
    c_matched = np.array([h_matched.GetBinContent(i)/h_matched.Integral() for i in range(1, h_matched.GetNbinsX()+1)])
    e_matched = np.array([h_matched.GetBinLowEdge(i) for i in range(1, h_matched.GetNbinsX()+2)])
    h_unmatched = df.Histo1D(ROOT.RDF.TH1DModel("h", "", 50, 0.0, 1.0), f"hltAK4PuppiJet_{score}_matched_to_non_b")
    c_unmatched = np.array([h_unmatched.GetBinContent(i)/h_unmatched.Integral() for i in range(1, h_unmatched.GetNbinsX()+1)])
    e_unmatched = np.array([h_unmatched.GetBinLowEdge(i) for i in range(1, h_unmatched.GetNbinsX()+2)])

    fig, ax = plt.subplots(figsize=(10, 10))
    hep.histplot(c_matched, bins=e_matched, label="HLT jets matched to b-jets", ax=ax, histtype='step', linewidth=2)
    hep.histplot(c_unmatched, bins=e_unmatched, label="HLT jets matched to non-b-jets", ax=ax, histtype='step', linewidth=2)
    ax.set_xlabel(f"b-tag {score}")
    ax.set_ylabel("Number of Jets")
    ax.legend()
    plt.grid()
    print(f" ### INFO: Saving plot in {output_dir}/Matched_b_hlt_jets_{score}.png")
    plt.savefig(f"{output_dir}/Matched_b_hlt_jets_{score}.png")
    plt.savefig(f"{output_dir}/Matched_b_hlt_jets_{score}.pdf")

#####################################################
### Plot ROC curve for b-tag vs all

tpr = []
fpr = []

for threshold in range(1, h_matched.GetNbinsX()+2):
    tp = h_matched.Integral(threshold, h_matched.GetNbinsX()+1)/(h_matched.Integral()) # True Positive Rate
    fp = h_unmatched.Integral(threshold, h_unmatched.GetNbinsX()+1)/(h_unmatched.Integral()) # False Positive Rate
    tpr.append(tp)
    fpr.append(fp)

fig, ax = plt.subplots(figsize=(10, 10))
ax.plot(tpr, fpr, label='ROC Curve DeepFlavour')
ax.plot([0, 1], [0, 1], linestyle='--', color='red', label='Random Guessing')
ax.set_xlabel('Tagging Efficency')
ax.set_ylabel('Mistagging Rate')
ax.set_title('ROC Curve for DeepJet B-Tagging')
ax.legend()
plt.grid()
print(f" ### INFO: Saving plot in {output_dir}/ROC.png")
plt.savefig(f"{output_dir}/ROC.png")
plt.savefig(f"{output_dir}/ROC.pdf")

#####################################################
### Plot ROC curve for b-tag vs all (Uttiya's method)

arrs = df.AsNumpy(["hltAK4PuppiJet_DeepFlavour_prob_b_vs_all", "hltAK4PuppiJet_mask_matched_to_b"])
scores = np.concatenate(arrs["hltAK4PuppiJet_DeepFlavour_prob_b_vs_all"])
labels = np.concatenate(arrs["hltAK4PuppiJet_mask_matched_to_b"]).astype(int)

# Compute ROC curves
fpr_dj, tpr_dj, _ = roc_curve(labels, scores)
auc_dj = auc(fpr_dj, tpr_dj)

# Plot ROC: TPR vs FPR with log y-axis
plt.figure()
plt.plot(tpr_dj, fpr_dj, label=f"DeepJet (AUC = {auc_dj:.3f})", color='blue')
plt.axhline(1, color='gray', linestyle='--', linewidth=0.8)

# Axes and log-scale
plt.xlabel("True Positive Rate")
plt.ylabel("False Positive Rate")
plt.yscale("log")
plt.xlim(0, 1)
plt.ylim(1e-3, 1)
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.legend()

# CMS Phase-II Label
hep.cms.label("Simulation Preliminary", data=True, year="Phase-II", com=14)

plt.tight_layout()
plt.show()

print(f" ### INFO: Saving plot in {output_dir}/ROC_Uttiya.png")
plt.savefig(f"{output_dir}/ROC_Uttiya.png")
plt.savefig(f"{output_dir}/ROC_Uttiya.pdf")