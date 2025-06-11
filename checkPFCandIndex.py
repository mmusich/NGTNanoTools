import uproot
import awkward as ak
import numpy as np

# Load NanoAOD file
file = uproot.open("your_custom_nanoaod.root")
tree = file["Events"]

# Load PF candidate and track branches
branches = [
    "hltPFCandidate_pt", "hltPFCandidate_eta", "hltPFCandidate_phi", "hltPFCandidate_charge", "hltPFCandidate_trackIndex",
    "hltGeneralTrack_pt", "hltGeneralTrack_eta", "hltGeneralTrack_phi", "hltGeneralTrack_charge"
]

arrays = tree.arrays(branches)

# Set tolerances
pt_tol = 0.1   # 10%
eta_tol = 0.01
phi_tol = 0.01

# Loop over events
for i in range(len(arrays["hltPFCandidate_pt"])):
    pf_pts = arrays["hltPFCandidate_pt"][i]
    pf_etas = arrays["hltPFCandidate_eta"][i]
    pf_phis = arrays["hltPFCandidate_phi"][i]
    pf_charges = arrays["hltPFCandidate_charge"][i]
    pf_trkIdx = arrays["hltPFCandidate_trackIndex"][i]

    trk_pts = arrays["hltGeneralTrack_pt"][i]
    trk_etas = arrays["hltGeneralTrack_eta"][i]
    trk_phis = arrays["hltGeneralTrack_phi"][i]
    trk_charges = arrays["hltGeneralTrack_charge"][i]

    for j in range(len(pf_pts)):
        idx = pf_trkIdx[j]
        if idx < 0 or idx >= len(trk_pts):
            print(f"[Event {i}] PF candidate {j}: Invalid track index {idx}")
            continue

        # PF candidate values
        pf_pt, pf_eta, pf_phi, pf_charge = pf_pts[j], pf_etas[j], pf_phis[j], pf_charges[j]

        # Matched track values
        trk_pt, trk_eta, trk_phi, trk_charge = trk_pts[idx], trk_etas[idx], trk_phis[idx], trk_charges[idx]

        # Compare
        pt_match = abs(pf_pt - trk_pt) / pf_pt < pt_tol if pf_pt > 0 else False
        eta_match = abs(pf_eta - trk_eta) < eta_tol
        phi_match = abs(pf_phi - trk_phi) < phi_tol
        charge_match = pf_charge == trk_charge

        if not (pt_match and eta_match and phi_match and charge_match):
            print(f"[Event {i}] PF cand {j} does not match track {idx}")
            print(f"  PF  : pt={pf_pt:.2f}, eta={pf_eta:.3f}, phi={pf_phi:.3f}, q={pf_charge}")
            print(f"  Trk : pt={trk_pt:.2f}, eta={trk_eta:.3f}, phi={trk_phi:.3f}, q={trk_charge}")
