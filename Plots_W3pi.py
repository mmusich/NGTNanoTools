import ROOT, os, math, sys, glob
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import mplhep as hep
plt.style.use(hep.style.CMS)

files = ['/data/evernazz/2025_05_15/CMSSW_15_1_0_pre3/src/WTo3Pi/NoL1Filter/Phase2_HLT_inNANOAODSIM.root']

dataframe_files = ROOT.vector(str)()
for f in files:
    file = ROOT.TFile.Open(f)
    if file.Get("Events"):
        dataframe_files.push_back(f)
df = ROOT.RDataFrame("Events", dataframe_files)

os.system("mkdir -p Plots_W3pi")

#####################################################
# Print summary of the event, including hard scatter particles and W boson decay products

ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    int event_summary(Vint GenPart_pdgId, Vint GenPart_status, Vint GenPart_genPartIdxMother, Vfloat GenPart_mass, Vfloat GenPart_pt) {
        std::cout << std::endl;
        std::cout << "=== Event Summary ===" << std::endl;
        // W boson decay
        for (size_t i = 0; i < GenPart_pdgId.size(); ++i) {
            if (std::abs(GenPart_pdgId[i]) == 24) {
                std::cout << "W boson found (idx = " << i << ", pdgId = " << GenPart_pdgId[i] 
                    <<  ", GenPart_status = " << GenPart_status[i] 
                    << ", GenPart_mass = " << GenPart_mass[i] << ")" 
                    <<  ", GenPart_pt = " << GenPart_pt[i] 
                    << std::endl;
                std::cout << "  -> Decay products:" << std::endl;
                for (size_t j = 0; j < GenPart_pdgId.size(); ++j) {
                    if (GenPart_genPartIdxMother[j] == i) {
                        std::cout << "     idx " << j << ": pdgId = " << GenPart_pdgId[j] << ", GenPart_status = " << GenPart_status[j] << std::endl;
                    }
                }
            }
        }
        std::cout << "=====================" << std::endl;
        return 1;
    }
""")

df = df.Define("event_summary", "event_summary(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother, GenPart_mass, GenPart_pt)")
print(df.AsNumpy(["event_summary"]))

# sys.exit(0)  # Exit after printing the event summary

### Get 3 pions from W decay
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;
    ROOT::RVec<int> get_gen_pion_idices(Vfloat GenPart_pdgId, Vfloat GenPart_status, Vint GenPart_genPartIdxMother) {
        int i_gen_pi1 = -1;
        int i_gen_pi2 = -1;
        int i_gen_pi3 = -1;
        for (size_t i_gen = 0; i_gen < GenPart_pdgId.size(); i_gen ++) {
            if ((GenPart_genPartIdxMother.at(i_gen) >= 0) && 
                    (std::abs(GenPart_pdgId.at(GenPart_genPartIdxMother.at(i_gen))) == 24) && 
                    GenPart_status.at(GenPart_genPartIdxMother.at(i_gen)) == 62) {
                if (std::abs(GenPart_pdgId.at(i_gen)) == 211) {
                    if (i_gen_pi1 == -1) {
                        i_gen_pi1 = i_gen;
                    } else if (i_gen_pi2 == -1) {
                        i_gen_pi2 = i_gen;
                    } else if (i_gen_pi3 == -1) {
                        i_gen_pi3 = i_gen;
                    } else {
                        std::cout << "More than three pions found, only first three will be used." << std::endl;
                        break;
                    }
                }
            }
        }
        return {i_gen_pi1, i_gen_pi2, i_gen_pi3};
    }
""")

### Compute W mass from 3 pions
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;
    float get_W_mass (Vint pion_indices, Vfloat Part_pt, Vfloat Part_eta, Vfloat Part_phi, Vfloat Part_mass) {
        if (pion_indices[0] == -1 || pion_indices[1] == -1 || pion_indices[2] == -1) return 0;
        auto pion0_tlv = TLorentzVector();
        pion0_tlv.SetPtEtaPhiM(Part_pt.at(pion_indices[0]), Part_eta.at(pion_indices[0]), 
            Part_phi.at(pion_indices[0]), Part_mass.at(pion_indices[0]));
        auto pion1_tlv = TLorentzVector();
        pion1_tlv.SetPtEtaPhiM(Part_pt.at(pion_indices[1]), Part_eta.at(pion_indices[1]), 
            Part_phi.at(pion_indices[1]), Part_mass.at(pion_indices[1]));
        auto pion2_tlv = TLorentzVector();
        pion2_tlv.SetPtEtaPhiM(Part_pt.at(pion_indices[2]), Part_eta.at(pion_indices[2]), 
            Part_phi.at(pion_indices[2]), Part_mass.at(pion_indices[2]));
        auto W_tlv = pion0_tlv + pion1_tlv + pion2_tlv;
        return W_tlv.M();
    }
""")

df = df.Define("gen_pion_indices", "get_gen_pion_idices(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") 
df = df.Define("gen_pion_idx_0", "gen_pion_indices[0]").Define("gen_pion_idx_1", "gen_pion_indices[1]").Define("gen_pion_idx_2", "gen_pion_indices[2]")
df = df.Define("gen_W_mass", "get_W_mass(gen_pion_indices, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")

### Get 3 pions from hlt PF candidates
# ROOT.gInterpreter.Declare("""
#     using Vfloat = const ROOT::RVec<float>&;
#     using Vint   = const ROOT::RVec<int>&;
#     ROOT::RVec<int> get_hlt_PFcandidate_idices(Vfloat hltPFCandidate_pdgId) {
#         int i_hlt_pi1 = -1;
#         int i_hlt_pi2 = -1;
#         int i_hlt_pi3 = -1;
#         for (size_t i_hlt = 0; i_hlt < hltPFCandidate_pdgId.size(); i_hlt ++) {
#             if (std::abs(hltPFCandidate_pdgId.at(i_hlt)) == 211) {
#                 if (i_hlt_pi1 == -1) {
#                     i_hlt_pi1 = i_hlt;
#                 } else if (i_hlt_pi2 == -1) {
#                     i_hlt_pi2 = i_hlt;
#                 } else if (i_hlt_pi3 == -1) {
#                     i_hlt_pi3 = i_hlt;
#                 } else {
#                     std::cout << "More than three pions found, only first three will be used." << std::endl;
#                     break;
#                 }
#             }
#         }
#         return {i_hlt_pi1, i_hlt_pi2, i_hlt_pi3};
#     }
# """)
# df = df.Define("hlt_pion_indices", "get_hlt_PFcandidate_idices(hltPFCandidate_pdgId)") 

ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;
    ROOT::RVec<int> get_best_hlt_pion_triplet(Vfloat hltPFCandidate_pdgId, Vfloat hltPFCandidate_pt,
                                              Vfloat hltPFCandidate_eta, Vfloat hltPFCandidate_phi,
                                              Vfloat hltPFCandidate_mass) {

        std::vector<int> pion_indices;
        for (size_t i = 0; i < hltPFCandidate_pdgId.size(); ++i) {
            if (std::abs(hltPFCandidate_pdgId[i]) == 211) {
                pion_indices.push_back(i);
            }
        }

        double target_mass = 80.4;  // W boson mass in GeV
        double closest_mass_diff = std::numeric_limits<double>::max();
        ROOT::RVec<int> best_triplet = {-1, -1, -1};

        if (pion_indices.size() < 3)
            return best_triplet;

        std::cout << "Found " << pion_indices.size() << " pions in HLT candidates" << std::endl;
        float min_pt = 5.0;
        int count = 0;
        for (size_t i = 0; i < pion_indices.size(); ++i) {
            if (hltPFCandidate_pt[pion_indices[i]] < min_pt) continue;
            count++;
            for (size_t j = i + 1; j < pion_indices.size(); ++j) {
                if (hltPFCandidate_pt[pion_indices[j]] < min_pt) continue;
                for (size_t k = j + 1; k < pion_indices.size(); ++k) {
                    if (hltPFCandidate_pt[pion_indices[j]] < min_pt) continue;

                    int idx1 = pion_indices[i];
                    int idx2 = pion_indices[j];
                    int idx3 = pion_indices[k];

                    TLorentzVector p1, p2, p3;
                    p1.SetPtEtaPhiM(hltPFCandidate_pt[idx1], hltPFCandidate_eta[idx1],
                                    hltPFCandidate_phi[idx1], hltPFCandidate_mass[idx1]);
                    p2.SetPtEtaPhiM(hltPFCandidate_pt[idx2], hltPFCandidate_eta[idx2],
                                    hltPFCandidate_phi[idx2], hltPFCandidate_mass[idx2]);
                    p3.SetPtEtaPhiM(hltPFCandidate_pt[idx3], hltPFCandidate_eta[idx3],
                                    hltPFCandidate_phi[idx3], hltPFCandidate_mass[idx3]);

                    TLorentzVector total = p1 + p2 + p3;
                    double mass_diff = std::abs(total.M() - target_mass);

                    if (mass_diff < closest_mass_diff) {
                        closest_mass_diff = mass_diff;
                        best_triplet = {idx1, idx2, idx3};
                    }
                }
            }
        }

        std::cout << "Found " << count << " pions with pt > " << min_pt << std::endl;
        return best_triplet;
    }
""")

df = df.Define("hlt_pion_indices", 
               "get_best_hlt_pion_triplet(hltPFCandidate_pdgId, hltPFCandidate_pt, hltPFCandidate_eta, hltPFCandidate_phi, hltPFCandidate_mass)")
df = df.Define("hlt_pion_idx_0", "hlt_pion_indices[0]").Define("hlt_pion_idx_1", "hlt_pion_indices[1]").Define("hlt_pion_idx_2", "hlt_pion_indices[2]")
df = df.Define("hlt_W_mass", "get_W_mass(hlt_pion_indices, hltPFCandidate_pt, hltPFCandidate_eta, hltPFCandidate_phi, hltPFCandidate_mass)")

h_gen = df.Histo1D("gen_W_mass")
c_gen = np.array([h_gen.GetBinContent(i) for i in range(1, h_gen.GetNbinsX()+1)])
e_gen = np.array([h_gen.GetBinLowEdge(i) for i in range(1, h_gen.GetNbinsX()+2)])

h_hlt = df.Histo1D("hlt_W_mass")
c_hlt = np.array([h_hlt.GetBinContent(i) for i in range(1, h_hlt.GetNbinsX()+1)])
e_hlt = np.array([h_hlt.GetBinLowEdge(i) for i in range(1, h_hlt.GetNbinsX()+2)])

fig, ax = plt.subplots(figsize=(10, 10))
hep.histplot(c_gen, bins=e_gen, label="Gen-level", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_hlt, bins=e_hlt, label="HLT-level", ax=ax, histtype='step', linewidth=2)
ax.set_xlabel("W mass [GeV]")
ax.set_ylabel("Number of events")
ax.legend()
plt.grid()
plt.savefig("Plots_W3pi/W_mass_80.png")
plt.savefig("Plots_W3pi/W_mass_80.pdf")

df = df.Define("gen_pion_pt0", "GenPart_pt[gen_pion_idx_0]").Define("gen_pion_pt1", "GenPart_pt[gen_pion_idx_1]").Define("gen_pion_pt2", "GenPart_pt[gen_pion_idx_2]") 

h_gen_pt0 = df.Histo1D("gen_pion_pt0")
h_gen_pt1 = df.Histo1D("gen_pion_pt1")
h_gen_pt2 = df.Histo1D("gen_pion_pt2")
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
plt.savefig("Plots_W3pi/gen_pion_pt.png")
plt.savefig("Plots_W3pi/gen_pion_pt.pdf")
