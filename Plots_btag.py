import ROOT, os, math, sys, glob
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import mplhep as hep
plt.style.use(hep.style.CMS)

# files = ['/eos/user/e/evernazz/www/NGT/NanoAOD/TTbar/Phase2_HLT_inNANOAODSIM.root']
files = ['/data/evernazz/2025_05_15/CMSSW_15_1_0_pre3/src/TTbar/Phase2_HLT_inNANOAODSIM.root']

dataframe_files = ROOT.vector(str)()
for f in files:
    file = ROOT.TFile.Open(f)
    if file.Get("Events"):
        dataframe_files.push_back(f)
df = ROOT.RDataFrame("Events", dataframe_files)

output_dir = "Plots_btag"
os.system(f"mkdir -p {output_dir}")

### For each hlt jet, find index of the matched gen jet
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    ROOT::RVec<int> match_jets_hlt_to_gen (Vfloat GenJet_pt, Vfloat GenJet_eta, Vfloat GenJet_phi, Vfloat GenJet_mass,
        Vfloat hltAK4PuppiJet_pt, Vfloat hltAK4PuppiJet_eta, Vfloat hltAK4PuppiJet_phi, Vfloat hltAK4PuppiJet_mass) {

        const double dR_max = 0.4; // dR < dR_max
        const double dPt_max_factor = 0.5; // dPt < dPt_max_factor * pt_resolution
        const double pt_resolution = 5.0;

        ROOT::RVec<int> matched_indices(hltAK4PuppiJet_pt.size(), -1);

        for (size_t i_hlt = 0; i_hlt < hltAK4PuppiJet_pt.size(); i_hlt ++) {
                          
            // Skip jets with pt < 20 GeV ?
            // if (hltAK4PuppiJet_pt.at(i_hlt) < 20.0) continue;
            
            auto hlt_tlv = TLorentzVector();
            hlt_tlv.SetPtEtaPhiM(hltAK4PuppiJet_pt.at(i_hlt), hltAK4PuppiJet_eta.at(i_hlt), hltAK4PuppiJet_phi.at(i_hlt), hltAK4PuppiJet_mass.at(i_hlt));
            double best_dR = std::numeric_limits<double>::infinity();
            int best_match = -1;
            
            for (size_t i_gen = 0; i_gen < GenJet_pt.size(); i_gen ++) {
                auto gen_tlv = TLorentzVector();
                gen_tlv.SetPtEtaPhiM(GenJet_pt.at(i_gen), GenJet_eta.at(i_gen), GenJet_phi.at(i_gen), GenJet_mass.at(i_gen));
                double dR = gen_tlv.DeltaR(hlt_tlv);
                if (dR > best_dR) continue;
                if (dR < dR_max) {
                    double dPt = std::abs(gen_tlv.Pt() - hlt_tlv.Pt());
                    if (dPt > dPt_max_factor * pt_resolution) {
                        best_dR = dR;
                        best_match = static_cast<int>(i_gen);
                    }
                }
            }
            matched_indices[i_hlt] = best_match;
            // std::cout << "HLT Jet index " << i_hlt << ": matched to gen jet index " << best_match << " with dR = " << best_dR << std::endl;
        }
        return matched_indices;
    }
""")

df = df.Define("matched_gen_jet_indices",
               "match_jets_hlt_to_gen(GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, " \
               "hltAK4PuppiJet_pt, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi, hltAK4PuppiJet_mass)") 

# h = df.Histo1D("matched_gen_jet_indices")
# c = np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)])

### For each hlt jet matched to gen b-jet (hadronFlavour = 5), get the b tagging scores
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;
    ROOT::RVec<float> get_btag_scores_for_hlt_bjets (Vint matched_gen_jet_indices, Vint GenJet_hadronFlavour,
            Vfloat hltAK4PuppiJet_DeepFlavour_prob_b, Vfloat hltAK4PuppiJet_DeepFlavour_prob_bb, Vfloat hltAK4PuppiJet_DeepFlavour_prob_lepb,
            Vfloat hltAK4PuppiJet_DeepFlavour_prob_c, Vfloat hltAK4PuppiJet_DeepFlavour_prob_g, Vfloat hltAK4PuppiJet_DeepFlavour_prob_uds) {
        ROOT::RVec<float> btag_scores(matched_gen_jet_indices.size(), -1.0);
        for (size_t i_hlt = 0; i_hlt < matched_gen_jet_indices.size(); ++i_hlt) {
            if ((matched_gen_jet_indices.at(i_hlt) != -1) && std::abs(GenJet_hadronFlavour.at(matched_gen_jet_indices.at(i_hlt))) == 5) {
                float b_vs_all_score = (hltAK4PuppiJet_DeepFlavour_prob_b.at(i_hlt) + hltAK4PuppiJet_DeepFlavour_prob_bb.at(i_hlt) +
                                    hltAK4PuppiJet_DeepFlavour_prob_lepb.at(i_hlt)) /
                                   (hltAK4PuppiJet_DeepFlavour_prob_b.at(i_hlt) + hltAK4PuppiJet_DeepFlavour_prob_bb.at(i_hlt) +
                                    hltAK4PuppiJet_DeepFlavour_prob_lepb.at(i_hlt) + hltAK4PuppiJet_DeepFlavour_prob_c.at(i_hlt) +
                                    hltAK4PuppiJet_DeepFlavour_prob_g.at(i_hlt) + hltAK4PuppiJet_DeepFlavour_prob_uds.at(i_hlt));
                btag_scores[i_hlt] = static_cast<float>(b_vs_all_score);
                // std::cout << "HLT Jet index " << i_hlt << ": matched to gen jet index " << matched_gen_jet_indices.at(i_hlt)
                    // << " with hadronFlavour = " << GenJet_hadronFlavour.at(matched_gen_jet_indices.at(i_hlt))
                    // << " and b-tagging score = " << btag_scores.at(i_hlt) << std::endl;
            }
        }
        return btag_scores;
    }
""")

df = df.Define("btag_scores_signal",
               "get_btag_scores_for_hlt_bjets(matched_gen_jet_indices, GenJet_hadronFlavour, " \
               "hltAK4PuppiJet_DeepFlavour_prob_b, hltAK4PuppiJet_DeepFalvour_prob_bb, hltAK4PuppiJet_DeepFalvour_prob_lepb, " \
               "hltAK4PuppiJet_DeepFalvour_prob_c, hltAK4PuppiJet_DeepFalvour_prob_g, hltAK4PuppiJet_DeepFalvour_prob_uds)") 

# h = df.Histo1D("btag_scores_signal")
# c = np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)])

### For each hlt jet matched to gen non-b-jet (hadronFlavour != 5), get the b tagging scores
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;
    ROOT::RVec<float> get_btag_scores_for_hlt_nonbjets (Vint matched_gen_jet_indices, Vint GenJet_hadronFlavour,
            Vfloat hltAK4PuppiJet_DeepFlavour_prob_b, Vfloat hltAK4PuppiJet_DeepFlavour_prob_bb, Vfloat hltAK4PuppiJet_DeepFlavour_prob_lepb,
            Vfloat hltAK4PuppiJet_DeepFlavour_prob_c, Vfloat hltAK4PuppiJet_DeepFlavour_prob_g, Vfloat hltAK4PuppiJet_DeepFlavour_prob_uds) {
        ROOT::RVec<float> btag_scores(matched_gen_jet_indices.size(), -1.0);
        for (size_t i_hlt = 0; i_hlt < matched_gen_jet_indices.size(); ++i_hlt) {
            if ((matched_gen_jet_indices.at(i_hlt) != -1) && std::abs(GenJet_hadronFlavour.at(matched_gen_jet_indices.at(i_hlt))) != 5) {
                float b_vs_all_score = (hltAK4PuppiJet_DeepFlavour_prob_b.at(i_hlt) + hltAK4PuppiJet_DeepFlavour_prob_bb.at(i_hlt) +
                                    hltAK4PuppiJet_DeepFlavour_prob_lepb.at(i_hlt)) /
                                   (hltAK4PuppiJet_DeepFlavour_prob_b.at(i_hlt) + hltAK4PuppiJet_DeepFlavour_prob_bb.at(i_hlt) +
                                    hltAK4PuppiJet_DeepFlavour_prob_lepb.at(i_hlt) + hltAK4PuppiJet_DeepFlavour_prob_c.at(i_hlt) +
                                    hltAK4PuppiJet_DeepFlavour_prob_g.at(i_hlt) + hltAK4PuppiJet_DeepFlavour_prob_uds.at(i_hlt));
                btag_scores[i_hlt] = static_cast<float>(b_vs_all_score);
                // std::cout << "HLT Jet index " << i_hlt << ": matched to gen jet index " << matched_gen_jet_indices.at(i_hlt)
                    // << " with hadronFlavour = " << GenJet_hadronFlavour.at(matched_gen_jet_indices.at(i_hlt))
                    // << " and b-tagging score = " << btag_scores.at(i_hlt) << std::endl;
            }
        }
        return btag_scores;
    }
""")

df = df.Define("btag_scores_background",
               "get_btag_scores_for_hlt_nonbjets(matched_gen_jet_indices, GenJet_hadronFlavour, " \
               "hltAK4PuppiJet_DeepFlavour_prob_b, hltAK4PuppiJet_DeepFalvour_prob_bb, hltAK4PuppiJet_DeepFalvour_prob_lepb, " \
               "hltAK4PuppiJet_DeepFalvour_prob_c, hltAK4PuppiJet_DeepFalvour_prob_g, hltAK4PuppiJet_DeepFalvour_prob_uds)") 

# h = df.Histo1D("btag_scores_background")
# c = np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)])

### Plot b tagging scores

df = df.Define("btag_scores_signal_clean", "ROOT::RVec<float>(btag_scores_signal) [btag_scores_signal != -1]")
h_signal = df.Histo1D("btag_scores_signal_clean")
# c_signal = np.array([h_signal.GetBinContent(i)/h_signal.Integral() for i in range(1, h_signal.GetNbinsX()+1)])
c_signal = np.array([h_signal.GetBinContent(i) for i in range(1, h_signal.GetNbinsX()+1)])
e_signal = np.array([h_signal.GetBinLowEdge(i) for i in range(1, h_signal.GetNbinsX()+2)])

df = df.Define("btag_scores_background_clean", "ROOT::RVec<float>(btag_scores_background) [btag_scores_background != -1]")
h_background = df.Histo1D("btag_scores_background_clean")
# c_background = np.array([h_background.GetBinContent(i)/h_background.Integral() for i in range(1, h_background.GetNbinsX()+1)])
c_background = np.array([h_background.GetBinContent(i) for i in range(1, h_background.GetNbinsX()+1)])
e_background = np.array([h_background.GetBinLowEdge(i) for i in range(1, h_background.GetNbinsX()+2)])

fig, ax = plt.subplots(figsize=(10, 10))
hep.histplot(c_signal, bins=e_signal, label="HLT jets matched to gen b-jets", ax=ax, histtype='step', linewidth=2)
hep.histplot(c_background, bins=e_background, label="HLT jets matched to gen non-b-jets", ax=ax, histtype='step', linewidth=2)
ax.set_xlabel("BTag Score")
ax.set_ylabel("Number of Jets")
ax.legend()
plt.grid()
plt.savefig(f"{output_dir}/b_vs_all_scores.png")
plt.savefig(f"{output_dir}/b_vs_all_scores.pdf")
