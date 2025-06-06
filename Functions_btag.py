import ROOT

#####################################################
### For each hlt jet, find index of the matched gen jet 
# from https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L40

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