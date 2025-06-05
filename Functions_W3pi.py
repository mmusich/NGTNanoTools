import ROOT

#####################################################
# Print summary of the event, including hard scatter particles and W boson decay products

ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    int event_summary(Vint GenPart_pdgId, Vint GenPart_status, Vint GenPart_genPartIdxMother, Vfloat GenPart_mass, Vfloat GenPart_pt, Vfloat GenPart_eta) {
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
                        std::cout << "     idx " << j << ": pdgId = " << GenPart_pdgId[j] 
                        << ", GenPart_pt = " << GenPart_pt[j]
                        << ", GenPart_eta = " << GenPart_eta[j]
                        << " (mass = " << GenPart_mass[j]
                        << ", GenPart_status = " << GenPart_status[j] << ")" << std::endl;
                    }
                }
            }
        }
        std::cout << "=====================" << std::endl;
        return 1;
    }
""")

#####################################################
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

#####################################################
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

#####################################################
### For each gen pion, among the hlt PF candidates with dR < 0.4 and dPt < 10 GeV, find the closest one in dR
ROOT.gInterpreter.Declare("""
    using Vfloat = const ROOT::RVec<float>&;
    using Vint   = const ROOT::RVec<int>&;

    ROOT::RVec<int> match_gen_to_hlt (Vint gen_pion_indices, Vfloat GenPart_pt, Vfloat GenPart_eta, Vfloat GenPart_phi, Vfloat GenPart_mass,
        Vfloat hltPFCandidate_pdgId, Vfloat hltPFCandidate_pt, Vfloat hltPFCandidate_eta, Vfloat hltPFCandidate_phi, Vfloat hltPFCandidate_mass) {

        const double dR_max = 0.4; // dR < dR_max
        ROOT::RVec<int> matched_hlt_pion_indices (gen_pion_indices.size(), -1);

        std::cout << "Start event" << std::endl;
        for (size_t i_gen_pion = 0; i_gen_pion < gen_pion_indices.size(); i_gen_pion ++) {
            if (gen_pion_indices.at(i_gen_pion) == -1) continue; // Skip if no gen pion found
            auto gen_tlv = TLorentzVector();
            gen_tlv.SetPtEtaPhiM(GenPart_pt.at(gen_pion_indices.at(i_gen_pion)), GenPart_eta.at(gen_pion_indices.at(i_gen_pion)), GenPart_phi.at(gen_pion_indices.at(i_gen_pion)), GenPart_mass.at(gen_pion_indices.at(i_gen_pion)));
            double matched_dR = std::numeric_limits<double>::infinity();
            int matched_i_hlt = -1;
            for (size_t i_hlt = 0; i_hlt < hltPFCandidate_pt.size(); i_hlt ++) {
                if (std::abs(hltPFCandidate_pdgId.at(i_hlt)) != 211) continue; // Only consider pions
                // if (hltPFCandidate_pt.at(i_hlt) < 1) continue; // Reduce combinatorics
                auto hlt_tlv = TLorentzVector();
                hlt_tlv.SetPtEtaPhiM(hltPFCandidate_pt.at(i_hlt), hltPFCandidate_eta.at(i_hlt), hltPFCandidate_phi.at(i_hlt), hltPFCandidate_mass.at(i_hlt));
                double dR = gen_tlv.DeltaR(hlt_tlv);
                if (dR > dR_max) continue;
                std::cout << "Gen Jet index " << gen_pion_indices.at(i_gen_pion) 
                    << ": HLT jet index " << i_hlt << " with dR = " << dR 
                    << " (Gen pt = " << GenPart_pt.at(gen_pion_indices.at(i_gen_pion))
                    << ", HLT pt = " << hltPFCandidate_pt.at(i_hlt) << ")" << std::endl;
                double dPt = std::abs(gen_tlv.Pt() - hlt_tlv.Pt());
                if (dPt > 10) continue;
                if (dR < matched_dR) {
                    matched_dR = dR;
                    matched_i_hlt = static_cast<int>(i_hlt);
                }
            }
            matched_hlt_pion_indices[i_gen_pion] = matched_i_hlt;
            std::cout << "Gen Jet index " << gen_pion_indices.at(i_gen_pion) << ": matched to hlt jet index " << matched_i_hlt << " with dR = " << matched_dR << std::endl;
        }
        return matched_hlt_pion_indices;
    }
""")