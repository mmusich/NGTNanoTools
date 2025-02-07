import ROOT, sys

fileName = sys.argv[1]
file = ROOT.TFile.Open(fileName)
if not file or file.IsZombie():
    print(f" ### ERROR: Error opening file: {fileName}!")
    exit()

print(f' | {"Sample":<20} | {"Eta":<9} | At   | {"On":<50} | {"Efficiency":<15} | ')
print(f' | {"-"*20} | {"-"*9} | {"-"*4} | {"-"*50} | {"-"*15} | ')

EtaLimits = fileName.split("Eta")[1].split(".root")[0].split("_")
LowEta = EtaLimits[0].replace("p", ".")
HighEta = EtaLimits[1].replace("p", ".")
for dir in ["hltPhase2Pixel", "hltInitialStepTrackSelectionHighPurity", "hltHighPtTripletStepTrackSelectionHighPurity", "hltGeneral"]:
    mainDir = f"DQMData/Run 1/HLT/Run summary/Tracking/ValidationWRTtp"
    hist = file.Get(f"{mainDir}/{dir}_hltAssociatorByHits/globalEfficiencies")
    if not hist:
        print(f" ### WARNING: Histogram '{hist}' not found!")
        continue
    eff = hist.GetBinContent(1)
    eff_err = hist.GetBinError(1)
    print(f' | {"SingleMuPt15":<20} | [{LowEta:<3},{HighEta:<3}] | HLT  | {dir:<50} | {eff:.2f} +/- {eff_err:.4f} | ')
    