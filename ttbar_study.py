
## -------------------------------------------------------------------------- ##
##    Author:    A.M.M Elsayed                                                ##
##    Email:     ahmedphysica@outlook.com                                     ##
##    Institute: University of Science and Technology of China                ##
## -------------------------------------------------------------------------- ##
##                                                                            ##
##    Updates On:                                                             ##
##    https://github.com/ammelsayed/LHEAnalysis                               ##
## -------------------------------------------------------------------------- ##

from LHEAnalysis import LHEAnalysis, UseBranch

## ==================================
## Generate the truth.root file
## ==================================

event_num = "1M"

# Run only on first time
lhe = LHEAnalysis([1,-1,2])
lhe.LoadLHE(f"./data/ttbar_{event_num}.lhe.gz")
lhe.SaveAsROOT(f"./data/ttbar_truth_{event_num}.root")

## ==================================
## Analyze the truth.root file
## ==================================

import ROOT

file = ROOT.TFile.Open(f"./data/ttbar_truth_{event_num}.root")
tree = file.Get("LHEAnalysis;1")

# event weight
lum = 400 #/fb
weight = 3.743 * lum / (1e4 if event_num == "10k" else 1e6)

# Branches we are going to use. 
# The UseBranch function is helpful for writing a clean code
WBosons = UseBranch("W", tree)
BottomQuarks = UseBranch("Bottom", tree)
TopQuarks = UseBranch("Top", tree)

# booking histograms
hist_Wb_mass = ROOT.TH1D("hist_Wb_mass", "M(W^{#pm}b); M(W^{#pm}b); Events", 50, 100, 250)
hist_Wb_mass.SetLineColor(ROOT.kRed)
hist_top_mass = ROOT.TH1D("hist_top_mass", "m(t); m(t); Events", 50, 100, 250)
hist_top_mass.SetLineColor(ROOT.kBlue)
hist_Wb_PT = ROOT.TH1D("hist_Wb_PT", "; Total p_{T}; Events", 50, 400, 1000)
hist_Wb_PT.SetLineColor(ROOT.kRed)
hist_PT = ROOT.TH2D("hist_PT", "; p_{T}(top); p_{T}(W) + p_{T}(b); Events", 50, 400, 1000,  50, 400, 1000)
hist_DeltaR = ROOT.TH1D("hist_DeltaR", "; #DeltaR(W, b); Events", 50, 0, 6)
hist_DeltaR.SetLineColor(ROOT.kRed)
hist_DeltaR_Pt = ROOT.TH2D("hist_DeltaR_Pt", "; p_{T}(top); #DeltaR(W, b); Events", 50, 400, 1000,  50, 0, 6)
hist_DeltaPhi = ROOT.TH1D("hist_DeltaPhi", "; |#Delta#phi(W, b)|; Events", 50, 0, 3.2)
hist_DeltaPhi.SetLineColor(ROOT.kRed)
hist_DeltaPhi_Pt = ROOT.TH2D("hist_DeltaPhi_Pt", "; p_{T}(top); |#Delta#phi(W, b)|; Events", 50, 400, 1000,  50, 0, 3.2)
hist_Asy = ROOT.TH1D("hist_Asy", "; Asy[p_{T}(W) - p_{T}(b)]; Events", 50, -1, 1)
hist_Asy.SetLineColor(ROOT.kRed)
hist_Asy_Pt = ROOT.TH2D("hist_Asy_Pt", "; p_{T}(top); Asy[p_{T}(W) - p_{T}(b)]; Events", 50, 400, 1000,  50, -1, 1)

def From(pdgID, p) -> None:
    """A function to help you identify the mother particle"""
    return (int(p.pdgId_Mother1) == pdgID and int(p.pdgId_Mother2) == pdgID)

def Asymmetry(Var, obj1, obj2):
    Sum = getattr(obj1, Var) + getattr(obj2, Var)
    Dif = getattr(obj1, Var) - getattr(obj2, Var)
    return Dif/Sum if Sum > 0 else 0

for event in range(tree.GetEntries()):

    tree.GetEntry(event)
    
    nW = WBosons.GetEntriesFast()
    nBottom = BottomQuarks.GetEntriesFast()
    nTop= TopQuarks.GetEntriesFast()

    if nTop >= 0 and nW > 0 and nBottom > 0:

        Cand_Top = TopQuarks.At(0)
        hist_top_mass.Fill(Cand_Top.P4().M(), weight)

        TopPT = Cand_Top.PT
        
        found = False
        for i in range(nW):
            for j in range(nBottom):

                Cand_W = WBosons.At(i)
                Cand_Bottom = BottomQuarks.At(j)

                # Make sure they both come from same top quark (6)
                if not (From(6, Cand_W) and From(6, Cand_Bottom)):
                    continue

                # Fill invariant mass histogram
                InvMass = (Cand_W.P4() + Cand_Bottom.P4()).M()
                hist_Wb_mass.Fill(InvMass, weight)

                # Fill Total PT histogram
                TotalPT = Cand_W.PT + Cand_Bottom.PT
                hist_Wb_PT.Fill(TotalPT, weight)
                hist_PT.Fill(TopPT, TotalPT, weight)

                # Fill DeltaR histogram
                DeltaR = Cand_W.P4().DeltaR(Cand_Bottom.P4(), useRapidity = False)
                hist_DeltaR.Fill(DeltaR, weight)
                hist_DeltaR_Pt.Fill(TopPT, DeltaR, weight)

                # Fill DeltaPhi histogram
                DeltaPhi = Cand_W.P4().DeltaPhi(Cand_Bottom.P4())
                hist_DeltaPhi.Fill(DeltaPhi, weight)
                hist_DeltaPhi_Pt.Fill(TopPT, DeltaPhi, weight)

                # Fill transvere momentum assymetry data:
                pTAsy = Asymmetry("PT", Cand_W, Cand_Bottom)
                hist_Asy.Fill(pTAsy, weight)
                hist_Asy_Pt.Fill(TopPT, pTAsy, weight)

                found = True; break
            if found: break


# Saving plots
ROOT.gROOT.SetStyle("ATLAS")
ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
canvas = ROOT.TCanvas()
canvas.SetCanvasSize(800*2,600*5)
canvas.Divide(2,5)
canvas.cd(1).SetLogy(True)
hist_top_mass.Draw("hist")
canvas.cd(2).SetLogy(True)
hist_Wb_mass.Draw("hist")
canvas.cd(3).SetLogy(True)
hist_Wb_PT.Draw("hist")
canvas.cd(4).SetLogy(False)
hist_PT.Draw("colz")
canvas.cd(5).SetLogy(True)
hist_DeltaR.Draw("hist")
canvas.cd(6).SetLogy(False)
hist_DeltaR_Pt.Draw("colz")
canvas.cd(7).SetLogy(True)
hist_DeltaPhi.Draw("hist")
canvas.cd(8).SetLogy(False)
hist_DeltaPhi_Pt.Draw("colz")
canvas.cd(9).SetLogy(True)
hist_Asy.Draw("hist")
canvas.cd(10).SetLogy(False)
hist_Asy_Pt.Draw("colz")
canvas.SaveAs(f"./data/Characterstics_{event_num}.pdf")
canvas.SaveAs(f"./data/Characterstics_{event_num}.png")
file.Close()
canvas.Close()