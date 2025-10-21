#include "ActSRIM.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "THStack.h"
#include "TStyle.h"
#include "TVirtualPad.h"

#include <iostream>
void compBeamELoss()
{
    // 20Mg
    auto* file {new TFile {"./Outputs/beam_eloss_20Mg.root"}};
    auto* hmg {file->Get<TH1D>("hELoss20Mg")};

    auto* otherFile {new TFile {"./Outputs/beam_eloss_20Na.root"}};
    auto* hna {otherFile->Get<TH1D>("hELoss20Na")};

    // Determine 20Na scaling factor
    ActPhysics::SRIM srim;
    srim.ReadTable("950", "../Calibrations/SRIM/20Na_950mbar_95-5.txt");
    srim.ReadTable("800", "../Calibrations/SRIM/20Na_800mbar_95-5.txt");
    double length {8 * 2.};
    double Eini {84.8}; // 4.24 MeV/u * 20
    auto eloss950 {Eini - srim.Slow("950", Eini, length)};
    auto eloss800 {Eini - srim.Slow("800", Eini, length)};
    auto transPressure {eloss800 / eloss950};
    std::cout << "20Na DeltaE 950: " << eloss950 << '\n';
    std::cout << "20Na DeltaE 800: " << eloss800 << '\n';
    std::cout << "20Na scaling   : " << transPressure << '\n';

    // Rough calculation of gain trans
    auto gainNa {5420.}; // V
    auto gainMg {4600.}; // V
    auto transGain {gainMg / gainNa};

    // Transform histogram with this scale
    auto* hTransNa {(TH1D*)hna->Clone()};
    hTransNa->Reset();
    for(int b = 1; b <= hTransNa->GetNbinsX(); b++)
    {
        auto x {hTransNa->GetBinCenter(b)};
        // Transform 20Na pressure
        auto x950 {x / transPressure};
        // Transform 20Na gain
        auto xGainNa {x950 / transGain};
        // Eval y
        auto y {hna->Interpolate(xGainNa)};
        hTransNa->AddBinContent(b, y);
    }

    // Scale the 20na
    hTransNa->Scale(6.1e4 / hTransNa->Integral());

    // Stack
    auto* stack {new THStack};
    stack->SetTitle("^{20}Na vs ^{20}Mg;E_{beam} [MeV];Couns");
    hmg->SetLineColor(1);
    stack->Add(hmg);
    hTransNa->SetLineColor(8);
    stack->Add(hTransNa, "hist");
    // Fit
    hmg->Fit("gaus", "0Q", "", 10000, 15000);
    hmg->GetFunction("gaus")->ResetBit(TF1::kNotDraw);

    // Draw
    gStyle->SetOptFit(true);
    auto* c0 {new TCanvas {"c0", "Beam ELoss canvas"}};
    hmg->Draw();
    hTransNa->Draw("same");
    // stack->Draw("nostack");
    gPad->BuildLegend();
}
