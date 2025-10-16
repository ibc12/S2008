#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadedObject.hxx>

#include "TCanvas.h"
#include "TColor.h"
#include "TExec.h"
#include "THStack.h"
#include "TStyle.h"
#include "TVirtualPad.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

#include <vector>

#include "../PostAnalysis/HistConfig.h"

void plotECM_intervalsThetaCM()
{
    // Read data
    auto filename {TString::Format("../PostAnalysis/Outputs/tree_ex_20Na_p_p.root")};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", filename};

    // Get histograms of Ecm on intervals of RP.x()
    std::vector<ROOT::TThreadedObject<TH1D>*> hECMs;
    double step {10}; // deg
    double Thetamin {40};
    double Thetamax {170};
    int idx {};
    for(double theta = Thetamin; theta < Thetamax; theta += step)
    {
        hECMs.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hECM%d", idx),
            TString::Format("#theta_{CM} [%.2f, %.2f);E_{CM} [MeV];Counts / 10 keV", theta, theta + step),
            HistConfig::ECM.fNbinsX, HistConfig::ECM.fXLow, HistConfig::ECM.fXUp));
        idx++;
    }
    // Initialize slot 0 to not crash
    for(auto& h : hECMs)
        h->GetAtSlot(0)->GetEntries();
    // Fill histograms
    auto hECM2d {df.Histo2D(HistConfig::ThetaCMECM, "ThetaCM", "ECM")};
    df.ForeachSlot(
        [&](unsigned int slot, double thetaCM, double ecm)
        {
            // get the hstogram we have to fill
            for(size_t i = 0; i < hECMs.size(); i++)
            {
                double theta = Thetamin + i * step;
                if(thetaCM >= theta && thetaCM < theta + step)
                {
                    hECMs[i]->GetAtSlot(slot)->Fill(ecm);
                }
            }
        },
        {"ThetaCM", "ECM"});


    // Styling options
    gStyle->SetPalette(kRainBow);

    auto* stack {new THStack};
    stack->SetTitle("Stacked E_{CM};E_{CM} [MeV];Counts / 50 keV");
    // Plot them in canvas
    auto* c0 {new TCanvas("c0", "Ecm intervals canvas")};
    c0->DivideSquare(hECMs.size());
    int p {1};
    for(auto& h : hECMs)
    {
        c0->cd(p);
        auto merged {h->Merge()};
        merged->Rebin(5);
        // merged->GetXaxis()->SetRangeUser(0.6, 4.1);
        merged->DrawClone();
        auto* clone {(TH1D*)merged->Clone()};
        clone->SetLineWidth(2);
        stack->Add(clone);
        p++;
    }

    auto* c1 {new TCanvas {"c1", "Ecm 2d canvas"}};
    c1->DivideSquare(2);
    c1->cd(1);
    hECM2d->DrawClone("colz");
    c1->cd(2);
    stack->Draw("plc pmc");
    gPad->BuildLegend();
}
