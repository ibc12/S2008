#include "ActMergerData.h"

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
    auto nodeSil {df.Filter([](ActRoot::MergerData& mer) { return !mer.fLight.IsL1(); }, {"MergerData"})};
    auto nodeL1 {df.Filter([](ActRoot::MergerData& mer) { return mer.fLight.IsL1(); }, {"MergerData"})};

    // Get histograms of Ecm on intervals of RP.x()
    std::vector<ROOT::TThreadedObject<TH1D>*> hsSil, hsL1;
    double step {10}; // deg
    double Thetamin {40};
    double Thetamax {180};
    int idx {};
    for(double theta = Thetamin; theta < Thetamax; theta += step)
    {
        hsSil.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hECM%d", idx),
            TString::Format("Silicon #theta_{CM} [%.2f, %.2f);E_{CM} [MeV];Counts / 10 keV", theta, theta + step),
            HistConfig::ECM.fNbinsX, HistConfig::ECM.fXLow, HistConfig::ECM.fXUp));
        hsL1.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hECM%d", idx),
            TString::Format("L1 #theta_{CM} [%.2f, %.2f);E_{CM} [MeV];Counts / 10 keV", theta, theta + step),
            HistConfig::ECM.fNbinsX, HistConfig::ECM.fXLow, HistConfig::ECM.fXUp));
        idx++;
    }
    // Initialize slot 0 to not crash
    for(auto& h : hsSil)
        h->GetAtSlot(0);
    for(auto& h : hsL1)
        h->GetAtSlot(0);

    // Fill histograms
    auto hECM2d {df.Histo2D(HistConfig::ThetaCMECM, "ThetaCM", "ECM")};
    // Silicons
    nodeSil.ForeachSlot(
        [&](unsigned int slot, double thetaCM, double ecm)
        {
            // get the hstogram we have to fill
            for(size_t i = 0; i < hsSil.size(); i++)
            {
                double theta = Thetamin + i * step;
                if(thetaCM >= theta && thetaCM < theta + step)
                {
                    hsSil[i]->GetAtSlot(slot)->Fill(ecm);
                }
            }
        },
        {"ThetaCM", "Rec_ECM"});
    // L1
    nodeL1.ForeachSlot(
        [&](unsigned int slot, double thetaCM, double ecm)
        {
            // get the hstogram we have to fill
            for(size_t i = 0; i < hsSil.size(); i++)
            {
                double theta = Thetamin + i * step;
                if(thetaCM >= theta && thetaCM < theta + step)
                {
                    hsL1[i]->GetAtSlot(slot)->Fill(ecm);
                }
            }
        },
        {"ThetaCM", "Rec_ECM"});


    // Styling options
    gStyle->SetPalette(kRainBow);

    // Silicon stack
    auto* stack {new THStack};
    stack->SetTitle("Stacked E_{CM};E_{CM} [MeV];Counts / 50 keV");

    // Plot them in canvas
    auto* c0 {new TCanvas("c0", "Silicon ECM")};
    c0->DivideSquare(hsSil.size());
    int p {1};
    for(auto& h : hsSil)
    {
        c0->cd(p);
        auto merged {h->Merge()};
        // merged->GetXaxis()->SetRangeUser(0.6, 4.1);
        merged->DrawClone();
        auto* clone {(TH1D*)merged->Clone()};
        clone->SetLineWidth(2);
        stack->Add(clone);
        p++;
    }

    // L1 canvas
    auto* c1 {new TCanvas("c1", "L1 ECM")};
    c1->DivideSquare(hsL1.size());
    p = 1;
    for(auto& h : hsL1)
    {
        c1->cd(p);
        auto merged {h->Merge()};
        merged->DrawClone();
        p++;
    }

    auto* c2 {new TCanvas {"c2", "Ecm 2d canvas"}};
    c2->DivideSquare(2);
    c2->cd(1);
    hECM2d->DrawClone("colz");
    c2->cd(2);
    stack->Draw("plc pmc");
    gPad->BuildLegend();
}
