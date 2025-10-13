#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadedObject.hxx>

#include "TCanvas.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

#include <map>
#include <vector>

#include "../PostAnalysis/HistConfig.h"

void plotECM_intervalsThetaCM()
{
    // Read data
    auto filename {TString::Format("../PostAnalysis/Outputs/tree_ex_20Mg_p_p.root")};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", filename};

    // Get histograms of Ecm on intervals of RP.x()
    std::vector<ROOT::TThreadedObject<TH1D>*> hECMs;
    double step {25}; // deg
    double Thetamin {0};
    double Thetamax {180};
    int idx {};
    for(double theta = Thetamin; theta < Thetamax; theta += step)
    {
        hECMs.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hECM%d", idx),
            TString::Format("E_{CM} for #theta_{CM} [%.2f, %.2f);E_{CM} [MeV];Counts", theta, theta + step),
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

    // Plot them in canvas
    auto* c0 {new TCanvas("c0", "Ecm intervals canvas")};
    c0->DivideSquare(hECMs.size());
    int p {1};
    for(auto& h : hECMs)
    {
        c0->cd(p);
        h->Merge()->DrawClone();
        p++;
    }

    auto* c1 {new TCanvas {"c1", "Ecm 2d canvas"}};
    hECM2d->DrawClone("colz");
}
