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


void plotECM_intervalsThetaCM()
{
    // Read data
    auto filename {TString::Format("../PostAnalysis/Outputs/tree_ex_20Mg_p_p.root")};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", filename};
    auto df_filtered = df.Filter([](ActRoot::MergerData& m) { return m.fLightIdx != -1; }, {"MergerData"});

    // Get histograms of Ecm on intervals of RP.x()
    std::vector<ROOT::TThreadedObject<TH1D>*> hECMs;
    double step {30}; // mm
    double Thetamin {0};
    double Thetamax {180};

    for(double theta = Thetamin; theta < Thetamax; theta += step)
    {
        hECMs.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hECM_%.0f-%.0f", theta, theta + step),
            TString::Format("Ecm (Theta in [%.0f, %.0f]);Ecm [MeV];Counts", theta, theta + step), 200, 0, 15));
    }
    // Initialize slot 0 to not crash
    for(auto& h : hECMs)
        h->GetAtSlot(0)->GetEntries();
    df_filtered.ForeachSlot(
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

    // Plot them all in a canvas
    auto* c {new TCanvas("c", "Ecm intervals canvas")};
    c->DivideSquare(hECMs.size());
    int p {1};
    for(auto& h : hECMs)
    {
        c->cd(p);
        h->Merge()->DrawClone();
        p++;
    }
}