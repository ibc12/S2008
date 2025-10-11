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


void plotECM_intervalsRP()
{
    // Read data
    auto filename {TString::Format("../PostAnalysis/Outputs/tree_ex_20Mg_p_p.root")};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", filename};
    auto df_filtered = df.Filter([](ActRoot::MergerData& m) { return m.fLightIdx != -1; }, {"MergerData"});

    // Get histograms of Ecm on intervals of RP.x()
    std::vector<ROOT::TThreadedObject<TH1D>*> hECMs;
    double step {20}; // mm
    double RPmin {20};
    double RPmax {256};

    for(double rp = RPmin; rp < RPmax; rp += step)
    {
        hECMs.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hECM_%.0f-%.0f", rp, rp + step),
            TString::Format("Ecm (RP.x in [%.0f, %.0f]);Ecm [MeV];Counts", rp, rp + step), 200, 0, 15));
    }
    // Initialize slot 0 to not crash
    for(auto& h : hECMs)
        h->GetAtSlot(0)->GetEntries();
    df_filtered.ForeachSlot(
        [&](unsigned int slot, ActRoot::MergerData& m, double ecm)
        {
            auto rp_x {m.fRP.X()};
            // get the hstogram we have to fill
            for(size_t i = 0; i < hECMs.size(); i++)
            {
                double rp = RPmin + i * step;
                if(rp_x >= rp && rp_x < rp + step)
                {
                    hECMs[i]->GetAtSlot(slot)->Fill(ecm);
                }
            }
        },
        {"MergerData", "ECM"});

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