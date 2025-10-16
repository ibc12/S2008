#ifndef PIPE3_RPCUTS_H
#define PIPE3_RPCUTS_H
#include "ActCutsManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TColor.h"
#include "THStack.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

#include <map>
#include <string>
#include <vector>

#include "../HistConfig.h"

void Pipe3_RPCuts(const std::string& beam, const std::string& target, const std::string& light)
{
    auto infile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", infile.Data()};


    // Define intervals and histograms
    double xmin {0};
    double xmax {200};
    double step {15};
    std::vector<std::pair<double, double>> ivs;
    std::map<int, ROOT::TThreadedObject<TH2D>> hs, hsEpR;
    std::map<int, ROOT::TThreadedObject<TH1D>> hsebeam;
    auto* hmodel {new TH2D {"hKin", "kin", 300, 0, 90, 300, 0, 14}};
    auto* hebeam {new TH1D {"hEBeam", "E beam", 300, 0, 90}};
    auto hEpR {HistConfig::EpRMg.GetHistogram()};
    int idx {};
    for(double x = xmin; x <= xmax; x += step)
    {
        // Compute interval
        std::pair<double, double> iv {x, x + step};
        ivs.push_back(iv);
        // Histogram
        hs.emplace(idx, *hmodel);
        hs[idx]->SetTitle(TString::Format("Kin for RP.X in [%.2f, %.2f)", iv.first, iv.second));
        hsebeam.emplace(idx, *hebeam);
        // Ep vs Range
        hsEpR.emplace(idx, *hEpR);
        hsEpR[idx]->SetTitle(TString::Format("RP.X in [%.2f, %.2f)", iv.first, iv.second));

        idx++;
    }

    // Enable slot 0
    for(auto& [_, h] : hs)
        h.GetAtSlot(0);
    for(auto& [_, h] : hsebeam)
        h.GetAtSlot(0);

    // Fill histograms
    df.Foreach(
        [&](ActRoot::MergerData& mer, double elab, double ebeam, float range)
        {
            auto rpx {mer.fRP.X()};
            for(int i = 0; i < ivs.size(); i++)
            {
                if(ivs[i].first <= rpx && rpx < ivs[i].second)
                {
                    hs[i].Get()->Fill(mer.fThetaLight, elab);
                    hsebeam[i].Get()->Fill(ebeam);
                    hsEpR[i].Get()->Fill(range, elab);
                }
            }
        },
        {"MergerData", "EVertex", "EBeam", "RangeHeavy"});

    // Get kinematics
    std::vector<ActPhysics::Kinematics> kins;
    for(auto& [i, h] : hsebeam)
    {
        h.Merge();
        auto ebeam {h.GetAtSlot(0)->GetMean()};
        kins.push_back(ActPhysics::Kinematics(TString::Format("20Na(p,p)@%.2f", ebeam).Data()));
    }

    // Plot
    auto* c0 {new TCanvas {"c0", "EBeam canvas"}};
    c0->DivideSquare(hsebeam.size());
    for(int i = 0; i < hsebeam.size(); i++)
    {
        c0->cd(i + 1);
        hsebeam[i].GetAtSlot(0)->DrawClone();
    }

    auto* c1 {new TCanvas {"c2", "Kin canvas"}};
    c1->DivideSquare(hs.size());
    for(int i = 0; i < hs.size(); i++)
    {
        c1->cd(i + 1);
        hs[i].Merge()->DrawClone("colz");
        auto* theo {kins[i].GetKinematicLine3()};
        theo->Draw("l");
    }

    auto* c2 {new TCanvas {"c2", "Ep vs range"}};
    c2->DivideSquare(hsEpR.size());
    for(int i = 0; i < hsEpR.size(); i++)
    {
        c2->cd(i + 1);
        hsEpR[i].Merge()->DrawClone("colz");
    }
}

#endif
