#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActTPCData.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TMath.h"

void plotELossBeam_20Ne()
{
    ROOT::EnableImplicitMT();
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chain2.get());
    ROOT::RDataFrame df {*chain};

    // Create columns for E_Loss and E_Beam
    auto df_gated1 = df.Filter([](ActRoot::ModularData& m) { return m.Get("GATCONF") == 64; }, {"ModularData"});

    auto df_gated = df_gated1.Filter([](ActRoot::TPCData& d) { return d.fClusters.size() != 0; }, {"TPCData"});

    auto df_final = df_gated
                        .Define("E_Loss",
                                [](ActRoot::TPCData& d)
                                {
                                    auto voxelsBeam {d.fClusters[0].GetRefToVoxels()};
                                    double DE {0};
                                    for(auto v : voxelsBeam)
                                    {
                                        if(v.GetPosition().X() < 10)
                                        {
                                            DE += v.GetCharge();
                                        }
                                    }
                                    return DE;
                                },
                                {"TPCData"})
                        .Define("E_Beam",
                                [](ActRoot::TPCData& d)
                                {
                                    auto voxelsBeam {d.fClusters[0].GetRefToVoxels()};
                                    double E {0};
                                    for(auto v : voxelsBeam)
                                    {
                                        E += v.GetCharge();
                                    }
                                    return E;
                                },
                                {"TPCData"});

    auto hDE_E {df_final.Histo1D({"hDE_20Ne", "DE;DE [u.a.]", 100, 0, 1e5}, "E_Loss")};
    // Save histos in files
    auto f_out {std::make_unique<TFile>("./Outputs/DE_20Ne_10pads.root", "recreate")};
    f_out->WriteObject(hDE_E.GetPtr(), "hDE_20Ne");
    auto g {df_final.Graph("E_Beam", "E_Loss")};

    // Plot
    auto c {new TCanvas("c")};
    hDE_E->DrawClone();

    // Draw
    auto c1 = new TCanvas("c1", "E_Loss vs E_Beam", 800, 600);
    g->SetMarkerColor(kRed);
    auto* clone = g->DrawClone("AP");

    // Legend
    auto leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(clone, "^{20}Ne", "p");
    leg->Draw();
}
