#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include <memory>

void rangeFromCFA()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain {dataman.GetChain()};
    auto chainTPC {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(chainTPC.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // Gate on CFAdiv = 64
    // and define beam track length
    auto gated {df.Filter([](ActRoot::ModularData& mod) { return mod.Get("GATCONF") == 64; }, {"ModularData"})
                    .Filter([](ActRoot::TPCData& tpc) { return tpc.fClusters.size() == 1; }, {"TPCData"})
                    .Define("LastPad", [](ActRoot::TPCData& tpc)
                            { return 2 * tpc.fClusters.front().GetXRange().second; }, {"TPCData"})};

    auto hLast {gated.Histo1D({"hLast", "Range from CFA;Range [mm]", 100, 0, 250}, "LastPad")};

    // Draw
    auto* c0 {new TCanvas {"c0", "Range from CFA"}};
    hLast->DrawClone();

    // auto fout {std::make_unique<TFile>("./Outputs/range_before_after.root", "update")};
    // fout->cd();
    // hLast->Write("hAfter");
    // fout->Close();
}
