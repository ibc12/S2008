#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include <fstream>
#include <memory>

void beta27P()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain {dataman.GetChain()};
    auto chainTpc {dataman.GetChain(ActRoot::ModeType::EReadTPC)};
    auto chainMerger {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chainTpc.get());
    chain->AddFriend(chainMerger.get());

    ROOT::RDataFrame d {*chain};
    auto df {d.Define("GATCONF", [](ActRoot::ModularData& mod) { return mod.Get("GATCONF"); }, {"ModularData"})};
    // Gate on silicon events
    // And validate with thresholds
    auto specs {std::make_shared<ActPhysics::SilSpecs>()};
    specs->ReadFile("../configs/silspecs.conf");
    auto dfGat {df.Filter(
        [&](float gatconf, ActRoot::SilData& sil)
        {
            // Validate sil data
            sil.ApplyFinerThresholds(specs);
            bool isVal {false};
            for(const auto& [key, es] : sil.fSiE)
                if(sil.GetMult(key) == 1)
                    isVal = true;
            return (gatconf == 512) && isVal;
        },
        {"GATCONF", "SilData"})};

    // Gate on events with a track = proton!
    auto dfProton {dfGat.Filter([](ActRoot::TPCData& d) { return d.fClusters.size() == 1; }, {"TPCData"})};

    // Write events
    std::ofstream streamer {"./Outputs/maybe_betas.dat"};
    dfProton.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    streamer.close();

    std::cout << "% proton events: " << (double)(*dfProton.Count()) / (*dfGat.Count()) * 100 << '\n';
}
