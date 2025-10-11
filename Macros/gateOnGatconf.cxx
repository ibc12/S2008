#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>

void gateOnGatconf()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain {dataman.GetChain()};
    auto chainMerger {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chainMerger.get());

    ROOT::RDataFrame df {*chain};

    auto dfFilter {df.Filter(
        [](ActRoot::SilData& sil, ActRoot::ModularData& mod)
        {
            if(mod.Get("GATCONF") == 4) // front silicons
                return true;
            return false;
        },
        {"SilData", "ModularData"})};

    // Stream entry number
    std::ofstream streamer {"./Outputs/gatconf_f0.dat"};
    dfFilter.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    streamer.close();
}
