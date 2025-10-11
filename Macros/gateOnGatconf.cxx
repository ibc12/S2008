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
    auto dfE {dfFilter.Filter(
        [](ActRoot::SilData& sil, ActRoot::MergerData& mer)
        {
            if(sil.fSiE.count("f0"))
            {
                if(sil.fSiE["f0"].size() == 1)
                    if(sil.fSiE["f0"].front() >= 12)
                        // if(mer.fSP.X() > 0 && mer.fSP.X() < 45)
                        return true;
            }
            return false;
        },
        {"SilData", "MergerData"})};

    // Stream entry number
    std::ofstream streamer {"./Outputs/gatconf_f0.dat"};
    dfFilter.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    streamer.close();
    std::ofstream streamer1 {"./Outputs/gatconf_f0_true.dat"};
    dfE.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer1); }, {"MergerData"});
    streamer1.close();
}
