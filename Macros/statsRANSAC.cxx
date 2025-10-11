#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include <atomic>
#include <fstream>

#include "../PostAnalysis/HistConfig.h"
void statsRANSAC()
{
    // Get data
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    ROOT::RDataFrame df {*chain};

    // Read file
    std::map<int, std::set<int>> file;
    std::ifstream streamer {"./Outputs/gatconf_f0_true.dat"};
    int run {};
    int entry {};
    int total {};
    while(streamer >> run >> entry)
    {
        file[run].insert(entry);
        total++;
    }

    // Get recovered
    auto rec {df.Filter(
        [&](ActRoot::MergerData& d)
        {
            if(file.count(d.fRun))
            {
                auto& set {file[d.fRun]};
                if(set.find(d.fEntry) != set.end())
                {
                    if(d.fLightIdx != -1)
                        return true;
                }
            }
            return false;
        },
        {"MergerData"})};

    // Plot recovered kinematics
    auto hKin {rec.Define("y", [](ActRoot::MergerData& d) { return d.fLight.fEs.front(); }, {"MergerData"})
                   .Histo2D(HistConfig::Kin, "fThetaLight", "y")};

    // Plot
    auto* c1 {new TCanvas {"c1", "RANSAC recovered events"}};
    hKin->DrawClone("colz");

    auto recCount {*rec.Count()};
    std::cout << "% rec events: " << (double)recCount / total * 100 << '\n';
}
