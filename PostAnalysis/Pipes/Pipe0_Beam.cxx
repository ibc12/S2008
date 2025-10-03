#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActTypes.h"
#include "ActTPCData.h"
#include "ActMergerData.h"
#include "ActCluster.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <atomic>
#include <stdexcept>

void Pipe0_Beam(const std::string &beam)
{
    std::string dataconf{"./../configs/data.conf"};

    ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman{dataconf, ActRoot::ModeType::EReadSilMod};
    auto chain{datman.GetJoinedData()};
    auto chain2{datman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3{datman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    ROOT::RDataFrame df{*chain};

    // Get GATCONF values
    auto def{df.Define("GATCONF", [](ActRoot::ModularData &mod)
                       { return static_cast<int>(mod.fLeaves["GATCONF"]); },
                       {"ModularData"})};

    // Get plots for DE/E for beam, first 10 pads, until x = 10.
    auto df1{def.Define("DE",
                        [](ActRoot::TPCData &d, ActRoot::MergerData &m)
                        {
                            int beamIdx{m.fBeamIdx};
                            auto beamCluster{d.fClusters[beamIdx]};
                            auto voxels{beamCluster.GetRefToVoxels()};
                            double DE{0};
                            for (auto &v : voxels)
                            {
                                if (v.GetPosition().X() < 10)
                                    DE += v.GetCharge();
                            }
                            return DE;
                        },
                        {"TPCData", "MergerData"})
                 .Define("E",
                         [](ActRoot::TPCData &d, ActRoot::MergerData &m)
                         {
        int beamIdx {m.fBeamIdx};
        auto beamCluster {d.fClusters[beamIdx]};
        auto voxels {beamCluster.GetRefToVoxels()};
        double E {0};
        for (auto& v : voxels)
            E += v.GetCharge();
        return E; },
                         {"TPCData", "MergerData"})};

    // Book histograms
    auto hGATCONF{def.Histo1D("GATCONF")};
    auto hDE_E{df1.Histo2D({"hDE_E", "DE vs E;E;DE", 100, 0, 1e5, 100, 0, 1e5}, "E", "DE")};

    // And cound CFA triggers
    std::atomic<unsigned long int> cfa{};
    def.Foreach(
        [&](const int &gatconf)
        {
            if (gatconf == 64)
                cfa++;
        },
        {"GATCONF"});

    // Draw
    auto *c0{new TCanvas{"c00", "Pipe 0 canvas 0"}};
    hGATCONF->DrawClone();

    // Print report
    std::cout << "===== GATCONF report =====" << '\n';
    std::cout << "-> CFA/div = " << cfa << '\n';
    std::cout << "==========================" << '\n';
}
