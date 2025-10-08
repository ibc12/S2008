#include "ActCluster.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <atomic>
#include <stdexcept>
#include <utility>

void Pipe0_Beam(const std::string& beam)
{
    std::string dataconf {"./../configs/data.conf"};

    ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman {dataconf, ActRoot::ModeType::EReadSilMod};
    auto chain {datman.GetJoinedData()};
    auto chain2 {datman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3 {datman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    ROOT::RDataFrame df {*chain};

    // Get GATCONF values
    auto defGat {df.Define("GATCONF", [](ActRoot::ModularData& mod)
                           { return static_cast<int>(mod.fLeaves["GATCONF"]); }, {"ModularData"})};

    // Get plots for DE/E for beam, first 10 pads, until x = 10.
    auto defBeam {defGat.Filter("fBeamIdx != -1")
                      .Define("Pair",
                              [](ActRoot::TPCData& d, ActRoot::MergerData& m)
                              {
                                  int beamIdx {m.fBeamIdx};
                                  auto beamCluster {d.fClusters[beamIdx]};
                                  auto voxels {beamCluster.GetRefToVoxels()};
                                  double dE {0}; // DeltaE in first 10 pads
                                  double E {};   // Total energy
                                  for(auto& v : voxels)
                                  {
                                      auto q {v.GetCharge()};
                                      if(v.GetPosition().X() < 10)
                                          dE += q;
                                      E += q;
                                  }
                                  return std::make_pair(dE, E);
                              },
                              {"TPCData", "MergerData"})
                      .Define("dE", "Pair.first")
                      .Define("E", "Pair.second")};

    // Book histograms
    auto hGATCONF {defGat.Histo1D("GATCONF")};
    auto hdEE {
        defBeam.Histo2D({"hdEE", "Beam ID;Q_{total} [au];Q_{10 pads} [au]", 300, 0, 1e5, 300, 0, 1e5}, "E", "dE")};

    // And cound CFA triggers
    std::atomic<unsigned long int> cfa {};
    defGat.Foreach(
        [&](const int& gatconf)
        {
            if(gatconf == 64)
                cfa++;
        },
        {"GATCONF"});

    // Draw
    auto* c0 {new TCanvas {"c00", "Pipe 0 canvas 0"}};
    c0->DivideSquare(4);
    c0->cd(1);
    hdEE->DrawClone("colz");
    c0->cd(2);
    hGATCONF->DrawClone();

    // Print report
    std::cout << "===== GATCONF report =====" << '\n';
    std::cout << "-> CFA/div = " << cfa << '\n';
    std::cout << "==========================" << '\n';
}
