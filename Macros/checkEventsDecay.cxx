#include "ActCluster.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TMath.h"

void checkEventsDecay()
{
    // Read the data using the data.conf file
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain {dataman.GetChain()}; // Get all Merge files for Runs in a single TChain
    // Add friends if necessary
    auto friend1 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(friend1.get());
    auto friend2 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(friend2.get());

    // Build the RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto gated {df.Filter(
        [](ActRoot::MergerData& m)
        {
            if(m.fLightIdx != -1)
                return true;
            return false;
        },
        {"MergerData"})};

    // Get distance from
    auto def {gated.Define("DeltaZ_RP_ClusterLight",
                           [](ActRoot::MergerData& m, ActRoot::TPCData& tpc)
                           {
                               // Assign the highest angle cluster to the light
                               int lightIdx {m.fLightIdx};
                               auto& clusters {tpc.fClusters};
                               auto& clusterLight {clusters[lightIdx]};

                               clusterLight.GetRefToLine().AlignUsingPoint(m.fRP);
                               // Now check if there is a shift in z coordinate respect the rp
                               clusterLight.SortAlongDir();
                               auto rp {m.fRP};

                               double zFirstClusterLight {clusterLight.GetRefToVoxels().front().GetPosition().Z()};
                               double distance {(zFirstClusterLight * 2.208 - rp.Z())};

                               // std::cout << "Z first voxel light cluster: " << distance << '\n';

                               return distance;
                           },
                           {"MergerData", "TPCData"})};

    // Plot the z distance
    auto hDZ {def.Histo1D(
        {"hDZ", "Delta Z RP - first voxel light cluster;distance in Z [mm];Counts", 100,
         -200, 200},
        "DeltaZ_RP_ClusterLight")};
    auto c {new TCanvas("c", "Delta Z RP - first voxel light cluster")};
    hDZ->GetEntries();
    hDZ->DrawClone();
}