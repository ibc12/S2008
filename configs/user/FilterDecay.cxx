#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActRANSAC.h"
#include "ActTPCData.h"

#include "FilterDecay.h"

#include <memory>

void ActAlgorithm::FilterDecay::FilterDecay::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
    if(block->CheckTokenExists("MinLength"))
        fMinLength = block->GetDouble("MinLength");
}

void ActAlgorithm::FilterDecay::Run()
{
    if(!fIsEnabled)
        return;

    if(!fTPCData) {
        if(fIsVerbose) std::cerr << "FilterDecay: fTPCData is null\n";
        return;
    }

    // copy/ref clusters
    auto &clusters = fTPCData->fClusters;
    if(clusters.empty()) {
        if(fIsVerbose) std::cout << "FilterDecay: no clusters\n";
        return;
    }

    // Assign the highest angle cluster to the light
    ActRoot::Cluster clusterLight = clusters.front();
    double angleLight = -1.0;
    for(const auto &c : clusters)
    {
        // defensively get direction
        auto dir = c.GetLine().GetDirection();
        double angleCluster = std::abs(dir.Dot(ROOT::Math::XYZVectorF(1, 0, 0)));
        if(angleCluster > angleLight)
        {
            angleLight = angleCluster;
            clusterLight = c;
        }
    }

    clusterLight.SortAlongDir();

    // guard RPs
    if(fTPCData->fRPs.empty()) {
        if(fIsVerbose) std::cout << "FilterDecay: no RPs available\n";
        return;
    }

    // safe copy of rp
    auto rp = fTPCData->fRPs[0];

    // guard voxels
    const auto &voxels = clusterLight.GetRefToVoxels();
    if(voxels.empty()) {
        if(fIsVerbose) std::cout << "FilterDecay: selected cluster has no voxels\n";
        return;
    }

    double zFirstClusterLight = voxels.front().GetPosition().Z();

    if(std::abs(zFirstClusterLight - rp.Z()) > fMinLength)
    {
        if(fIsVerbose)
        {
            std::cout << BOLDGREEN << "-- Decay filter --\n";
            std::cout << "First Z cluster light: " << zFirstClusterLight << '\n';
            std::cout << "Z RP: " << rp.Z() << '\n';
            std::cout << "Delta Z: " << std::abs(zFirstClusterLight - rp.Z()) << RESET << '\n';
        }
        // Clear all clusters and push only the light one
        clusters.clear();
    }
    else
    {
        if(fIsVerbose)
        {
            std::cout << BOLDRED << "-- Decay filter --\n";
            std::cout << "First Z cluster light: " << zFirstClusterLight << '\n';
            std::cout << "Z RP: " << rp.Z() << '\n';
            std::cout << "Delta Z: " << std::abs(zFirstClusterLight - rp.Z()) << RESET << '\n';
            std::cout << "Event rejected by decay filter\n";
        }
    }
}


void ActAlgorithm::FilterDecay::Print() const
{
    std::cout << BOLDCYAN << "····· " << GetActionID() << " ·····" << '\n';
    if(!fIsEnabled)
    {
        std::cout << "······························" << RESET << '\n';
        return;
    }
    std::cout << "  MinLength      : " << fMinLength << '\n';
}

// Create symbol to load class from .so
extern "C" ActAlgorithm::FilterDecay* CreateUserAction()
{
    return new ActAlgorithm::FilterDecay;
}
