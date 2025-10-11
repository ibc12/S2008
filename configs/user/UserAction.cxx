#include "UserAction.h"

#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActRANSAC.h"
#include "ActTPCData.h"

#include <memory>

void ActAlgorithm::UserAction::UserAction::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
    // if(block->CheckTokenExists("MaxAngle"))
    //     fMaxAngle = block->GetDouble("MaxAngle");
    // if(block->CheckTokenExists("MinLength"))
    //     fMinLength = block->GetDouble("MinLength");
    // if(block->CheckTokenExists("CylinderR"))
    //     fCylinderR = block->GetDouble("CylinderR");
}

void ActAlgorithm::UserAction::Run()
{
    if(!fIsEnabled)
        return;

    // Enable only if only one cluster
    if(fTPCData->fClusters.size() > 1)
        return;

    const auto& noise {fTPCData->fRaw};
    int iter {120};
    int minVoxels {7};
    double distThresh {2.};
    ActAlgorithm::RANSAC ransac {iter, minVoxels, distThresh};
    auto [clusters, back] {ransac.Run(noise)};
    if(fIsVerbose)
    {
        std::cout << BOLDGREEN << "-- RecRANSAC --" << '\n';
        std::cout << "Noise size: " << noise.size() << '\n';
        std::cout << "New clusters: " << clusters.size() << RESET << '\n';
    }
    if(clusters.size())
    {
        // Get clusters with best chi2
        std::sort(clusters.begin(), clusters.end(), [](const ActRoot::Cluster& a, const ActRoot::Cluster& b)
                  { return a.GetLine().GetChi2() < b.GetLine().GetChi2(); });
        // Push the first one
        fTPCData->fClusters.push_back(clusters.front());
    }
}

void ActAlgorithm::UserAction::Print() const
{
    std::cout << BOLDCYAN << "····· " << GetActionID() << " ·····" << '\n';
    if(!fIsEnabled)
    {
        std::cout << "······························" << RESET << '\n';
        return;
    }
    std::cout << "  MaxAngle       : " << fMaxAngle << '\n';
    std::cout << "  MinLength      : " << fMinLength << '\n';
    std::cout << "  CylinderRadius : " << fCylinderR << RESET << '\n';
}

// Create symbol to load class from .so
extern "C" ActAlgorithm::UserAction* CreateUserAction()
{
    return new ActAlgorithm::UserAction;
}
