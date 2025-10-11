#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActRANSAC.h"
#include "ActTPCData.h"

#include "UserAction.h"

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
    // Assign the highest angle cluster to the light
    ActRoot::Cluster clusterLight {};
    auto clusters {fTPCData->fClusters};
    double angleLight {-1};
    for(auto c : clusters)
    {
        double angleCluster {c.GetLine().GetDirection().Dot(ROOT::Math::XYZVectorF(1, 0, 0))};
        angleCluster = TMath::Abs(angleCluster);
        if(angleCluster > angleLight)
        {
            angleLight = angleCluster;
            clusterLight = c;
        }
    }

    // Now check if there is a shift in z coordinate respect the rp
    clusterLight.SortAlongDir();
    auto rp {fTPCData->fRPs[0]};

    double zFirstClusterLight {clusterLight.GetRefToVoxels().front().GetPosition().Z()};

    double distMinDecay {5.}; // pads
    if(std::abs(zFirstClusterLight - rp.Z()) > distMinDecay)
    {
        if(fIsVerbose)
        {
            std::cout << BOLDGREEN << "-- Decay filter --" << '\n';
            std::cout << "First Z cluster light: " << zFirstClusterLight << '\n';
            std::cout << "Z RP: " << rp.Z() << '\n';
            std::cout << "Delta Z: " << std::abs(zFirstClusterLight - rp.Z()) << RESET << '\n';
        }
        // Clear all clusters and push only the light one
        fTPCData->fClusters.clear();
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
