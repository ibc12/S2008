#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "../PostAnalysis/HistConfig.h"

struct twoAngles
{
    float theta1;
    float theta2;
};

void lookFor3pDecay()
{
    std::string dataconf {"./../configs/data.conf"};

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    chain->AddFriend(chain3.get());

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto df_filtered {df.Filter([](ActRoot::TPCData& d) { return d.fClusters.size() == 5; }, {"TPCData"})};

    std::ofstream streamer {"./Outputs/debug_2p_decay.dat"};
    df_filtered.Foreach([&](ActRoot::MergerData& d) { d.Stream(streamer); }, {"MergerData"});
    streamer.close();

    // srim files
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("20Mg", "../Calibrations/SRIM/20Mg_800mbar_95-5.txt");

    // Init particles
    ActPhysics::Particle pb {"20Mg"};
    auto mbeam {pb.GetMass()};
    ActPhysics::Particle pt {"p"};
    auto mtarget {pt.GetMass()};
    ActPhysics::Particle pl {"p"};

    double EBeam {4.24}; // MeV/u from SRIM interpolation of exp. 20Mg range
    // All the above beam energies include energy losses in CFA, window, etc...
    // They're given at X = 0 of the pad plane
    ActPhysics::Kinematics kin {pb, pt, pl, EBeam * pb.GetAMU()};
    // Vector of kinematics as one object is needed per
    // processing slot (since we are changing EBeam in each entry)
    std::vector<ActPhysics::Kinematics> vkins {df.GetNSlots()};
    for(auto& k : vkins)
        k = kin;

    // Output 

    // Beam energy calculation and ECM
    auto def {df_filtered
                  .Define("EBeam", [&](const ActRoot::TPCData& d)
                          { return srim->Slow("20Mg", EBeam * pb.GetAMU(), d.fRPs[0].X()); }, {"TPCData"})
                  .Define("ECM", [&](double EBeam) { return (mtarget / (mbeam + mtarget)) * EBeam; }, {"EBeam"})
                  .Define("RPx", [&](ActRoot::TPCData& d) { return d.fRPs[0].X(); }, {"TPCData"})};
    // Create node to gate on different conditions: silicon layer, l1, etc
    //// L0 trigger
    // auto nodel0 {def.Filter([](ActRoot::MergerData& d) { return d.fLight.IsL1() == false; }, {"MergerData"})};
    //// L0 -> side silicons
    // auto nodeLat {nodel0.Filter([](ActRoot::MergerData& d)
    //                             { return d.fLight.fLayers.front() == "l0" || d.fLight.fLayers.front() == "r0"; },
    //                             {"MergerData"})};
    //// L0 -> front silicons
    // auto nodeFront {
    //     nodel0.Filter([](ActRoot::MergerData& d) { return d.fLight.fLayers.front() == "f0"; }, {"MergerData"})};
    //
    //// L1 trigger
    // auto nodel1 {def.Filter([](ActRoot::MergerData& d) { return d.fLight.IsL1() == true; }, {"MergerData"})};

    // // Plot Ecm
    // std::vector<std::string> labels {"All", "Lat", "Front", "L1"};
    // std::vector<ROOT::RDF::RNode> nodes {def, nodeLat, nodeFront, nodel1};
    // // ECM histos
    // std::vector<ROOT::RDF::RResultPtr<TH1D>> hsECM;
    // for(int i = 0; i < labels.size(); i++)
    // {
    //     auto h {nodes[i].Histo1D(HistConfig::ECM, "ECM")};
    //     hsECM.push_back(h);
    // }

    // Get Angles of light and heavy
    auto df_angles {def.Define(
        "threeAngles",
        [&](ActRoot::TPCData& d)
        {
            twoAngles angles {};
            auto clusters {d.fClusters};
            int beamIdx {};
            int idx {};
            for(auto c : clusters)
            {
                auto isBeam {c.GetIsBeamLike()};
                if(isBeam)
                {
                    beamIdx = idx;
                    break;
                }
                idx++;
            }
            // Loop over rest of clusters to get heavy idx
            int heavyIdx {};
            idx = 0;
            float angle {1000.};
            for(auto c : clusters)
            {
                if(idx != beamIdx)
                {
                    auto line {c.GetLine()};
                    auto anglePreliminar {line.GetDirection().Unit().Dot(ROOT::Math::XYZVector(1, 0, 0))};
                    if(anglePreliminar < angle)
                    {
                        heavyIdx = idx;
                        angle = anglePreliminar;
                        break;
                    }
                }
                idx++;
            }
            // Get angles of light respect to heavy
            idx = 0;
            for(auto c : clusters)
            {
                if(idx != beamIdx && idx != heavyIdx)
                {
                    auto line {c.GetLine()};
                    auto direction {line.GetDirection().Unit()};
                    auto theta {std::acos(direction.Dot(clusters[heavyIdx].GetLine().GetDirection().Unit())) * 180. /
                                TMath::Pi()};
                    if(angles.theta1 == 0.)
                        angles.theta1 = theta;
                    else if(angles.theta2 == 0.)
                        angles.theta2 = theta;
                }
                idx++;
            }

            return angles;
        },
        {"TPCData"})};

    std::string outfile {"./Outputs/tree_ex_20Mg_p_3p.root"};
    // Save the tree with angles and epsilon
    df_angles.Snapshot("Final_Tree", outfile);

}