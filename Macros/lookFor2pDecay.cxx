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

void lookFor2pDecay()
{
    std::string dataconf {"./../configs/data.conf"};

    // RDataFrame
    ROOT::EnableImplicitMT();

    // std::ofstream streamer {"./Outputs/debug_3p_decay.dat"};
    // df_filtered.Foreach([&](ActRoot::MergerData& d) { d.Stream(streamer); }, {"MergerData"});
    // streamer.close();

    // srim files
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("20Na", "../Calibrations/SRIM/20Na_950mbar_95-5.txt");

    // Init particles
    ActPhysics::Particle pb {"20Na"};
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
    
    // Get data from file
    std::string outfile {"./Outputs/tree_ex_20Na_p_2p.root"};
    ROOT::RDataFrame df_angles {"Final_Tree", outfile};

    // Draw them
    auto* c1 {new TCanvas("c1", "2p decay ECM")};
    c1->DivideSquare(2);
    c1->cd(1);
    auto hECM {df_angles.Histo1D(HistConfig::ECM, "ECM")};
    hECM->DrawClone();
    c1->cd(2);
    auto hRpx {df_angles.Histo1D(HistConfig::RPx, "RPx")};
    hRpx->DrawClone();
    // Plot angles
    // auto* c2 {new TCanvas("c2", "2p decay angles")};
    // c2->DivideSquare(3);
    // c2->cd(1);
    // const ROOT::RDF::TH1DModel ThetaLab {"hThetaLab", "ThetaLab;#theta_{lab} [#circ]", 600, 0, 180};
    // auto hTheta1 {df_angles.Histo1D(ThetaLab, "threeAngles.theta1")};
    // hTheta1->DrawClone();
    // c2->cd(2);
    // auto hTheta2 {df_angles.Histo1D(ThetaLab, "threeAngles.theta2")};
    // hTheta2->DrawClone();
    // c2->cd(3);
    // auto hTheta3 {df_angles.Histo1D(ThetaLab, "threeAngles.theta3")};

    // Plot the Dalizt plots
    // auto cDalitz {new TCanvas("cDalitz", "3p decay Dalitz")};
    // const ROOT::RDF::TH2DModel Dalitz {"hDalitz", "Dalitz;#epsilon_{1};#epsilon_{2}", 200, 0, 1, 200, 0, 1};
    // cDalitz->DivideSquare(3);
    // cDalitz->cd(1);
    // auto hDalitz {df_angles.Histo2D(Dalitz, "epsilon1", "epsilon2")};
    // hDalitz->DrawClone("colz");
    // cDalitz->cd(2);
    // auto hDalitz2 {df_angles.Histo2D(Dalitz, "epsilon2", "epsilon3")};
    // hDalitz2->GetXaxis()->SetTitle("#epsilon_{2}");
    // hDalitz2->GetYaxis()->SetTitle("#epsilon_{3}");
    // hDalitz2->DrawClone("colz");
    // cDalitz->cd(3);
    // auto hDalitz3 {df_angles.Histo2D(Dalitz, "epsilon3", "epsilon1")};
    // hDalitz3->GetXaxis()->SetTitle("#epsilon_{3}");
    // hDalitz3->GetYaxis()->SetTitle("#epsilon_{1}");
    // hDalitz3->DrawClone("colz");
    
}