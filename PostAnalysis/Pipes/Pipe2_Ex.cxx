#ifndef Pipe2_Ex_cxx
#define Pipe2_Ex_cxx

#include "ActCutsManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TAttLine.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TVirtualPad.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../HistConfig.h"

void Pipe2_Ex(const std::string& beam, const std::string& target, const std::string& light)
{
    // Read data
    auto filename {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"PID_Tree", filename};

    // Init SRIM
    auto* srim {new ActPhysics::SRIM};
    // Correct SRIM names
    std::string srimName {};
    if(light == "d")
        srimName = "2H";
    else if(light == "p")
        srimName = "1H";
    else if(light == "t")
        srimName = "3H";
    else if(light == "3He")
        srimName = "3He";
    else if(light == "4He")
        srimName = "4He";
    int pressure {800}; // 20Me beam
    if(beam == "20Ne")
        pressure = 950;
    srim->ReadTable(light,
                    TString::Format("../Calibrations/SRIM/%s_%dmbar_95-5.txt", srimName.c_str(), pressure).Data());
    srim->ReadTable(beam, TString::Format("../Calibrations/SRIM/%s_%dmbar_95-5.txt", beam.c_str(), pressure).Data());


    // Read cuts
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("ep_range", "./Cuts/ep_range_p_20Mg.root");
    cuts.ReadCut("debug", "./Cuts/debug_l1.root");

    // Build energy at vertex
    auto dfVertex = df.Define("EVertex",
                              [&](const ActRoot::MergerData& d)
                              {
                                  double ret {};
                                  if(d.fLight.IsFilled() && std::isfinite(d.fLight.fTL))
                                      ret = srim->EvalInitialEnergy(light, d.fLight.fEs.front(), d.fLight.fTL);
                                  else if(d.fLight.IsL1()) // L1 trigger
                                      ret = srim->EvalEnergy(light, d.fLight.fTL);
                                  return ret;
                              },
                              {"MergerData"});

    // Init particles
    ActPhysics::Particle pb {beam};
    auto mbeam {pb.GetMass()};
    ActPhysics::Particle pt {target};
    auto mtarget {pt.GetMass()};
    ActPhysics::Particle pl {light};
    // Declare kinematics
    // double EBeamIni {4.05235}; // AMeV at X = 0 of pad plane; energy meassure 5.44 before cfa
    double EBeamIni {4.17}; // AMeV at X = 0. WARNING: SUSPICTION THIS IS NOT 20Mg
    // That is, including elosses in CFA, entrance window, etc
    ActPhysics::Kinematics kin {pb, pt, pl, EBeamIni * pb.GetAMU()};
    // Vector of kinematics as one object is needed per
    // processing slot (since we are changing EBeam in each entry)
    std::vector<ActPhysics::Kinematics> vkins {df.GetNSlots()};
    for(auto& k : vkins)
        k = kin;

    // Beam energy calculation and ECM
    auto def {dfVertex
                  .Define("EBeam", [&](const ActRoot::MergerData& d)
                          { return srim->Slow(beam, EBeamIni * pb.GetAMU(), d.fRP.X()); }, {"MergerData"})
                  .DefineSlot("Rec_EBeam", // assuming Ex = 0 using outgoing light particle kinematics
                              [&](unsigned int slot, double EVertex, const ActRoot::MergerData& d)
                              {
                                  // no need for slots here but for the sake of consistency with next calculations...
                                  return vkins[slot].ReconstructBeamEnergyFromLabKinematics(
                                      EVertex, d.fThetaLight * TMath::DegToRad());
                              },
                              {"EVertex", "MergerData"})
                  .Define("ECM", [&](double EBeam) { return (mtarget / (mbeam + mtarget)) * EBeam; }, {"EBeam"})
                  .Define("Rec_ECM", [&](double rec_EBeam) { return (mtarget / (mbeam + mtarget)) * rec_EBeam; },
                          {"Rec_EBeam"})};

    def =
        def.DefineSlot("Ex",
                       [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                       {
                           vkins[slot].SetBeamEnergy(EBeam);
                           return vkins[slot].ReconstructExcitationEnergy(EVertex, (d.fThetaLight) * TMath::DegToRad());
                       },
                       {"MergerData", "EVertex", "EBeam"})
            .DefineSlot("ThetaCM",
                        [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                        {
                            vkins[slot].SetBeamEnergy(EBeam);
                            return vkins[slot].ReconstructTheta3CMFromLab(EVertex,
                                                                          (d.fThetaLight) * TMath::DegToRad()) *
                                   TMath::RadToDeg();
                        },
                        {"MergerData", "EVertex", "EBeam"});

    // Define range of heavy particle
    def = def.Define("RangeHeavy", [&](ActRoot::MergerData& d) { return d.fHeavy.fTL; }, {"MergerData"});


    // Create two nodes
    auto nodel0 {def.Filter([](ActRoot::MergerData& d) { return d.fLight.IsL1() == false; }, {"MergerData"})};
    auto nodeLat {nodel0.Filter([](ActRoot::MergerData& d)
                                { return d.fLight.fLayers.front() == "l0" || d.fLight.fLayers.front() == "r0"; },
                                {"MergerData"})};
    auto nodeFront {
        nodel0.Filter([](ActRoot::MergerData& d) { return d.fLight.fLayers.front() == "f0"; }, {"MergerData"})};
    auto nodel1 {def.Filter([](ActRoot::MergerData& d) { return d.fLight.IsL1() == true; }, {"MergerData"})};

    // Filter in Ep vs Range
    auto nodeEpRange {nodel0.Filter([&](float range, double elab) { return cuts.IsInside("ep_range", range, elab); },
                                    {"RangeHeavy", "EVertex"})};

    // Kinematics and Ex
    auto hKin {def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};
    auto hKinCM {def.Histo2D(HistConfig::KinCM, "ThetaCM", "EVertex")};
    auto hEBeam {def.Histo1D(HistConfig::EBeam, "EBeam")};
    auto hRecEBeam {def.Histo1D(HistConfig::EBeam, "Rec_EBeam")};
    auto hExSilLat {nodeLat.Histo1D(HistConfig::Ex, "Ex")};
    hExSilLat->SetTitle("Side");
    auto hExSilFront {nodeFront.Histo1D(HistConfig::Ex, "Ex")};
    hExSilFront->SetTitle("Front");
    auto hExSil {nodel0.Histo1D(HistConfig::Ex, "Ex")};
    hExSil->SetTitle("All");

    auto hExL1 {nodel1.Histo1D(HistConfig::Ex, "Ex")};
    hExL1->SetTitle("Ex with L1");

    auto hThetaBeam {def.Histo2D(HistConfig::ThetaBeam, "fRP.fCoordinates.fX", "fThetaBeam")};
    auto hRP {def.Histo2D(HistConfig::RP, "fRP.fCoordinates.fX", "fRP.fCoordinates.fY")};
    auto hRPx {def.Histo1D(HistConfig::RPx, "fRP.fCoordinates.fX")};
    auto hRPxFront {nodeFront.Histo1D(HistConfig::RPx, "fRP.fCoordinates.fX")};
    hRPxFront->SetLineColor(46);
    hRPxFront->SetTitle("Front");
    auto hRPxSide {nodeLat.Histo1D(HistConfig::RPx, "fRP.fCoordinates.fX")};
    hRPxSide->SetLineColor(8);
    hRPxSide->SetTitle("Side");
    auto hRPxL1 {nodel1.Histo1D(HistConfig::RPx, "fRP.fCoordinates.fX")};
    hRPxL1->SetLineColor(kOrange);
    hRPxL1->SetTitle("L1");

    auto hThetaCMLab {def.Histo2D(HistConfig::ThetaCMLab, "fThetaLight", "ThetaCM")};

    // Ex dependences
    auto hExThetaCM {def.Histo2D(HistConfig::ExThetaCM, "ThetaCM", "Ex")};
    auto hExThetaLab {def.Histo2D(HistConfig::ExThetaLab, "fThetaLight", "Ex")};
    auto hExRP {def.Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "Ex")};

    // Heavy histograms
    auto hThetaHLLab {def.Histo2D(HistConfig::ChangeTitle(HistConfig::ThetaHeavyLight, "Lab correlations"),
                                  "fThetaLight", "fThetaHeavy")};

    // CM things
    auto hECM {def.Histo1D(HistConfig::ECM, "ECM")};
    hECM->SetTitle("All");
    auto hECMLat {nodeLat.Histo1D(HistConfig::ECM, "ECM")};
    hECMLat->SetLineColor(8);
    hECMLat->SetTitle("Side");
    auto hECMFront {nodeFront.Histo1D(HistConfig::ECM, "ECM")};
    hECMFront->SetLineColor(46);
    hECMFront->SetTitle("Front");
    auto hECML1 {nodel1.Histo1D(HistConfig::ECM, "ECM")};
    hECML1->SetLineColor(kOrange);
    hECML1->SetTitle("L1");

    auto hRecECM {def.Histo1D(HistConfig::ECM, "Rec_ECM")};
    hRecECM->SetTitle("Rec E_{CM} with E_{x} = 0");

    auto hECM2dFront {nodeFront.Histo2D(HistConfig::ThetaCMECM, "ThetaCM", "ECM")};
    hECM2dFront->SetTitle("Front");
    auto hECM2dLat {nodeLat.Histo2D(HistConfig::ThetaCMECM, "ThetaCM", "ECM")};
    hECM2dLat->SetTitle("Side");
    auto hECM2dL1 {nodel1.Histo2D(HistConfig::ThetaCMECM, "ThetaCM", "ECM")};
    hECM2dL1->SetTitle("L1");
    
    auto hRPxELab {
        nodel1.Filter("70 < fThetaLight && fThetaLight < 80")
        .Histo2D({"hRPxELab", "#theta_{Lab} in [70, 80];RP.X [mm];E_{Vertex} [#circ]", 400, 0, 260, 300, 0, 30},
                       "fRP.fCoordinates.fX", "EVertex")};

    auto hRPxThetaLab {
        nodel1.Histo2D({"hRPxThetaLab", "L1 exlusion zone;RP.X [mm];#theta_{Lab} [#circ]", 400, 0, 260, 300, 0, 90},
                       "fRP.fCoordinates.fX", "fThetaLight")};

    auto hECMRPx {nodeFront.Histo2D(HistConfig::RPxECM, "fRP.fCoordinates.fX", "ECM")};
    auto hRecECMRPx {def.Histo2D(HistConfig::RPxECM, "fRP.fCoordinates.fX", "Rec_ECM")};
    auto hEpRMg {def.Histo2D(HistConfig::EpRMg, "RangeHeavy", "EVertex")};

    auto hECMCut {nodeEpRange.Histo1D(HistConfig::ECM, "ECM")};

    // Save!
    auto outfile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    def.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << '\n';

    // std::ofstream streamer {"./debug_l1.dat"};
    // auto nodeStreamer {nodel1.Filter([&](double ECM, double thetaCM) { return cuts.IsInside("debug", thetaCM, ECM);
    // },
    //                                  {"ECM", "ThetaCM"})};
    // nodeStreamer.Foreach([&](ActRoot::MergerData& d) { d.Stream(streamer); }, {"MergerData"});
    // streamer.close();


    auto* c22 {new TCanvas("c22", "Pipe2 canvas 2")};
    c22->DivideSquare(4);
    c22->cd(1);
    hRP->DrawClone();
    c22->cd(2);
    hRPx->DrawClone();
    hRPxFront->DrawClone("same");
    hRPxSide->DrawClone("same");
    hRPxL1->DrawClone("same");
    gPad->BuildLegend();
    c22->cd(3);
    hThetaBeam->DrawClone("colz");
    c22->cd(4);
    hEBeam->DrawClone();
    hRecEBeam->SetLineColor(8);
    hRecEBeam->DrawClone("same");

    auto* c21 {new TCanvas("c21", "Pipe2 canvas 1")};
    c21->DivideSquare(6);
    c21->cd(1);
    hKin->DrawClone("colz");
    auto* theo {kin.GetKinematicLine3()};
    kin.SetEx(1.6);
    auto* theoIn {kin.GetKinematicLine3()};
    theoIn->SetLineColor(46);
    theoIn->SetLineStyle(kDashed);
    theo->Draw("same");
    theoIn->Draw("same");
    c21->cd(2);
    hExSil->DrawClone();
    hExSilLat->SetLineColor(8);
    hExSilLat->DrawClone("same");
    hExSilFront->SetLineColor(46);
    hExSilFront->DrawClone("same");
    gPad->BuildLegend();
    c21->cd(3);
    hKinCM->DrawClone("colz");
    c21->cd(4);
    hExL1->DrawClone();
    // hExThetaLab->DrawClone("colz");
    c21->cd(5);
    hExThetaCM->DrawClone("colz");
    c21->cd(6);
    hExRP->DrawClone("colz");

    auto* c23 {new TCanvas {"c23", "Pipe2 canvas 3"}};
    c23->DivideSquare(4);
    c23->cd(1);
    hThetaHLLab->DrawClone("colz");
    c23->cd(2);
    hThetaCMLab->DrawClone("colz");
    c23->cd(3);

    auto* c24 {new TCanvas {"c24", "Pipe2 canvas 4"}};
    c24->DivideSquare(6);
    c24->cd(1);
    hECM->DrawClone();
    hECMLat->DrawClone("same");
    hECMFront->DrawClone("same");
    hECML1->DrawClone("same");
    gPad->BuildLegend();
    c24->cd(2);
    // hRPxELab->DrawClone("colz");
    hECMRPx->DrawClone("colz");
    // hECM2dFront->DrawClone("colz");
    // hECM2dLat->DrawClone("colz same");
    c24->cd(3);
    hRPxThetaLab->DrawClone("colz");
    c24->cd(4);
    hECM2dL1->DrawClone("colz");
    c24->cd(5);
    hEpRMg->DrawClone("colz");
    cuts.DrawCut("ep_range");
    c24->cd(6);
    hECMCut->DrawClone("colz");
    // hECMRPx->DrawClone("colz");
}
#endif
