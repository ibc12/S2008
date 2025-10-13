#ifndef Pipe2_Ex_cxx
#define Pipe2_Ex_cxx

#include "ActCutsManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
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
#include <stdexcept>
#include <string>
#include <vector>

#include "../HistConfig.h"

void Pipe2_Ex(const std::string& beam, const std::string& target, const std::string& light)
{
    // Read data
    auto filename {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    ROOT::EnableImplicitMT();
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
    // Beam energies
    // First set of data, runs [31-40]
    double EBeamFirst {4.24}; // MeV/u from SRIM interpolation of exp. 20Mg range
    // 2nd set of data, runs [47, onwards] -> is the SAME
    double EBeamSecond {4.24}; // MeV/u from SRIM
    std::map<int, double> EBeams;
    for(int run = 31; run <= 40; run++)
        EBeams[run] = EBeamFirst;
    for(int run = 47; run <= 60; run++)
        EBeams[run] = EBeamSecond;
    // Allegedly wrong LISE calculation below
    // double EBeamIni {4.05235}; // AMeV at X = 0 of pad plane; energy meassure 5.44 before cfa
    // Energy of unidentified beam on Sunday 12
    // double EBeamIni {4.17}; // AMeV at X = 0. WARNING: SUSPICTION THIS IS NOT 20Mg

    // All the above beam energies include energy losses in CFA, window, etc...
    // They're given at X = 0 of the pad plane
    ActPhysics::Kinematics kin {pb, pt, pl, EBeamFirst * pb.GetAMU()};
    // Vector of kinematics as one object is needed per
    // processing slot (since we are changing EBeam in each entry)
    std::vector<ActPhysics::Kinematics> vkins {df.GetNSlots()};
    for(auto& k : vkins)
        k = kin;

    // Beam energy calculation and ECM
    auto def {dfVertex
                  .Define("EBeam",
                          [&](const ActRoot::MergerData& d)
                          {
                              double EBeam {};
                              if(EBeams.count(d.fRun))
                                  EBeam = EBeams[d.fRun];
                              else
                                  throw std::runtime_error("Defining EBeam: no initial beam energy for run " +
                                                           std::to_string(d.fRun));
                              return srim->Slow(beam, EBeam * pb.GetAMU(), d.fRP.X());
                          },
                          {"MergerData"})
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


    // Create node to gate on different conditions: silicon layer, l1, etc
    // L0 trigger
    auto nodel0 {def.Filter([](ActRoot::MergerData& d) { return d.fLight.IsL1() == false; }, {"MergerData"})};
    // L0 -> side silicons
    auto nodeLat {nodel0.Filter([](ActRoot::MergerData& d)
                                { return d.fLight.fLayers.front() == "l0" || d.fLight.fLayers.front() == "r0"; },
                                {"MergerData"})};
    // L0 -> front silicons
    auto nodeFront {
        nodel0.Filter([](ActRoot::MergerData& d) { return d.fLight.fLayers.front() == "f0"; }, {"MergerData"})};

    // L1 trigger
    auto nodel1 {def.Filter([](ActRoot::MergerData& d) { return d.fLight.IsL1() == true; }, {"MergerData"})};

    // Selection in Ep vs R20Mg plot
    auto nodeEpRange {nodel0.Filter([&](float range, double elab) { return cuts.IsInside("ep_range", range, elab); },
                                    {"RangeHeavy", "EVertex"})};

    // Kinematics and Ex
    auto hKin {def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};
    auto hKinCM {def.Histo2D(HistConfig::KinCM, "ThetaCM", "EVertex")};

    // Create vector of nodes and labels
    std::vector<std::string> labels {"All", "Lat", "Front", "L1"};
    std::vector<ROOT::RDF::RNode> nodes {def, nodeLat, nodeFront, nodel1};

    // EBeam
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsEBeam;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::EBeam, "EBeam")};
        hsEBeam.push_back(h);
    }
    // Rec EBeam
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsRecEBeam;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::EBeam, "Rec_EBeam")};
        hsRecEBeam.push_back(h);
    }

    // Ex
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsEx;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::Ex, "Ex")};
        hsEx.push_back(h);
    }
    // ECM
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsECM;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::ECM, "ECM")};
        hsECM.push_back(h);
    }
    std::vector<ROOT::RDF::RResultPtr<TH2D>> hsECM2d;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo2D(HistConfig::ThetaCMECM, "ThetaCM", "ECM")};
        hsECM2d.push_back(h);
    }
    // RPx
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hsRPx;
    for(int i = 0; i < labels.size(); i++)
    {
        auto h {nodes[i].Histo1D(HistConfig::RPx, "fRP.fCoordinates.fX")};
        hsRPx.push_back(h);
    }
    auto hThetaBeam {def.Histo2D(HistConfig::ThetaBeam, "fRP.fCoordinates.fX", "fThetaBeam")};
    auto hRP {def.Histo2D(HistConfig::RP, "fRP.fCoordinates.fX", "fRP.fCoordinates.fY")};
    auto hThetaCMLab {def.Histo2D(HistConfig::ThetaCMLab, "fThetaLight", "ThetaCM")};
    // Ex dependences
    auto hExThetaCM {def.Histo2D(HistConfig::ExThetaCM, "ThetaCM", "Ex")};
    auto hExThetaLab {def.Histo2D(HistConfig::ExThetaLab, "fThetaLight", "Ex")};
    auto hExRP {def.Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "Ex")};
    auto hExZ {nodel0.Histo2D(HistConfig::ExZ, "fSP.fCoordinates.fZ", "Ex")};
    // Heavy histograms
    auto hThetaHLLab {def.Histo2D(HistConfig::ChangeTitle(HistConfig::ThetaHeavyLight, "Lab correlations"),
                                  "fThetaLight", "fThetaHeavy")};
    auto hRecECM {def.Histo1D(HistConfig::ECM, "Rec_ECM")};
    hRecECM->SetTitle("Rec E_{CM} with E_{x} = 0");

    // Histograms for online analysis
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


    // Set styles for histograms
    std::vector<int> colors {-1, 8, 46, 9};
    for(int i = 0; i < hsEBeam.size(); i++)
    {
        auto& col {colors[i]};
        auto title {labels[i].c_str()};
        hsEBeam[i]->SetTitle(title);
        hsEBeam[i]->SetLineColor(col);

        hsRecEBeam[i]->SetTitle(title);
        hsRecEBeam[i]->SetLineColor(col);

        hsEx[i]->SetTitle(title);
        hsEx[i]->SetLineColor(col);

        hsECM[i]->SetTitle(title);
        hsECM[i]->SetLineColor(col);

        hsECM2d[i]->SetTitle(title);

        hsRPx[i]->SetTitle(title);
        hsRPx[i]->SetLineColor(col);
        if(i == 0)
        {
            hsEBeam[i]->SetTitle("E_{beam}");
            hsRecEBeam[i]->SetTitle("Rec E_{beam}");
            hsEx[i]->SetTitle("E_{x}");
            hsECM[i]->SetTitle("E_{CM}");
            hsECM2d[i]->SetTitle("E_{CM} vs #theta_{CM}");
            hsRPx[i]->SetTitle("RP_{x}");
        }
    }


    auto* c22 {new TCanvas("c22", "Pipe2 canvas 2")};
    c22->DivideSquare(4);
    c22->cd(1);
    hRP->DrawClone();
    c22->cd(2);
    for(int i = 0; i < hsRPx.size(); i++)
        hsRPx[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();
    c22->cd(3);
    for(int i = 0; i < hsRPx.size(); i++)
        hsEBeam[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();
    c22->cd(4);
    for(int i = 0; i < hsRPx.size(); i++)
        hsRecEBeam[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();

    auto* c21 {new TCanvas("c21", "Pipe2 canvas 1")};
    c21->DivideSquare(6);
    c21->cd(1);
    hKin->DrawClone("colz");
    auto* theo {kin.GetKinematicLine3()};
    kin.SetEx(1.6);
    auto* theoIne {kin.GetKinematicLine3()};
    theoIne->SetLineColor(46);
    theoIne->SetLineStyle(kDashed);
    theo->Draw("same");
    theoIne->Draw("same");
    c21->cd(2);
    for(int i = 0; i < hsRPx.size() - 1; i++)
        hsEx[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();
    c21->cd(3);
    hsEx.back()->SetTitle("E_{x} with L1");
    hsEx.back()->DrawClone();
    c21->cd(4);
    hExThetaLab->DrawClone("colz");
    c21->cd(5);
    hExZ->DrawClone("colz");
    c21->cd(6);
    hExRP->DrawClone("colz");

    // auto* c23 {new TCanvas {"c23", "Pipe2 canvas 3"}};
    // c23->DivideSquare(4);
    // c23->cd(1);
    // hThetaHLLab->DrawClone("colz");
    // c23->cd(2);
    // hThetaCMLab->DrawClone("colz");
    // c23->cd(3);

    auto* c24 {new TCanvas {"c24", "Pipe2 canvas 4"}};
    c24->DivideSquare(6);
    c24->cd(1);
    for(int i = 0; i < hsRPx.size(); i++)
        hsECM[i]->DrawClone(i == 0 ? "" : "same");
    gPad->BuildLegend();
    c24->cd(2);
    hECMRPx->DrawClone("colz");
    c24->cd(3);
    hRPxThetaLab->DrawClone("colz");
    c24->cd(4);
    hsECM2d.front()->DrawClone("colz");
    // hECM2dL1->DrawClone("colz");
    c24->cd(5);
    hEpRMg->DrawClone("colz");
    cuts.DrawCut("ep_range");
    c24->cd(6);
    hECMCut->SetTitle("E_{CM} with cut on elastic");
    hECMCut->DrawClone("colz");
    // hECMRPx->DrawClone("colz");
}
#endif
