#ifndef Pipe2_Ex_cxx
#define Pipe2_Ex_cxx

#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

#include <iostream>
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
    // Build energy at vertex
    auto dfVertex = df.Define("EVertex",
                              [&](const ActRoot::MergerData& d)
                              {
                                  double ret {};
                                  if(d.fLight.IsFilled())
                                      ret = srim->EvalInitialEnergy(light, d.fLight.fEs.front(), d.fLight.fTL);
                                  else // L1 trigger
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
    double EBeamIni {4.6}; // AMeV at X = 0 of pad plane
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
                          {"Rec_EBeam"})
    };

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


    // Kinematics and Ex
    auto hKin {def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};
    auto hKinCM {def.Histo2D(HistConfig::KinCM, "ThetaCM", "EVertex")};
    auto hEBeam {def.Histo1D("EBeam")};
    auto hExSil {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
                     .Histo1D(HistConfig::Ex, "Ex")};
    hExSil->SetTitle("Ex with silicons");
    auto hExL1 {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                    .Histo1D(HistConfig::Ex, "Ex")};
    hExL1->SetTitle("Ex with L1");

    auto hThetaBeam {def.Histo2D(HistConfig::ThetaBeam, "fRP.fCoordinates.fX", "fThetaBeam")};
    auto hRP {def.Histo2D(HistConfig::RP, "fRP.fCoordinates.fX", "fRP.fCoordinates.fY")};
    auto hRPx {def.Histo1D(HistConfig::RPx, "fRP.fCoordinates.fX")};
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
    auto hRecECM {def.Histo1D(HistConfig::ECM, "Rec_ECM")};
    auto hECMRPx {def.Histo2D(HistConfig::RPxECM, "fRP.fCoordinates.fX", "ECM")};
    auto hRecECMRPx {def.Histo2D(HistConfig::RPxECM, "fRP.fCoordinates.fX", "Rec_ECM")};
    auto hEpRMg {def.Histo2D(HistConfig::EpRMg, "EVertex", "RangeHeavy")};

    // Save!
    auto outfile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    def.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << '\n';


    auto* c22 {new TCanvas("c22", "Pipe2 canvas 2")};
    c22->DivideSquare(4);
    c22->cd(1);
    hRP->DrawClone();
    c22->cd(2);
    hRPx->DrawClone();
    c22->cd(3);
    hThetaBeam->DrawClone("colz");
    c22->cd(4);
    hEBeam->DrawClone();

    auto* c21 {new TCanvas("c21", "Pipe2 canvas 1")};
    c21->DivideSquare(6);
    c21->cd(1);
    hKin->DrawClone("colz");
    auto* theo {kin.GetKinematicLine3()};
    theo->Draw("same");
    c21->cd(2);
    hExSil->DrawClone();
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
    c24->cd(2);
    hRecECM->DrawClone();
    c24->cd(3);
    hECMRPx->DrawClone("colz");
    c24->cd(4);
    hRecECMRPx->DrawClone("colz");
    c24->cd(5);
    hEpRMg->DrawClone("colz");
}
#endif
